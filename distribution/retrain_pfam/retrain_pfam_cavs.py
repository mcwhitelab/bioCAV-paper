#!/usr/bin/env python
"""
Retrain the 27,472 Pfam CAVs in the shared reference scaler frame.

Why: the original library (Ahmad, /xdisk/.../shamail/tcav_outputs_esmplusplus_all)
fit a fresh StandardScaler per motif on (positives + 10k negatives).  The paper
methods describe a single shared scaler, and GO/EC now use one
(shared_protein/reference_population/scaler_v1.pkl, fit on the same 10k
negatives).  Retraining in that frame makes the Pfam CAVs directly comparable
to GO/EC and usable without shipping 27k per-motif scalers.

Everything except the scaler is Ahmad's code path, unchanged:
esm2_tcav.src.train_cav.train_layer_cavs, with each motif's original report
`config` (C=1, balanced, lbfgs, 5-fold CV, seed 42, no PCA).  Two patches:

  1. create_scaler_and_pca  -> returns the shared scaler (no fit).
  2. load_embeddings_for_layer -> drops positive rows whose window starts at
     or past residue 1998.  The embedder truncates at 2,000 tokens (BOS +
     residues 0..1997 + EOS) but collection allowed proteins < 2,048 aa, so
     those windows pool nothing (NaN; killed 13 families in the original run)
     or only the EOS token (PF06943 row 11, trained silently).  Also drops any
     remaining non-finite row.  Count recorded as `n_dropped_truncated_pos`.
     Windows that merely cross the cap are kept as-is (2,789 rows, 423
     families), so this run differs from the original only in the scaler.

Embeddings are read in place; nothing is copied.  Idempotent: a motif whose
report in OUT already carries the `retrain` block is skipped.

Usage: retrain_pfam_cavs.py --shard I --n-shards N [--out DIR] [--motifs ...]
"""
import argparse, json, logging, sys, time, traceback
from pathlib import Path

import joblib
import numpy as np

sys.path.insert(0, "/groups/clairemcwhite/ahmad_workspace/esm_c")
from esm2_tcav.src import train_cav as tc  # noqa: E402

SRC = Path("/xdisk/clairemcwhite/shamail/tcav_outputs_esmplusplus_all")
EMB = SRC / "embeddings" / "ESMplusplus_large"
SHARED_SCALER = Path("/groups/clairemcwhite/cav_workspace/bioCAV-paper/"
                     "shared_protein/reference_population/scaler_v1.pkl")
DEFAULT_OUT = Path("/xdisk/clairemcwhite/clairemcwhite/pfam_cavs_shared_scaler")
LAYER, VERSION = 25, "v1"
TOKEN_CAP_RESIDUES = 1998  # max_length=2000 tokens minus BOS/EOS

logging.basicConfig(level=logging.WARNING,
                    format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger("retrain_pfam")
log.setLevel(logging.INFO)

SCALER = joblib.load(SHARED_SCALER)
assert SCALER.n_samples_seen_ == 10000, "expected the negatives-only shared scaler"
_dropped = {}


def shared_scaler_and_no_pca(X, use_pca=False, pca_dim=None):
    assert not use_pca, "original Pfam config has use_pca=False"
    assert X.shape[1] == SCALER.mean_.shape[0]
    return SCALER, None, {
        "scaler_mean": SCALER.mean_.tolist(),
        "scaler_std": SCALER.scale_.tolist(),
        "input_dim": X.shape[1],
        "pca_applied": False,
        "output_dim": X.shape[1],
        "scaler_source": str(SHARED_SCALER),
    }


def load_finite(layer, embed_dir, pos_file=None, neg_file=None):
    pos = np.load(Path(embed_dir) / f"L{layer}_pos.npy")
    neg = np.load(Path(embed_dir) / f"L{layer}_neg.npy")
    assert np.isfinite(neg).all()
    meta = json.load(open(Path(embed_dir) / "pos" / f"L{layer}_meta.json"))
    assert len(meta) == len(pos)
    keep = np.isfinite(pos).all(axis=1) & np.array([m["window_span"][0] < TOKEN_CAP_RESIDUES for m in meta])
    _dropped["n"] = int((~keep).sum())
    pos = pos[keep]
    X = np.vstack([pos, neg])
    y = np.hstack([np.ones(len(pos)), np.zeros(len(neg))])
    return X, y, np.arange(len(X))


tc.create_scaler_and_pca = shared_scaler_and_no_pca
tc.load_embeddings_for_layer = load_finite


def original_config(motif):
    rep = SRC / "cavs" / motif / f"L{LAYER}_report_{VERSION}.json"
    if rep.exists():
        return json.load(open(rep))["config"]
    # the 13 families that failed originally have no report; use a sibling's
    return json.load(open(SRC / "cavs" / "PF00001" / f"L{LAYER}_report_{VERSION}.json"))["config"]


def train_one(motif, out):
    cav_dir = out / "cavs" / motif
    report = cav_dir / f"L{LAYER}_report_{VERSION}.json"
    if report.exists():
        try:
            if "retrain" in json.load(open(report)):
                return "skip"
        except json.JSONDecodeError:
            pass
    cav_dir.mkdir(parents=True, exist_ok=True)
    _dropped.clear()
    tc.train_layer_cavs(layer=LAYER, embed_dir=str(EMB / motif), output_dir=str(cav_dir),
                        config=original_config(motif), artifact_version=VERSION)
    r = json.load(open(report))
    r["concept_metrics"]["n_dropped_truncated_pos"] = _dropped.get("n", 0)
    r["retrain"] = {"scaler": str(SHARED_SCALER), "embeddings": str(EMB / motif),
                    "date": time.strftime("%Y-%m-%d")}
    json.dump(r, open(report, "w"), indent=2)
    return "ok"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--n-shards", type=int, default=1)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--motifs", nargs="*")
    a = ap.parse_args()

    motifs = a.motifs or sorted(p.name for p in EMB.iterdir() if p.name.startswith("PF"))
    motifs = motifs[a.shard::a.n_shards]
    log.info(f"shard {a.shard}/{a.n_shards}: {len(motifs)} motifs -> {a.out}")
    fails = a.out / "failed" / f"shard_{a.shard}.tsv"
    fails.parent.mkdir(parents=True, exist_ok=True)
    counts = {"ok": 0, "skip": 0, "fail": 0}
    t0 = time.time()
    for i, m in enumerate(motifs, 1):
        try:
            counts[train_one(m, a.out)] += 1
        except Exception as e:
            counts["fail"] += 1
            with open(fails, "a") as fh:
                fh.write(f"{m}\t{type(e).__name__}: {str(e).splitlines()[0] if str(e) else ''}\n")
            log.error(f"{m} failed:\n{traceback.format_exc()}")
        if i % 50 == 0:
            log.info(f"{i}/{len(motifs)} {counts} {(time.time() - t0) / i:.1f}s/motif")
    log.info(f"done {counts}")


if __name__ == "__main__":
    main()
