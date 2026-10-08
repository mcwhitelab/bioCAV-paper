#!/usr/bin/env python
"""
GO/EC retrain at C=1e-4 (2026-09-29), step 2: train every GO and EC CAV with the Pfam
decision (METHODS §1.6): L2 logistic regression, C=1e-4, class_weight='balanced',
lbfgs tol=1e-10, max_iter=100000, random_state=42, shared reference scaler
(scaler_v1.pkl, 10k negatives only), negatives = the shared neg_10000.pkl, coef_ unit-normalised.

Positives come from step 1's pooled unique sequences: for each term, every span accession
is looked up in THAT term's .span.fasta (sequences differ between sets for 131
accessions), hashed, and mapped to its pooled row; one row per span line, in span order,
as prepare_embeddings.py does.  Sequences the embedder dropped (> 2048 aa) and span IDs
with no fasta record are skipped, as in the per-term runs; n_positive is checked against
the shipped report.

Output mirrors the shipped layout under a new root (nothing shipped is touched):
  /xdisk/clairemcwhite/clairemcwhite/go_ec_cavs_C1e-4_sharedscaler/<GO_sets|EC_sets>/<part>/<term>/
      random_positive_train_max1000_cav/{concept_v1.npy, report_v1.json}
Usage: 02_train_go_ec_C1e-4.py --shard I --n-shards N   (idempotent)
"""
import argparse, glob, hashlib, json, sys, time, traceback
from pathlib import Path

import joblib
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

PAPER = Path("/groups/clairemcwhite/cav_workspace/bioCAV-paper")
sys.path.insert(0, str(PAPER.parent / "bioCAV" / "core"))
from src.utils.data_loader import load_sequence_embeddings  # noqa: E402

REF = PAPER / "shared_protein" / "reference_population"
POOLED = Path("/xdisk/clairemcwhite/clairemcwhite/go_ec_positive_embeddings/pooled")
OUT = Path("/xdisk/clairemcwhite/clairemcwhite/go_ec_cavs_C1e-4_sharedscaler")
PROFAB = Path("/xdisk/clairemcwhite/proFAB")
C, TOL, MAX_ITER = 1e-4, 1e-10, 100000


def acc_of(header):
    p = header[1:].split()[0].split("|")
    return p[1] if len(p) >= 2 and p[0] in ("sp", "tr") else p[0]


def fasta_seqs(path):
    recs, h, seq = {}, None, []
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            if h: recs[acc_of(h)] = "".join(seq)
            h, seq = line, []
        else:
            seq.append(line.strip())
    if h: recs[acc_of(h)] = "".join(seq)
    return recs


def load_pooled():
    ids, mats = [], []
    for ids_f in sorted(POOLED.glob("chunk_*.ids")):
        ids += [l.strip() for l in open(ids_f) if l.strip()]
        mats.append(np.load(ids_f.with_suffix(".pos.npy")))
    M = np.vstack(mats)
    assert M.shape[0] == len(ids)
    return {k: i for i, k in enumerate(ids)}, M


def all_spans():
    return sorted(glob.glob(str(PROFAB / "GO_sets/go_dataset_part*/*/random_positive_train_max1000.span"))) + \
           sorted(glob.glob(str(PROFAB / "EC_sets/ec_dataset_part*/ecNo_*/random_positive_train_max1000.span")))


def train_one(span, index, M, Xneg_s, scaler):
    span = Path(span)
    rel = span.parent.relative_to(PROFAB)
    d = OUT / rel / "random_positive_train_max1000_cav"
    if (d / "report_v1.json").exists():
        return "skip"
    recs = fasta_seqs(str(span) + ".fasta")
    rows, n_nofasta, n_dropped = [], 0, 0
    for line in open(span):
        if not line.strip():
            continue
        a = line.split()[0]
        if a not in recs:
            n_nofasta += 1; continue
        k = "U" + hashlib.md5(recs[a].encode()).hexdigest()[:20]
        if k not in index:
            n_dropped += 1; continue
        rows.append(index[k])
    pos = M[rows].astype(np.float64)
    X = np.vstack([scaler.transform(pos), Xneg_s])
    y = np.r_[np.ones(len(pos)), np.zeros(len(Xneg_s))]
    lr = LogisticRegression(C=C, class_weight="balanced", solver="lbfgs", tol=TOL,
                            max_iter=MAX_ITER, random_state=42).fit(X, y)
    w = lr.coef_.ravel() / (np.linalg.norm(lr.coef_) + 1e-10)
    shipped_dir = span.parent / "random_positive_train_max1000_cav"
    try:
        npos_shipped = json.load(open(shipped_dir / "report_v1.json"))["concept_metrics"]["n_positive"]
    except Exception:
        npos_shipped = None
    shipped = shipped_dir / "concept_v1.npy"
    cos_shipped = float(w @ np.load(shipped).ravel()) if shipped.exists() else None
    d.mkdir(parents=True, exist_ok=True)
    np.save(d / "concept_v1.npy", w.astype(np.float32))
    json.dump({
        "term_dir": str(rel), "span": str(span),
        "concept_metrics": {"n_positive": len(pos), "n_negative": len(Xneg_s),
                            "train_auroc": float(roc_auc_score(y, lr.decision_function(X)))},
        "n_positive_shipped": npos_shipped, "n_positive_matches_shipped": npos_shipped == len(pos),
        "n_span_ids_without_fasta": n_nofasta, "n_dropped_gt2048": n_dropped,
        "config": {"C": C, "tol": TOL, "max_iter": MAX_ITER, "class_weight": "balanced", "solver": "lbfgs",
                   "random_state": 42, "scaler": str(REF / "scaler_v1.pkl"), "negatives": str(REF / "neg_10000.pkl")},
        "n_iter": int(lr.n_iter_[0]), "hit_max_iter": int(lr.n_iter_[0]) >= MAX_ITER,
        "cos_to_shipped": cos_shipped, "date": time.strftime("%Y-%m-%d"),
    }, open(d / "report_v1.json", "w"), indent=1)
    return "ok"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--n-shards", type=int, default=1)
    ap.add_argument("--spans", nargs="*", help="restrict to these span files (testing)")
    a = ap.parse_args()
    scaler = joblib.load(REF / "scaler_v1.pkl")
    assert scaler.n_samples_seen_ == 10000
    neg, _ = load_sequence_embeddings(str(REF / "neg_10000.pkl"))
    Xneg_s = scaler.transform(np.asarray(neg, dtype=np.float64))
    index, M = load_pooled()
    spans = (a.spans or all_spans())[a.shard::a.n_shards]
    print(f"shard {a.shard}/{a.n_shards}: {len(spans)} terms; pooled rows {len(index)}", flush=True)
    fails = OUT / "failed" / f"shard_{a.shard}.tsv"
    fails.parent.mkdir(parents=True, exist_ok=True)
    counts, t0 = {"ok": 0, "skip": 0, "fail": 0}, time.time()
    for i, s in enumerate(spans, 1):
        try:
            counts[train_one(s, index, M, Xneg_s, scaler)] += 1
        except Exception as e:
            counts["fail"] += 1
            with open(fails, "a") as fh:
                fh.write(f"{s}\t{type(e).__name__}: {str(e).splitlines()[0] if str(e) else ''}\n")
            traceback.print_exc()
        if i % 50 == 0:
            print(f"{i}/{len(spans)} {counts} {(time.time() - t0) / i:.1f}s/term", flush=True)
    print(f"done {counts}", flush=True)


if __name__ == "__main__":
    main()
