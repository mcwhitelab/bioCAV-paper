#!/usr/bin/env python
"""
GO/EC C-scan (2026-09-29; subset terms_cscan.tsv, embedded by 01b_embed_cscan_cpu.sbatch).
Adapted from 02_fit_cavs.py.  Arms: default (C=1, tol 1e-4, the shipped recipe), conv_C<c>
(shared scaler, tol 1e-10) for C in 1, 0.1, 0.01, 0.001, 0.0001, and own_conv_C1 / own_conv_C0.0001
(own per-term scaler).  Writes cavs_cscan/<arm>/<NAME>.npy and fits_cscan.tsv.

Original docstring follows.
GO/EC convergence pilot, step 2: refit each sampled term's CAV three ways.

Inputs: pooled train positives from 01_embed_positives.sbatch, the shared 10k
background (neg_10000.pkl) and the shared scaler (scaler_v1.pkl) -- the same
inputs train_cav_from_span.sh used for the shipped CAVs.

Arms (all C=1, lbfgs, class_weight='balanced', random_state=42, coef_ unit-normalised
exactly as core/src/train_cav.py::extract_cav_vector):
  default     shared scaler, tol=1e-4 (sklearn default), max_iter=1000 -- the shipped
              recipe; its cosine to the shipped concept_v1.npy checks the reproduction
  converged   shared scaler, tol=1e-10, max_iter=100000
  own_conv    StandardScaler fit on this term's positives + the 10k background
              (create_scaler), tol=1e-10, max_iter=100000.  Its vector lives in its
              own frame; it is also saved with its scaler for scoring.
For the first 3 terms of each library the converged fit is checked against an exact
newton-cholesky solve.

Writes /xdisk/.../converge_pilot_go_ec/cavs/<arm>/<NAME>.npy (+ own scalers) and fits.tsv.
"""
import sys, time, warnings
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore", category=UserWarning)
sys.path.insert(0, "/groups/clairemcwhite/cav_workspace/bioCAV/core")
from src.utils.data_loader import load_sequence_embeddings  # noqa: E402

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[1]
REF = PAPER / "shared_protein" / "reference_population"
OUT = Path("/xdisk/clairemcwhite/clairemcwhite/converge_pilot_go_ec")


def name_of(lib, term):
    return f"{lib}_{term.replace(':', '_').replace('.', '_')}"


def fit(X, y, tol, max_iter, solver="lbfgs", C=1.0):
    lr = LogisticRegression(C=C, max_iter=max_iter, random_state=42, solver=solver,
                            class_weight="balanced", tol=tol).fit(X, y)
    w = lr.coef_[0]
    return w / (np.linalg.norm(w) + 1e-10), int(lr.n_iter_[0])


CS = [1.0, 0.1, 0.01, 0.001, 0.0001]
OWN_CS = [1.0, 0.0001]
CAVS = OUT / "cavs_cscan"


def main():
    terms = pd.read_csv(HERE / "terms_cscan.tsv", sep="\t")
    neg, _ = load_sequence_embeddings(str(REF / "neg_10000.pkl"))
    shared = joblib.load(REF / "scaler_v1.pkl")
    assert shared.n_samples_seen_ == 10000
    arms = ["default"] + [f"conv_C{c:g}" for c in CS] + [f"own_conv_C{c:g}" for c in OWN_CS]
    for arm in arms:
        (CAVS / arm).mkdir(parents=True, exist_ok=True)
    rows = []
    for r in terms.itertuples(index=False):
        name = name_of(r.library, r.term)
        pos_f = OUT / "pos" / f"{name}.pos.npy"
        if not pos_f.exists():
            print(f"MISSING {name}", flush=True)
            continue
        pos = np.load(pos_f)
        X = np.vstack([pos, neg]); y = np.hstack([np.ones(len(pos)), np.zeros(len(neg))])
        shipped = np.load(Path(r.cav_dir) / "concept_v1.npy").ravel()
        Xs = shared.transform(X)
        own = StandardScaler().fit(X); Xo = own.transform(X)
        rec = {"library": r.library, "ontology": r.ontology, "term": r.term, "name": name,
               "n_positive_shipped": r.n_positive, "n_positive_now": len(pos)}
        W = {}
        W["default"], rec["iter_default"] = fit(Xs, y, 1e-4, 1000)
        for c in CS:
            W[f"conv_C{c:g}"], rec[f"iter_conv_C{c:g}"] = fit(Xs, y, 1e-10, 100000, C=c)
        for c in OWN_CS:
            W[f"own_conv_C{c:g}"], rec[f"iter_own_conv_C{c:g}"] = fit(Xo, y, 1e-10, 100000, C=c)
        for arm, w in W.items():
            rec[f"cos_{arm}_vs_shipped"] = float(w @ shipped)
            np.save(CAVS / arm / f"{name}.npy", w.astype(np.float32))
            if arm.startswith("own_"):
                joblib.dump(own, CAVS / arm / f"{name}.scaler.pkl")
        rows.append(rec)
        print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in rec.items() if k == "name" or k.startswith("cos_")}, flush=True)
    pd.DataFrame(rows).to_csv(HERE / "fits_cscan.tsv", sep="\t", index=False)
    print(f"wrote {len(rows)} fits -> fits_cscan.tsv")


if __name__ == "__main__":
    main()
