#!/usr/bin/env python
"""
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


def fit(X, y, tol, max_iter, solver="lbfgs"):
    lr = LogisticRegression(C=1.0, max_iter=max_iter, random_state=42, solver=solver,
                            class_weight="balanced", tol=tol).fit(X, y)
    w = lr.coef_[0]
    return w / (np.linalg.norm(w) + 1e-10), int(lr.n_iter_[0])


def main():
    terms = pd.read_csv(HERE / "terms.tsv", sep="\t")
    neg, _ = load_sequence_embeddings(str(REF / "neg_10000.pkl"))
    shared = joblib.load(REF / "scaler_v1.pkl")
    assert shared.n_samples_seen_ == 10000
    for arm in ["default", "converged", "own_conv"]:
        (OUT / "cavs" / arm).mkdir(parents=True, exist_ok=True)
    rows, newton_done = [], {}
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
        t0 = time.time()
        w_def, it_def = fit(Xs, y, 1e-4, 1000)
        w_conv, it_conv = fit(Xs, y, 1e-10, 100000)
        own = StandardScaler().fit(X)
        w_own, it_own = fit(own.transform(X), y, 1e-10, 100000)
        rec = {"library": r.library, "ontology": r.ontology, "term": r.term, "name": name,
               "n_positive_shipped": r.n_positive, "n_positive_now": len(pos),
               "iter_default": it_def, "iter_converged": it_conv, "iter_own_conv": it_own,
               "cos_default_vs_shipped": float(w_def @ shipped),
               "cos_converged_vs_shipped": float(w_conv @ shipped),
               "cos_converged_vs_default": float(w_conv @ w_def),
               "cos_own_conv_vs_converged": float(w_own @ w_conv),
               "cos_own_conv_vs_shipped": float(w_own @ shipped),
               "cos_converged_vs_newton": np.nan}
        if newton_done.get(r.library, 0) < 3:
            w_nt, _ = fit(Xs, y, 1e-12, 1000, solver="newton-cholesky")
            rec["cos_converged_vs_newton"] = float(w_conv @ w_nt)
            newton_done[r.library] = newton_done.get(r.library, 0) + 1
        rec["seconds"] = round(time.time() - t0, 1)
        np.save(OUT / "cavs" / "default" / f"{name}.npy", w_def.astype(np.float32))
        np.save(OUT / "cavs" / "converged" / f"{name}.npy", w_conv.astype(np.float32))
        np.save(OUT / "cavs" / "own_conv" / f"{name}.npy", w_own.astype(np.float32))
        joblib.dump(own, OUT / "cavs" / "own_conv" / f"{name}.scaler.pkl")
        rows.append(rec)
        print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in rec.items()
               if k in ("name", "iter_default", "iter_converged", "cos_default_vs_shipped",
                        "cos_converged_vs_shipped", "cos_own_conv_vs_converged")}, flush=True)
    pd.DataFrame(rows).to_csv(HERE / "fits.tsv", sep="\t", index=False)
    print(f"wrote {len(rows)} fits -> fits.tsv")


if __name__ == "__main__":
    main()
