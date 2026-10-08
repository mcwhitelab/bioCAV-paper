#!/usr/bin/env python
"""Train one CAV from a pooled, truncated positive set with the step-2 settings
(02_train_go_ec_C1e-4.py). Used for EC 3.4.22.28; see 03_ec_3-4-22-28_truncated.sh.
Usage: 03_train_truncated.py WORK_DIR SPAN"""
import importlib.util, json, sys, time
from pathlib import Path

import joblib
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("t", HERE / "02_train_go_ec_C1e-4.py")
t = importlib.util.module_from_spec(spec); spec.loader.exec_module(t)

work, span = Path(sys.argv[1]), Path(sys.argv[2])
pos = np.load(work / "pos.npy").astype(np.float64)
n_ids = sum(1 for l in open(work / "pos.span") if l.strip())
assert pos.shape == (n_ids, 1152), (pos.shape, n_ids)

scaler = joblib.load(t.REF / "scaler_v1.pkl")
assert scaler.n_samples_seen_ == 10000
neg, _ = t.load_sequence_embeddings(str(t.REF / "neg_10000.pkl"))
Xneg_s = scaler.transform(np.asarray(neg, dtype=np.float64))
X = np.vstack([scaler.transform(pos), Xneg_s])
y = np.r_[np.ones(len(pos)), np.zeros(len(Xneg_s))]
lr = LogisticRegression(C=t.C, class_weight="balanced", solver="lbfgs", tol=t.TOL,
                        max_iter=t.MAX_ITER, random_state=42).fit(X, y)
w = lr.coef_.ravel() / (np.linalg.norm(lr.coef_) + 1e-10)

rel = span.parent.relative_to(t.PROFAB)
d = t.OUT / rel / "random_positive_train_max1000_cav"
d.mkdir(parents=True, exist_ok=True)
np.save(d / "concept_v1.npy", w.astype(np.float32))
report = {
    "term_dir": str(rel), "span": str(span),
    "concept_metrics": {"n_positive": len(pos), "n_negative": len(Xneg_s),
                        "train_auroc": float(roc_auc_score(y, lr.decision_function(X)))},
    "n_positive_shipped": None, "n_positive_matches_shipped": False,
    "truncated_to": 2046,
    "note": "all positives > 2048 aa; embedded truncated to first 2046 aa (hf_embed -t 2046). "
            "No shipped CAV exists for this term.",
    "config": {"C": t.C, "tol": t.TOL, "max_iter": t.MAX_ITER, "class_weight": "balanced",
               "solver": "lbfgs", "random_state": 42, "scaler": str(t.REF / "scaler_v1.pkl"),
               "negatives": str(t.REF / "neg_10000.pkl")},
    "n_iter": int(lr.n_iter_[0]), "hit_max_iter": int(lr.n_iter_[0]) >= t.MAX_ITER,
    "cos_to_shipped": None, "date": time.strftime("%Y-%m-%d"),
}
json.dump(report, open(d / "report_v1.json", "w"), indent=1)
print(json.dumps({k: report[k] for k in ["concept_metrics", "n_iter", "hit_max_iter"]}))
print(f"-> {d}")
