#!/usr/bin/env python
"""Classify every GO CAV by which scaler frame it was TRAINED in.

A CAV that fit its own scaler standardised (positives + 10k negatives), so its
StandardScaler carries n_samples_seen_ = 10000 + n_positive > 10000.  A CAV built
with --scaler-pkl reuses the shared reference scaler, n_samples_seen_ == 10000.
That test is exact; md5 is not (joblib writes differing bytes for identical content).
"""
import glob, json, os, sys
from concurrent.futures import ProcessPoolExecutor
import joblib

PATTERN = "/xdisk/clairemcwhite/proFAB/GO_sets/go_dataset_part*/*/random_positive_train_max1000_cav"

def classify(cav_dir):
    span = cav_dir[: -len("_cav")] + ".span"
    try:
        n_seen = float(joblib.load(os.path.join(cav_dir, "scaler_v1.pkl")).n_samples_seen_)
    except Exception as e:
        return (cav_dir, span, "ERROR", -1, -1, str(e)[:80])
    try:
        rep = json.load(open(os.path.join(cav_dir, "report_v1.json")))
        n_pos = rep["concept_metrics"]["n_positive"]
    except Exception:
        n_pos = -1
    frame = "own" if n_seen > 10000 else "shared"
    return (cav_dir, span, frame, n_seen, n_pos, "")

def main():
    dirs = sorted(glob.glob(PATTERN))
    print(f"{len(dirs)} GO CAV directories", file=sys.stderr)
    with ProcessPoolExecutor(max_workers=16) as ex:
        rows = list(ex.map(classify, dirs, chunksize=32))
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "go_scaler_frames.tsv")
    with open(out, "w") as fh:
        fh.write("cav_dir\tspan\tframe\tn_samples_seen\tn_positive\terror\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")
    from collections import Counter
    c = Counter(r[2] for r in rows)
    print(dict(c), file=sys.stderr)
    # consistency check: own-frame scalers must satisfy n_seen == 10000 + n_positive
    bad = [r for r in rows if r[2] == "own" and r[4] >= 0 and r[3] != 10000 + r[4]]
    print(f"own-frame rows violating n_seen == 10000 + n_positive: {len(bad)}", file=sys.stderr)
    print(out, file=sys.stderr)

if __name__ == "__main__":
    main()
