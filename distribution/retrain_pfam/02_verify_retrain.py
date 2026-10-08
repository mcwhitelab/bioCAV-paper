#!/usr/bin/env python
"""Verify the shared-scaler Pfam retrain.  Expect 27,472 complete, 0 problems."""
import json, sys
from pathlib import Path
import joblib, numpy as np
OUT = Path("/xdisk/clairemcwhite/clairemcwhite/pfam_cavs_shared_scaler")
EMB = Path("/xdisk/clairemcwhite/shamail/tcav_outputs_esmplusplus_all/embeddings/ESMplusplus_large")
motifs = sorted(p.name for p in EMB.iterdir() if p.name.startswith("PF"))
bad, dropped = {}, {}
for m in motifs:
    d = OUT / "cavs" / m
    try:
        r = json.load(open(d / "L25_report_v1.json"))
        assert "retrain" in r, "no retrain block"
        assert joblib.load(d / "L25_scaler_v1.pkl").n_samples_seen_ == 10000, "wrong scaler"
        assert np.isfinite(np.load(d / "L25_concept_v1.npy")).all(), "non-finite CAV"
        if r["concept_metrics"]["n_dropped_truncated_pos"]:
            dropped[m] = r["concept_metrics"]["n_dropped_truncated_pos"]
    except Exception as e:
        bad[m] = f"{type(e).__name__}: {e}"
fails = [l for f in (OUT / "failed").glob("*.tsv") for l in open(f)]
print(f"{len(motifs)} motifs, {len(motifs) - len(bad)} complete, {len(bad)} problems, "
      f"{len(fails)} failure-log lines")
print(f"truncated positives dropped: {sum(dropped.values())} rows in {len(dropped)} families")
for m, e in list(bad.items())[:20]:
    print("  ", m, e)
sys.exit(1 if bad or fails else 0)
