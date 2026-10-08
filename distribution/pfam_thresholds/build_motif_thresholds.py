#!/usr/bin/env python
"""Rebuild motif_thresholds.json for the 27,472-motif Pfam CAV library.

The 27k library ships no motif_thresholds.json — only the 7,692 library has one
(covering 7,111 motifs).  But the file was never computed separately: in
batch_train_cavs.py it is a pure aggregation of fields that train_cav.py already
wrote into each L25_report_v1.json.  The 27k reports all carry those fields, so
the file can simply be re-aggregated.  Nothing is retrained.
"""
import argparse, json
from pathlib import Path

KEYS = ("decision_threshold", "f1_optimal_threshold", "f1_optimal_f1_train",
        "intercept", "coef_norm", "pos_score_mean", "pos_score_std",
        "neg_score_mean", "neg_score_std")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cav-dir", default="/xdisk/clairemcwhite/shamail/tcav_outputs_esmplusplus_all/cavs")
    ap.add_argument("--layer", type=int, default=25)
    ap.add_argument("--model", default="ESMplusplus_large")
    ap.add_argument("--version", default="v1")
    ap.add_argument("--out", default="motif_thresholds.json")
    a = ap.parse_args()

    per_motif, incomplete, unreadable = {}, [], []
    for d in sorted(Path(a.cav_dir).iterdir()):
        rep = d / f"L{a.layer}_report_{a.version}.json"
        if not rep.exists():
            unreadable.append(d.name); continue
        try:
            cm = json.loads(rep.read_text())["concept_metrics"]
        except Exception:
            unreadable.append(d.name); continue
        if not all(k in cm for k in KEYS):
            incomplete.append(d.name); continue
        per_motif[d.name] = {f"L{a.layer}": {k: cm[k] for k in KEYS}}

    Path(a.out).write_text(json.dumps({
        "model": a.model, "layers": [a.layer],
        "n_motifs": len(per_motif), "per_motif": per_motif}, indent=2))

    print(f"wrote {a.out}: {len(per_motif)} motifs")
    print(f"  no/unreadable report: {len(unreadable)}  incomplete metrics: {len(incomplete)}")
    if unreadable:
        Path("motifs_without_thresholds.txt").write_text("\n".join(sorted(unreadable + incomplete)) + "\n")
        print("  listed in motifs_without_thresholds.txt")

if __name__ == "__main__":
    main()
