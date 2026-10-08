#!/usr/bin/env python3
"""
Refresh the clan_acc / clan_name columns of the motif tables to Pfam 38.2.

The tables were annotated from Pfam 38.1 (27,481 families, 2025-09); the current
release is 38.2 (30,134 families, 2026-01).  Only the two clan columns are
touched: the clustering, the cut columns and the family names are left exactly
as they were, so no repacking is needed -- this is a relabel.

Accessions retired in 38.2 (21 of them) keep their 38.1 clan, since 38.2 has
nothing to say about them.  Families that lost their clan get an empty one.
Every change is written to a log CSV.
"""
import csv, gzip, sys, shutil
from pathlib import Path

CLANS_38_2 = sys.argv[1]
ROOT = Path(__file__).resolve().parent
TARGETS = sorted(ROOT.glob("results/dendrogram_comparison/*.csv")) + [
    ROOT / "results/figures/circlepack_input_diffusion.csv",
    ROOT / "results/figures/pfam_clusters_27k.csv",
]

opener = gzip.open if CLANS_38_2.endswith(".gz") else open
new = {r[0]: (r[1], r[2]) for r in csv.reader(opener(CLANS_38_2, "rt"), delimiter="\t")}
print(f"Pfam 38.2: {len(new):,} families")

log, seen = [], set()
for path in TARGETS:
    with open(path, newline="") as fh:
        rows = list(csv.reader(fh))
    head = rows[0]
    if not {"accession", "clan_acc", "clan_name"} <= set(head):
        continue                           # summary tables (per_clan, scores, cuts)
    i_acc, i_clan, i_name = (head.index(c) for c in ("accession", "clan_acc", "clan_name"))

    n_chg = n_gone = 0
    for r in rows[1:]:
        acc = r[i_acc]
        if acc not in new:
            n_gone += 1
            continue
        clan, name = new[acc]
        if (r[i_clan], r[i_name]) != (clan, name):
            n_chg += 1
            if acc not in seen:            # the tables share one annotation
                seen.add(acc)
                log.append([acc, r[i_clan], r[i_name], clan, name])
            r[i_clan], r[i_name] = clan, name

    tmp = path.with_suffix(path.suffix + ".tmp")
    with open(tmp, "w", newline="") as fh:
        csv.writer(fh, lineterminator="\n").writerows(rows)
    shutil.move(tmp, path)                 # atomic, so a crash cannot half-write
    print(f"{path.relative_to(ROOT)}: {n_chg:,} rows relabelled, "
          f"{n_gone} retired accessions left on 38.1")

out = ROOT / "results/figures/clan_changes_38.1_to_38.2.csv"
with open(out, "w", newline="") as fh:
    w = csv.writer(fh, lineterminator="\n")
    w.writerow(["accession", "clan_acc_38_1", "clan_name_38_1",
                "clan_acc_38_2", "clan_name_38_2"])
    w.writerows(sorted(log))
print(f"\n{len(log):,} distinct families changed clan -> {out.relative_to(ROOT)}")
