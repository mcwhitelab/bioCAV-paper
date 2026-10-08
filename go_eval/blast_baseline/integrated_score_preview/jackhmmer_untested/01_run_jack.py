#!/usr/bin/env python
"""jackhmmer for the discordant (CAV rank 1, BLAST >10 / no hit) pairs that the strict callout
run never tested. Recipe unchanged from ../../callouts_C1e-4/02_run_jack.py (pyhmmer, 3 iterations,
E=10, incE=1e-3, db = blast_baseline/train_positives_all.fa). Three already-tested strict queries
are added as a calibration check against ../../callouts_C1e-4/jack_vs_train.tsv."""
import sys
from pathlib import Path
import pandas as pd, pyhmmer
HERE = Path(__file__).resolve().parent; BB = HERE.parents[1]
d = pd.read_csv(HERE.parent / "cav_rank1_blast_low.csv")
untested = d[d.homology_status.isna()]
untested[["ontology", "protein_id", "go_term"]].to_csv(HERE / "untested_pairs.csv", index=False)
CAL = ["Q8IEJ0", "Q9VTG7", "Q2RVS1"]
want = set(untested.protein_id) | set(CAL)
alpha = pyhmmer.easel.Alphabet.amino()
with pyhmmer.easel.SequenceFile(str(BB / "val_queries.fa"), digital=True, alphabet=alpha) as f:
    qs = [s for s in f if (s.name.decode() if isinstance(s.name, bytes) else s.name) in want]
missing = want - {(s.name.decode() if isinstance(s.name, bytes) else s.name) for s in qs}
print(len(untested), "pairs;", untested.protein_id.nunique(), "proteins;", len(qs), "queries; missing", missing, flush=True)
if "--dry" in sys.argv: sys.exit()
with pyhmmer.easel.SequenceFile(str(BB / "train_positives_all.fa"), digital=True, alphabet=alpha) as f:
    db = f.read_block()
rows = []
for q, res in zip([s.name for s in qs], pyhmmer.hmmer.jackhmmer(qs, db, cpus=1, max_iterations=3, E=10, incE=1e-3)):
    q = q.decode() if isinstance(q, bytes) else q
    for h in res.hits:
        n = h.name.decode() if isinstance(h.name, bytes) else h.name
        rows.append((q, n, 0.0, 0, 0, 0, h.evalue, h.score))
    print(q, "iter", res.iteration, "hits", len(res.hits), flush=True)
pd.DataFrame(rows).to_csv(HERE / "jack_vs_train.tsv", sep="\t", header=False, index=False)
