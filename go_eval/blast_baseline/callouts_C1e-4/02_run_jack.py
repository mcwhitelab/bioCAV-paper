#!/usr/bin/env python
"""jackhmmer of the callout candidates vs all GO CAV training positives: the recipe of the old
callouts, recovered from its session transcript and unchanged (pyhmmer, 3 iterations, E=10,
incE=1e-3, database = blast_baseline/train_positives_all.fa, which also already contains every
GO:0032436 positive). Writes DIAMOND-style 8-col hits (pident 0) so blast_ranks.py can rank
terms by jackhmmer score."""
import re
from pathlib import Path
import pandas as pd, pyhmmer
HERE = Path(__file__).resolve().parent; BB = HERE.parent
want = set((HERE / "jack_queries.txt").read_text().split())
alpha = pyhmmer.easel.Alphabet.amino()
with pyhmmer.easel.SequenceFile(str(BB / "val_queries.fa"), digital=True, alphabet=alpha) as f:
    qs = [s for s in f if s.name in want]
with pyhmmer.easel.SequenceFile(str(BB / "train_positives_all.fa"), digital=True, alphabet=alpha) as f:
    db = f.read_block()
print(len(qs), "queries;", len(db), "targets;", "pyhmmer", pyhmmer.__version__, flush=True)
rows = []
for q, res in zip([s.name for s in qs], pyhmmer.hmmer.jackhmmer(qs, db, cpus=16, max_iterations=3, E=10, incE=1e-3)):
    n = res.hmm.name
    n = n.decode() if isinstance(n, bytes) else n
    assert n is None or re.sub(r"-i\d+$", "", n) == q, (n, q)
    for h in res.hits:
        rows.append((q, h.name, 0.0, 0, 0, 0, h.evalue, h.score))
    print(q, "iter", res.iteration, "hits", len(res.hits), flush=True)
pd.DataFrame(rows).to_csv(HERE / "jack_vs_train.tsv", sep="\t", header=False, index=False)
