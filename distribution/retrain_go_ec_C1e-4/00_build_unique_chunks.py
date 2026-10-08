#!/usr/bin/env python
"""
GO/EC retrain at C=1e-4 (2026-09-29), step 0: deduplicate the positives.

Every GO/EC random_positive_train_max1000.span is whole-sequence (IDs only, no span
coordinates), so a protein's pooled positive embedding does not depend on which term
it came from.  The 2.30M positive rows over 7,631 terms are 303k unique proteins.
Deduplicate by SEQUENCE, not accession: 131 accessions carry different sequences in
different term files (UniProt versions differ between sets), so each term must use the
sequence in its own .span.fasta.  Each unique sequence gets the ID U<md5[:20]> (the
embedding depends only on the sequence, not the header); step 2 maps every term's span
accession -> that term's fasta sequence -> md5 -> pooled row.  Chunks of 500 records
plus a whole-sequence span file per chunk; step 1 embeds with the per-term flags unchanged.
Out: /xdisk/clairemcwhite/clairemcwhite/go_ec_positive_embeddings/chunks/chunk_NNNN.{fasta,span}
"""
import glob, hashlib
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

OUT = Path("/xdisk/clairemcwhite/clairemcwhite/go_ec_positive_embeddings/chunks")
CHUNK = 500
spans = sorted(glob.glob("/xdisk/clairemcwhite/proFAB/GO_sets/go_dataset_part*/*/random_positive_train_max1000.span")) + \
        sorted(glob.glob("/xdisk/clairemcwhite/proFAB/EC_sets/ec_dataset_part*/ecNo_*/random_positive_train_max1000.span"))

def acc_of(header):
    p = header[1:].split()[0].split("|")
    return p[1] if len(p) >= 2 and p[0] in ("sp", "tr") else p[0]

def read(span):
    fa = span + ".fasta"
    ids = [l.split()[0] for l in open(span) if l.strip()]
    recs, h, seq = {}, None, []
    for line in open(fa):
        line = line.rstrip("\n")
        if line.startswith(">"):
            if h: recs[acc_of(h)] = (h, "".join(seq))
            h, seq = line, []
        else:
            seq.append(line.strip())
    if h: recs[acc_of(h)] = (h, "".join(seq))
    return span, ids, recs

with ThreadPoolExecutor(32) as ex:
    R = list(ex.map(read, spans))
uniq, missing, n_rows = {}, 0, 0
for span, ids, recs in R:
    for a in ids:
        if a not in recs:
            missing += 1; continue
        n_rows += 1
        seq = recs[a][1]
        uniq.setdefault("U" + hashlib.md5(seq.encode()).hexdigest()[:20], seq)
OUT.mkdir(parents=True, exist_ok=True)
accs = sorted(uniq)
for c in range(0, len(accs), CHUNK):
    n = c // CHUNK
    with open(OUT / f"chunk_{n:04d}.fasta", "w") as fa, open(OUT / f"chunk_{n:04d}.span", "w") as sp:
        for a in accs[c:c + CHUNK]:
            fa.write(f">{a}\n{uniq[a]}\n"); sp.write(f"{a}\n")
lens = sorted(len(s) for s in uniq.values())
print(f"{len(spans)} terms; {n_rows} positive rows; {len(accs)} unique sequences; span IDs with no fasta record: {missing}; "
      f"{(len(accs) + CHUNK - 1) // CHUNK} chunks of {CHUNK}; median length {lens[len(lens)//2]}, "
      f"> 2048 aa: {sum(l > 2048 for l in lens)}; total residues {sum(lens):,}")
