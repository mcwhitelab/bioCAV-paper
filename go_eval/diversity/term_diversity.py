#!/usr/bin/env python
"""
Sequence diversity of each Fig 2 GO term's CAV training positives, and sequence
novelty of each validation protein relative to them (exploratory, 2026-09-29).

For every term in the known-positive-excluded Fig 2 set
(figures/figure_data_staging/noknown/temporal_tool_comparison_{ont}.csv):

  term level  — DIAMOND all-vs-all of the term's random_positive_train_max1000
                FASTA (the exact CAV training positives). Families = connected
                components of the graph with an edge wherever a hit has
                >= 30% identity over >= 50% of the shorter sequence.
                Reports n_train, n_families, families_per_100, largest_family_frac.
  pair level  — each kept validation protein searched against the same FASTA;
                nn_identity = best % identity among hits with >= 50% query
                coverage (0 when there is no such hit, e-value <= 1e-3).

Writes term_diversity_{ont}.csv and pair_novelty_{ont}.csv next to this script.
"""
import glob
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[1]
GO_EVAL = PAPER / "go_eval"
STAGE = PAPER / "figures" / "figure_data_staging" / "noknown"
GO_SETS = Path(os.environ.get("PROFAB_GO_SETS", "/xdisk/clairemcwhite/proFAB/GO_sets"))
DIAMOND = os.environ.get("DIAMOND", "/opt/ohpc/pub/apps/diamond/2.0.9/diamond")
TRAIN_FASTA = "random_positive_train_max1000.span.fasta"
MIN_ID, MIN_COV = 30.0, 50.0
SKIP = {"sp", "tr", "sw", "ref"}
FMT = ["qseqid", "sseqid", "pident", "qlen", "slen", "length", "evalue"]


def norm(pid):
    parts = [p for p in str(pid).split()[0].split("|") if p and p not in SKIP]
    return parts[0] if parts else str(pid)


def read_fasta(path):
    seqs, cur = {}, None
    for line in open(path):
        line = line.rstrip()
        if line.startswith(">"):
            cur = norm(line[1:])
            seqs[cur] = []
        elif cur:
            seqs[cur].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


def write_fasta(seqs, path):
    with open(path, "w") as fh:
        for k, s in seqs.items():
            fh.write(f">{k}\n{s}\n")


def diamond(query, db, out):
    subprocess.run([DIAMOND, "blastp", "-q", query, "-d", db, "-o", out,
                    "--more-sensitive", "--max-target-seqs", "0", "--evalue", "1e-3",
                    "--threads", "1", "--quiet", "--outfmt", "6", *FMT], check=True)
    if os.path.getsize(out) == 0:
        return pd.DataFrame(columns=FMT)
    return pd.read_csv(out, sep="\t", header=None, names=FMT)


def families(ids, hits):
    parent = {i: i for i in ids}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    cov = 100 * hits["length"] / hits[["qlen", "slen"]].min(axis=1)
    keep = hits[(hits.pident >= MIN_ID) & (cov >= MIN_COV) & (hits.qseqid != hits.sseqid)]
    for a, b in zip(keep.qseqid, keep.sseqid):
        if a in parent and b in parent:
            parent[find(a)] = find(b)
    sizes = pd.Series([find(i) for i in ids]).value_counts()
    return len(sizes), sizes.iloc[0] / len(ids)


def find_go_dir(go_id):
    for d in sorted(glob.glob(str(GO_SETS / "go_dataset_part*" / go_id.replace(":", "_")))):
        return Path(d)
    return None


def main():
    for ont in ["mf", "bp", "cc"]:
        terms = pd.read_csv(STAGE / f"temporal_tool_comparison_{ont}.csv")["go_term"]
        res = pd.read_csv(GO_EVAL / "results" / f"temporal_eval_{ont}_noknown" /
                          "eval_temporal_results.tsv", sep="\t")
        res["p"] = res["protein_id"].map(norm)
        val_seqs = read_fasta(GO_EVAL / "results" / f"temporal_eval_{ont}" /
                              "eval_embeddings" / "val_proteins.span.fasta")
        term_rows, pair_rows = [], []
        for k, go_id in enumerate(terms, 1):
            train = read_fasta(find_go_dir(go_id) / TRAIN_FASTA)
            with tempfile.TemporaryDirectory() as tmp:
                tfa, db = f"{tmp}/train.fa", f"{tmp}/train"
                write_fasta(train, tfa)
                subprocess.run([DIAMOND, "makedb", "--in", tfa, "-d", db, "--quiet"], check=True)
                nfam, largest = families(list(train), diamond(tfa, db, f"{tmp}/self.tsv"))
                term_rows.append({"go_term": go_id, "n_train": len(train), "n_families": nfam,
                                  "families_per_100": 100 * nfam / len(train),
                                  "largest_family_frac": largest})
                val = {p: val_seqs[p] for p in res.loc[res.go_term == go_id, "p"] if p in val_seqs}
                vfa = f"{tmp}/val.fa"
                write_fasta(val, vfa)
                h = diamond(vfa, db, f"{tmp}/val.tsv")
                h = h[100 * h["length"] / h["qlen"] >= MIN_COV]
                best = h.groupby("qseqid")["pident"].max()
                for p in val:
                    pair_rows.append({"protein_id": p, "go_term": go_id,
                                      "nn_identity": float(best.get(p, 0.0))})
            if k % 50 == 0:
                print(f"[{ont}] {k}/{len(terms)} terms", flush=True)
        pd.DataFrame(term_rows).to_csv(HERE / f"term_diversity_{ont}.csv", index=False)
        pd.DataFrame(pair_rows).to_csv(HERE / f"pair_novelty_{ont}.csv", index=False)
        print(f"[{ont}] done: {len(term_rows)} terms, {len(pair_rows)} pairs", flush=True)


if __name__ == "__main__":
    sys.exit(main())
