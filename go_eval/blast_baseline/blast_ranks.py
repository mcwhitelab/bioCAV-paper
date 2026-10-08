#!/usr/bin/env python
"""
Homology-transfer (DIAMOND nearest-neighbour) baseline for Fig 2C, 2026-09-29.

For each validation protein q and each GO term t with a CAV, the BLAST score is
the best DIAMOND bitscore of q against t's CAV training positives
(random_positive_train_max1000.span.fasta), self-hits removed, 0 if no hit
(e <= 1e-3, --more-sensitive). The true term is then ranked among the same
candidates as the CAV/DeepGo-SE *_rank_ont_filt_pf columns
(rank_all_terms/compute_ranks_all_terms.py): same-ontology CAV terms, minus the
true term's ancestors, the protein's other GPAD terms (ancestor-closed) and
the terms proFAB lists it as a positive for (ancestor-closed).
Obsolete GO terms are not candidates (parse_obo drops them; 2026-09-29).

Ties: a protein's best hit is usually a positive for many related terms at
once, so BLAST scores tie far more often than CAV scores. Two ranks are kept:
  blast_rank_opt = 1 + #{score > true}                 (optimistic, as CAV/DS)
  blast_rank_mid = 1 + #{score > true} + #{ties}/2     (tie-averaged)
"No prediction" for BLAST = no hit to any training positive of the true term.

Inputs: val_vs_train.tsv, term_members.tsv (written alongside), the staged
go_specificity_ranks_allterms_{ont}.csv. Writes blast_ranks_{ont}.csv here.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "rank_all_terms"))
from compute_ranks_all_terms import (  # noqa: E402
    GO_EVAL, OBO, ONTS, ancestors, load_profab_positives, parse_obo,
)

STAGE = HERE.parents[1] / "figures" / "figure_data_staging" / "noknown"
HITS = ["qseqid", "sseqid", "pident", "length", "qlen", "slen", "evalue", "bitscore"]


def main():
    ns, parents = parse_obo(OBO)
    members = pd.read_csv(HERE / "term_members.tsv", sep="\t", names=["go_term", "sseqid"])
    terms = sorted(members["go_term"].unique())
    col = {g: i for i, g in enumerate(terms)}
    term_ns = np.array([ns.get(g) for g in terms], dtype=object)

    import os
    hits = pd.read_csv(os.environ.get("HITS_FILE", HERE / "val_vs_train.tsv"), sep="\t", names=HITS)
    hits = hits[hits.qseqid != hits.sseqid]
    best = hits.groupby(["qseqid", "sseqid"], as_index=False)[["bitscore", "pident"]].max()
    qt = best.merge(members, on="sseqid")
    qt = qt.groupby(["qseqid", "go_term"], as_index=False).agg(
        bitscore=("bitscore", "max"), pident=("pident", "max"))
    print(f"{len(hits)} hits -> {len(qt)} (protein, term) scores")
    by_q = {q: d for q, d in qt.groupby("qseqid")}

    gpad = pd.concat([
        pd.read_csv(GO_EVAL / "static" / "go_annotations" / f"{o}.gpad", sep="\t",
                    comment="!", header=None, usecols=[0, 1, 3], names=["db", "pid", "go"])
        for o in ONTS])
    gpad = gpad[gpad["db"] == "UniProtKB"]
    gpad["pid"] = gpad["pid"].astype(str).str.split("-").str[0]
    known_true = {}
    for pid, gs in gpad.groupby("pid")["go"]:
        s = set(gs)
        for t in list(s):
            s |= ancestors(t, parents)
        known_true[pid] = s

    pairs = {o: pd.read_csv(STAGE / f"go_specificity_ranks_allterms_{o}.csv") for o in ONTS}
    profab_pos = load_profab_positives(set().union(*[set(p.protein_id) for p in pairs.values()]))

    for ont in ONTS:
        ont_mask = term_ns == ont
        rows = []
        for r in pairs[ont].itertuples(index=False):
            p, g = r.protein_id, r.go_term
            score = np.zeros(len(terms))
            pid_ = np.zeros(len(terms))
            if p in by_q:
                d = by_q[p]
                idx = d.go_term.map(col).values
                score[idx] = d.bitscore.values
                pid_[idx] = d.pident.values
            drop = ancestors(g, parents) | known_true.get(p, set())
            for t in profab_pos.get(p, ()):
                drop |= {t} | ancestors(t, parents)
            drop -= {g}
            keep = ont_mask.copy()
            keep[[col[t] for t in drop if t in col]] = False
            j = col[g]
            cand, true = score[keep], score[j]
            gt = int((cand > true).sum())
            ties = int((cand == true).sum()) - 1  # the true term is in cand
            rows.append({"protein_id": p, "go_term": g,
                         "blast_true_bits": true, "blast_true_pident": pid_[j],
                         "blast_predicted": true > 0,
                         "blast_rank_opt": gt + 1, "blast_rank_mid": gt + 1 + ties / 2,
                         "blast_n_ties": ties, "n_cand": int(keep.sum())})
        out = pd.DataFrame(rows)
        out.to_csv(HERE / f"blast_ranks_{ont}.csv", index=False)
        pr = out.blast_predicted
        print(f"[{ont}] {len(out)} pairs | BLAST predicted {100*pr.mean():.1f}% | "
              f"rank1 opt {100*((out.blast_rank_opt == 1) & pr).mean():.1f}% "
              f"mid {100*((out.blast_rank_mid < 1.5) & pr).mean():.1f}% | "
              f"top10 opt {100*((out.blast_rank_opt <= 10) & pr).mean():.1f}% "
              f"mid {100*((out.blast_rank_mid <= 10) & pr).mean():.1f}%")


if __name__ == "__main__":
    main()
