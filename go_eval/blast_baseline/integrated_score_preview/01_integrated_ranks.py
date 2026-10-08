#!/usr/bin/env python
"""
Preview (2026-10-06): an integrated CAV x BLAST score for the Fig 2C GO ranking.

Same pairs, candidates and filters as the live Fig 2C (C=1e-4 library, known-positive-
excluded pairs, non-obsolete same-ontology terms, strict _pf filter), rebuilt here with
the full per-candidate score vectors so the true term can be re-ranked under a combined
score. CAV and BLAST ranks are recomputed and asserted equal to the live
figures/figure_data/go_rank_families_{ont}.csv before anything new is computed.

Per protein, over its candidate terms:
  CAV   c_t  = shared-scaler-standardised embedding . unit concept vector (raw projection)
  BLAST b_t  = best DIAMOND bitscore to term t's CAV training positives, 0 if no hit

Integrated scores (both are products; neither is tuned on these pairs, see eps note):
  rankprod : pct(c_t) * pct(b_t), pct = tie-averaged percentile among the candidates
             (1 = best). Parameter-free. Terms with no hit share the bottom BLAST
             percentile block, so any hit outranks no hit at equal CAV percentile.
  soft     : exp(z(c_t)) * (eps + b_t / max_t b_t), z = per-protein z of c over the
             candidates. i.e. log score = z(c) + log(eps + b/bmax). eps = 0.1 fixed a priori
             (a best-hit term can gain at most log(11) = 2.4 CAV SDs). eps 0.01 / 1 are
             also written as a sensitivity check; they are NOT used to pick the variant.
  For a protein with no BLAST hit to any candidate, both reduce to the CAV ranking.

Ranks are tie-averaged ("mid", as BLAST in Fig 2C). "No prediction" for the integrated
score = CAV gate fails (LLR <= 0) AND no BLAST hit to the true term, i.e. neither
component makes a call.
Writes integrated_ranks_{ont}.csv here.
"""
import glob
import os
import sys
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from scipy.stats import rankdata

PAPER = Path("/groups/clairemcwhite/cav_workspace/bioCAV-paper")
GE = PAPER / "go_eval"
HERE = Path(__file__).resolve().parent
STAGE = PAPER / "figures" / "figure_data_staging" / "C1e-4"
LIVE = PAPER / "figures" / "figure_data"
BASE = Path("/xdisk/clairemcwhite/clairemcwhite/go_ec_cavs_C1e-4_overlay/GO_sets")
EXTRA = GE / "blast_baseline" / "C1e-4_extra"
SUB = "random_positive_train_max1000_cav"
HITS = ["qseqid", "sseqid", "pident", "length", "qlen", "slen", "evalue", "bitscore"]
EPS = [0.1, 0.01, 1.0]

os.environ["FIG_DATA_DIR"] = str(STAGE)
sys.path.insert(0, str(GE / "rank_all_terms"))
import compute_ranks_all_terms as C  # noqa: E402
C.GO_SETS = BASE


def mid_rank_desc(x):
    """1 = largest; ties averaged."""
    return rankdata(-x, method="average")


def main():
    ns, parents = C.parse_obo(C.OBO)
    dirs = sorted(glob.glob(str(BASE / "go_dataset_part*" / "GO_*" / SUB)))
    go_ids = [Path(d).parent.name.replace("GO_", "GO:") for d in dirs]
    W = np.stack([np.load(Path(d) / "concept_v1.npy").ravel().astype(np.float64) for d in dirs])
    col = {g: i for i, g in enumerate(go_ids)}
    cav_ns = np.array([ns.get(g) for g in go_ids], dtype=object)
    live = cav_ns != None  # noqa: E711
    print(f"{len(go_ids)} CAVs, {live.sum()} live")

    gpad = pd.concat([
        pd.read_csv(GE / "static" / "go_annotations" / f"{o}.gpad", sep="\t", comment="!",
                    header=None, usecols=[0, 1, 3], names=["db", "pid", "go"]) for o in C.ONTS])
    gpad = gpad[gpad["db"] == "UniProtKB"]
    gpad["pid"] = gpad["pid"].astype(str).str.split("-").str[0]
    known_true = {}
    for pid, gs in gpad.groupby("pid")["go"]:
        s = set(gs)
        for t in list(s):
            s |= C.ancestors(t, parents)
        known_true[pid] = s

    members = pd.read_csv(EXTRA / "term_members.tsv", sep="\t", names=["go_term", "sseqid"])
    hits = pd.read_csv(EXTRA / "val_vs_train.tsv", sep="\t", names=HITS)
    hits = hits[hits.qseqid != hits.sseqid]
    best = hits.groupby(["qseqid", "sseqid"], as_index=False)[["bitscore", "pident"]].max()
    qt = best.merge(members, on="sseqid").groupby(["qseqid", "go_term"], as_index=False).agg(
        bitscore=("bitscore", "max"), pident=("pident", "max"))
    qt = qt[qt.go_term.isin(col)]
    by_q = {q: d for q, d in qt.groupby("qseqid")}

    pairs = {o: pd.read_csv(STAGE / f"go_specificity_ranks_{o}.csv") for o in C.ONTS}
    profab_pos = C.load_profab_positives(set().union(*[set(p.protein_id) for p in pairs.values()]))
    scaler = joblib.load(C.SCALER)

    for ont in C.ONTS:
        fam = pd.read_csv(LIVE / f"go_rank_families_{ont}.csv")
        prs = pairs[ont].drop_duplicates(["protein_id", "go_term"])
        embs, sids = C.load_sequence_embeddings(
            GE / "results" / f"temporal_eval_{ont}" / "eval_embeddings" / "val_proteins.span.fasta.pkl")
        idx = {}
        for i, s in enumerate(sids):
            idx.setdefault(C.normalise_pid(s), i)
        X = scaler.transform(embs)
        ont_mask = cav_ns == ont
        rows = []
        for r in prs.itertuples(index=False):
            p, g = r.protein_id, r.go_term
            j = col[g]
            c = X[idx[p]] @ W.T
            b = np.zeros(len(go_ids)); pi = np.zeros(len(go_ids))
            if p in by_q:
                d = by_q[p]
                ii = d.go_term.map(col).values
                b[ii] = d.bitscore.values; pi[ii] = d.pident.values
            drop = (C.ancestors(g, parents) | known_true.get(p, set())) - {g}
            for t in profab_pos.get(p, ()):
                drop |= {t} | C.ancestors(t, parents)
            drop -= {g}
            keep = live & ont_mask
            keep[[col[t] for t in drop if t in col]] = False
            cand = np.where(keep)[0]
            jj = int(np.searchsorted(cand, j))
            assert cand[jj] == j
            cc, bb = c[cand], b[cand]
            n = len(cand)
            row = {"protein_id": p, "go_term": g, "n_cand": n,
                   "cav_rank": int((cc > cc[jj]).sum()) + 1,
                   "blast_rank_mid": float(mid_rank_desc(bb)[jj]),
                   "blast_true_bits": bb[jj], "blast_true_pident": pi[j],
                   "blast_max_bits": bb.max(), "blast_n_hit_terms": int((bb > 0).sum())}
            # rank product of tie-averaged percentiles
            pc = 1 - (mid_rank_desc(cc) - 1) / n
            pb = 1 - (mid_rank_desc(bb) - 1) / n
            row["rankprod_rank"] = float(mid_rank_desc(pc * pb)[jj])
            # soft product
            z = (cc - cc.mean()) / cc.std()
            bn = bb / bb.max() if bb.max() > 0 else bb
            for e in EPS:
                s = z + np.log(e + bn)
                row[f"soft_eps{e:g}_rank"] = float(mid_rank_desc(s)[jj])
            # oracle: best of the two components (ceiling for any combination)
            row["oracle_min_rank"] = min(row["cav_rank"], row["blast_rank_mid"])
            rows.append(row)
        out = pd.DataFrame(rows).merge(
            fam[["protein_id", "go_term", "n_families", "llr", "cav_rank_ont_filt_pf",
                 "tool_predicted", "tool_rank_ont_filt_pf", "blast_predicted", "blast_rank_mid"]]
            .rename(columns={"blast_rank_mid": "blast_rank_mid_live"}),
            on=["protein_id", "go_term"], validate="one_to_one")
        assert len(out) == len(fam), (len(out), len(fam))
        ok_c = (out.cav_rank == out.cav_rank_ont_filt_pf).mean()
        ok_b = np.isclose(out.blast_rank_mid, out.blast_rank_mid_live).mean()
        ok_p = ((out.blast_true_bits > 0) == out.blast_predicted).mean()
        print(f"[{ont}] {len(out)} pairs; reproduce live CAV rank {ok_c:.4f}, BLAST mid {ok_b:.4f}, "
              f"BLAST hit flag {ok_p:.4f}", flush=True)
        assert ok_c == 1 and ok_b == 1 and ok_p == 1
        out.drop(columns=["cav_rank_ont_filt_pf", "blast_rank_mid_live"]).to_csv(
            HERE / f"integrated_ranks_{ont}.csv", index=False)


if __name__ == "__main__":
    main()
