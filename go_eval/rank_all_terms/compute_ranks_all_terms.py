#!/usr/bin/env python
"""
Fig 2C alternative: rank the true GO term against ALL trained GO CAVs.

The shipped panel C (compare_tool_temporal.py::make_rank_scatter) ranks each
(validation protein, true term) pair only among the terms *evaluated* in that
ontology (MF 226 / BP 395 / CC 258). This script re-ranks the same pairs
against wider candidate sets, identically for CAV and DeepGo-SE:

  rank_all  — all trained GO CAVs, pooled across MF/BP/CC
  rank_ont  — all trained GO CAVs in the pair's own namespace (per go.obo)

Every candidate set excludes CAVs whose term is obsolete in go.obo (484 of 6,677)
or absent from it (122), i.e. 6,071 live terms. Obsolete terms have no is_a/part_of
parents, so the ancestor filter cannot remove them, and many are broad grouping terms
(e.g. "obsolete cellular metabolic process") that DeepGo-SE, trained on an earlier
release, scores highly for most proteins. Excluded 2026-09-29; before that they were
candidates (MF 1,161 / BP 4,687 / CC 707 -> MF 1,114 / BP 4,297 / CC 660).
  *_filt_pf — rank_ont_filt, additionally removing every term for which proFAB
              already lists the protein as a positive (plus ancestors); added
              2026-09-29 because the GPADs hold only IDA/IPI-descendant evidence
  *_filt    — rank_all / rank_ont with the true term's is_a/part_of ancestors
              and every other term the protein is annotated with in the GPADs
              (any date, ancestor-closed) removed from the candidates
              ("filtered rank"). Needed because DeepGo-SE scores are
              hierarchy-consistent (ancestor >= descendant), so unfiltered it is
              outranked by the true term's own ancestors by construction.
  rank_eval — the original evaluated-term set, recomputed with the CAVs now on
              disk (sanity check against the shipped go_specificity_ranks_*.csv)

CAV score = shared-scaler-standardised embedding · unit concept vector (as in
make_rank_scatter). DeepGo-SE score = its emitted score, 0.0 where it emitted
nothing; predictions for all three ontologies are pooled for each protein (every
validation FASTA was predicted in all three).

The pair list and the "no prediction" flags (llr, tool_predicted) are taken
unchanged from the shipped go_specificity_ranks_{ont}.csv.
"""
import glob
import gzip
import os
import sys
from pathlib import Path

import joblib
import numpy as np
import pandas as pd

PAPER = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PAPER.parent / "bioCAV" / "core"))
from src.utils.data_loader import load_sequence_embeddings  # noqa: E402

GO_SETS = Path("/xdisk/clairemcwhite/proFAB/GO_sets")
GO_EVAL = PAPER / "go_eval"
# FIG_DATA_DIR overrides where go_specificity_ranks_{ont}.csv is read and the
# allterms CSVs are written (e.g. a staging dir for a preview run).
FIG_DATA = Path(os.environ.get("FIG_DATA_DIR", PAPER / "figures" / "figure_data"))
SCALER = PAPER / "shared_protein" / "reference_population" / "scaler_v1.pkl"
OBO = GO_EVAL / "static" / "go_annotations" / "go.obo"
CAV_SUBDIR = "random_positive_train_max1000_cav"
ONTS = ["mf", "bp", "cc"]
NS = {"molecular_function": "mf", "biological_process": "bp", "cellular_component": "cc"}
SKIP = {"sp", "tr", "sw", "ref"}


def normalise_pid(pid):
    parts = [p for p in str(pid).split() [0].split("|") if p and p not in SKIP]
    return parts[0] if parts else str(pid)


def parse_obo(path, drop_obsolete=True):
    """Namespaces and is_a/part_of parents (the relations GO annotations propagate over).
    With drop_obsolete, obsolete terms get namespace None, which takes them out of every
    per-ontology candidate mask."""
    ns, parents, cur, alts, obsolete = {}, {}, None, [], set()
    for line in open(path):
        line = line.strip()
        if line == "is_obsolete: true" and cur:
            obsolete.add(cur)
        elif line.startswith("is_a: GO:") and cur:
            parents.setdefault(cur, set()).add(line[6:16])
        elif line.startswith("relationship: part_of GO:") and cur:
            parents.setdefault(cur, set()).add(line[22:32])
        elif line == "[Term]" or line.startswith("["):
            cur, alts = None, []
        elif line.startswith("id: GO:"):
            cur = line[4:]
        elif line.startswith("alt_id: GO:"):
            alts.append(line[8:])
        elif line.startswith("namespace:") and cur:
            n = NS.get(line.split(":", 1)[1].strip())
            ns[cur] = n
            for a in alts:
                ns.setdefault(a, n)
                parents.setdefault(a, set()).add(cur)
    if drop_obsolete:
        for t in obsolete:
            ns[t] = None
    return ns, parents


def ancestors(term, parents, _cache={}):
    if term not in _cache:
        out = set()
        for p in parents.get(term, ()):
            out |= {p} | ancestors(p, parents)
        _cache[term] = out
    return _cache[term]


def load_cavs():
    ids, vecs, frames = [], [], []
    for d in sorted(glob.glob(str(GO_SETS / "go_dataset_part*" / "GO_*"))):
        cav = Path(d) / CAV_SUBDIR / "concept_v1.npy"
        if not cav.exists():
            continue
        ids.append(Path(d).name.replace("GO_", "GO:"))
        vecs.append(np.load(cav))
        sc = joblib.load(Path(d) / CAV_SUBDIR / "scaler_v1.pkl")
        frames.append("own" if sc.n_samples_seen_ > 10000 else "shared")
    return ids, np.stack(vecs), frames


def load_tool(ont):
    """DeepGo-SE predictions for this ontology's validation proteins, all 3 aspects."""
    files = [GO_EVAL / "outside_tools" / f"val_proteins.span_preds_{ont}.tsv"]
    files += [GO_EVAL / "results" / f"temporal_eval_{ont}" / "eval_embeddings" /
              f"val_proteins.span_preds_{o}.tsv.gz" for o in ONTS if o != ont]
    dfs = []
    for f in files:
        dfs.append(pd.read_csv(f, sep="\t", header=None,
                               names=["protein_id", "go_term", "tool_score"]))
    df = pd.concat(dfs, ignore_index=True)
    df["protein_id"] = df["protein_id"].map(normalise_pid)
    df = df.groupby(["protein_id", "go_term"], as_index=False)["tool_score"].max()
    return df


def load_profab_positives(proteins):
    """GO terms (with a CAV) for which each protein is a proFAB random train/test
    positive. proFAB builds positives from full Swiss-Prot annotation with
    true-path propagation, so these are known-true terms our GPADs can miss."""
    out = {}
    for d in sorted(glob.glob(str(GO_SETS / "go_dataset_part*" / "GO_*"))):
        if not (Path(d) / CAV_SUBDIR).is_dir():
            continue
        g = Path(d).name.replace("GO_", "GO:")
        for f in ["random_positive_train_max1000.span", "random_positive_test_max1000.span"]:
            fp = Path(d) / f
            if fp.exists():
                for x in fp.read_text().split():
                    x = normalise_pid(x)
                    if x in proteins:
                        out.setdefault(x, set()).add(g)
    return out


def rank_of(scores, true_score, mask=None):
    if mask is not None:
        scores = scores[mask]
    return int((scores > true_score).sum()) + 1


def main():
    ns, parents = parse_obo(OBO)
    # every GO term each protein is annotated with in the three GPADs, any date
    # (pre- and post-cutoff), closed over ancestors
    gpad = pd.concat([
        pd.read_csv(GO_EVAL / "static" / "go_annotations" / f"{o}.gpad", sep="\t",
                    comment="!", header=None, usecols=[0, 1, 3], names=["db", "pid", "go"])
        for o in ONTS])
    gpad = gpad[gpad["db"] == "UniProtKB"]
    gpad["pid"] = gpad["pid"].astype(str).str.split("-").str[0]
    known_true = {}
    for pid, terms in gpad.groupby("pid")["go"]:
        s = set(terms)
        for t in list(s):
            s |= ancestors(t, parents)
        known_true[pid] = s
    go_ids, cav_matrix, frames = load_cavs()
    pair_proteins = set()
    for o in ONTS:
        pair_proteins |= set(pd.read_csv(FIG_DATA / f"go_specificity_ranks_{o}.csv")["protein_id"])
    profab_pos = load_profab_positives(pair_proteins)
    print(f"proFAB positives loaded for {len(profab_pos)} of {len(pair_proteins)} pair proteins")
    frames = pd.Series(frames)
    print(f"Loaded {len(go_ids)} GO CAVs; scaler frames: {frames.value_counts().to_dict()}")
    cav_ns = np.array([ns.get(g) for g in go_ids], dtype=object)
    print("CAV namespaces:", pd.Series(cav_ns).value_counts(dropna=False).to_dict())
    live = cav_ns != None  # noqa: E711  obsolete (484) or absent from go.obo (122)
    col = {g: i for i, g in enumerate(go_ids)}
    scaler = joblib.load(SCALER)

    for ont in ONTS:
        old = pd.read_csv(FIG_DATA / f"go_specificity_ranks_{ont}.csv")
        eval_terms = set(old["go_term"])
        eval_mask = np.array([g in eval_terms for g in go_ids])
        ont_mask = cav_ns == ont

        embs, sids = load_sequence_embeddings(
            GO_EVAL / "results" / f"temporal_eval_{ont}" / "eval_embeddings" /
            "val_proteins.span.fasta.pkl")
        idx = {}
        for i, s in enumerate(sids):
            idx.setdefault(normalise_pid(s), i)
        S = scaler.transform(embs) @ cav_matrix.T

        tool = load_tool(ont)
        tool = tool[tool["go_term"].isin(col)]
        T = np.zeros_like(S)
        for p, g, v in tool.itertuples(index=False):
            if p in idx:
                T[idx[p], col[g]] = v

        rows = []
        for r in old.itertuples(index=False):
            i, j = idx[r.protein_id], col[r.go_term]
            cs, ts = S[i], T[i]
            # filtered: drop the true term's ancestors and every other known-true term
            drop = (ancestors(r.go_term, parents) |
                    known_true.get(r.protein_id, set())) - {r.go_term}
            keep = live.copy()
            keep[[col[g] for g in drop if g in col]] = False
            # _pf: additionally drop terms proFAB already lists this protein as a
            # positive for (and their ancestors), keeping the true term
            pf = set()
            for g in profab_pos.get(r.protein_id, ()):
                pf |= {g} | ancestors(g, parents)
            pf -= {r.go_term}
            keep_pf = keep.copy()
            keep_pf[[col[g] for g in pf if g in col]] = False
            rows.append({
                "protein_id": r.protein_id, "go_term": r.go_term,
                "llr": r.llr, "tool_predicted": r.tool_predicted,
                "cav_rank_all":  rank_of(cs, cs[j], live),
                "tool_rank_all": rank_of(ts, ts[j], live),
                "cav_rank_ont":  rank_of(cs, cs[j], ont_mask),
                "tool_rank_ont": rank_of(ts, ts[j], ont_mask),
                "cav_rank_all_filt":  rank_of(cs, cs[j], keep),
                "tool_rank_all_filt": rank_of(ts, ts[j], keep),
                "cav_rank_ont_filt":  rank_of(cs, cs[j], keep & ont_mask),
                "tool_rank_ont_filt": rank_of(ts, ts[j], keep & ont_mask),
                "n_dropped_filt": int((live & ~keep).sum()),
                "cav_rank_ont_filt_pf":  rank_of(cs, cs[j], keep_pf & ont_mask),
                "tool_rank_ont_filt_pf": rank_of(ts, ts[j], keep_pf & ont_mask),
                "n_dropped_filt_pf": int((live & ~keep_pf).sum()),
                "cav_rank_eval":  rank_of(cs, cs[j], eval_mask),
                "tool_rank_eval": rank_of(ts, ts[j], eval_mask),
                "cav_rank_eval_shipped":  r.cav_rank,
                "tool_rank_eval_shipped": r.tool_rank,
                "n_terms_all": int(live.sum()), "n_terms_ont": int(ont_mask.sum()),
                "n_terms_eval": int(eval_mask.sum()),
            })
        out = pd.DataFrame(rows)
        dest = FIG_DATA / f"go_specificity_ranks_allterms_{ont}.csv"
        out.to_csv(dest, index=False)

        print(f"\n[{ont}] {len(out)} pairs -> {dest.name}  "
              f"(candidates: all {live.sum()}, ont {ont_mask.sum()}, eval {eval_mask.sum()})")
        print(f"  tool_rank_eval == shipped: {(out.tool_rank_eval == out.tool_rank_eval_shipped).mean():.3f}; "
              f"cav_rank_eval == shipped: {(out.cav_rank_eval == out.cav_rank_eval_shipped).mean():.3f}")
        for k in ["eval_shipped", "eval", "ont", "ont_filt", "ont_filt_pf", "all", "all_filt"]:
            c, t = out[f"cav_rank_{k}"], out[f"tool_rank_{k}"]
            print(f"  {k:13s} CAV rank1 {100*(c==1).mean():5.1f}%  top10 {100*(c<=10).mean():5.1f}%  median {c.median():6.0f} | "
                  f"DeepGo-SE rank1 {100*(t==1).mean():5.1f}%  top10 {100*(t<=10).mean():5.1f}%  median {t.median():6.0f}")


if __name__ == "__main__":
    main()
