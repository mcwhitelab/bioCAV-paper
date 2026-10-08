#!/usr/bin/env python
"""
GO/EC C-scan version (2026-09-29) of 03_evaluate.py: arms from 02b_fit_cscan.py, inputs
terms_cscan.tsv / fits_cscan.tsv, outputs per_term_cscan.tsv / per_pair_cscan.tsv.  own_conv_* arms are
scored in their own frame and left out of the ranks.  Caveat as below: ranks swap only the true
term's CAV (competitors stay shipped/early-stopped).

Original docstring follows.
GO/EC convergence pilot, step 3: do the headline metrics move?

For every sampled term and every arm (shipped concept_v1.npy, default refit,
converged, own_conv) this scores, as the eval scripts do (standardise, dot with the
unit CAV; own_conv is standardised with its own scaler, the others with the shared
scaler_v1.pkl):

  per-term  AUROC / AUPR of the term's validation positives vs its proFAB
            random_negative_test proteins
              GO: leak-free Fig 2 pairs (temporal_eval_{ont}_noknown), embeddings
                  temporal_eval_{ont}/eval_embeddings/val_proteins.span.fasta.pkl
              EC: the 209 Fig 3A pairs (ec_specificity_ranks_samelevel.csv),
                  embeddings EC-Bench/test_ec.fasta.pkl
  per-pair  LLR against Gaussians fit to the term's random test positives and
            negatives (eval_go_temporal.py / eval_ec.py; numpy std, ddof 0), and
            whether it is a prediction (LLR > 0)
  per-pair  rank of the true term (GO: cav_rank_ont_filt_pf via
            go_eval/rank_all_terms/compute_ranks_all_terms.py; EC:
            cav_rank_level_filt as compute_ec_ranks_same_level.py), with ONLY the
            true term's CAV swapped for the refit, every competitor shipped.

The swap is the limit of a pilot: a full retrain moves the competitors too, so the
rank effect of retraining everything is not measured here.  own_conv is left out of
the ranks because the rank scripts score every candidate in the shared frame.

Writes per_term.tsv, per_pair.tsv and summary.txt.
"""
import glob, os, sys
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[1]
OUT = Path("/xdisk/clairemcwhite/clairemcwhite/converge_pilot_go_ec")
sys.path.insert(0, str(PAPER.parent / "bioCAV" / "core"))
from src.utils.data_loader import load_sequence_embeddings  # noqa: E402

SHARED = joblib.load(PAPER / "shared_protein" / "reference_population" / "scaler_v1.pkl")
CS = ["1", "0.1", "0.01", "0.001", "0.0001"]
ARMS = ["shipped", "default"] + [f"conv_C{c}" for c in CS] + ["own_conv_C1", "own_conv_C0.0001"]
RANK_ARMS = ["shipped", "default"] + [f"conv_C{c}" for c in CS]


SKIP = {"sp", "tr", "sw", "ref"}


def norm_pid(pid):   # as normalise_pid in both rank scripts
    parts = [p for p in str(pid).split()[0].split("|") if p and p not in SKIP]
    return parts[0] if parts else str(pid)


def load_pkl(path):
    e, ids = load_sequence_embeddings(str(path))
    return e, [norm_pid(i) for i in ids]


def vec(arm, r):
    if arm == "shipped":
        return np.load(Path(r.cav_dir) / "concept_v1.npy").ravel(), SHARED
    w = np.load(OUT / "cavs_cscan" / arm / f"{r.name}.npy").ravel()
    sc = joblib.load(OUT / "cavs_cscan" / arm / f"{r.name}.scaler.pkl") if arm.startswith("own_") else SHARED
    return w, sc


def logn(x, mu, sd):
    sd = max(sd, 1e-6)
    return -0.5 * np.log(2 * np.pi * sd ** 2) - (x - mu) ** 2 / (2 * sd ** 2)


def main():
    terms = pd.read_csv(HERE / "terms_cscan.tsv", sep="\t")
    fits = pd.read_csv(HERE / "fits_cscan.tsv", sep="\t")
    terms = terms.merge(fits[["term", "name"]], on="term")

    # validation positives per term
    val = {}
    for ont in ["mf", "bp", "cc"]:
        res = pd.read_csv(PAPER / f"go_eval/results/temporal_eval_{ont}_noknown/eval_temporal_results.tsv", sep="\t")
        for t, g in res.groupby("go_term"):
            val[t] = sorted(set(map(norm_pid, g.protein_id)))
    ec_pairs = pd.read_csv(PAPER / "figures/figure_data/ec_specificity_ranks_samelevel.csv", dtype={"ec_number": str})
    for t, g in ec_pairs.groupby("ec_number"):
        val[t] = sorted(set(g.protein_id))
    val_emb = {ont: load_pkl(PAPER / f"go_eval/results/temporal_eval_{ont}/eval_embeddings/val_proteins.span.fasta.pkl")
               for ont in ["mf", "bp", "cc"]}
    ec_emb = load_pkl(PAPER / "ec_eval/static/ec_annotations/EC-Bench/test_ec.fasta.pkl")

    term_rows, pair_rows = [], []
    for r in terms.itertuples(index=False):
        base = Path(r.span).parent
        tp_e, _ = load_pkl(base / "random_positive_test_max1000.span.fasta.pkl")
        tn_e, _ = load_pkl(base / "random_negative_test_max1000.span.fasta.pkl")
        ve, vids = val_emb[r.ontology] if r.library == "GO" else ec_emb
        idx = {}
        for i, p in enumerate(vids):
            idx.setdefault(p, i)
        prots = [p for p in val[r.term] if p in idx]
        V = ve[[idx[p] for p in prots]]
        for arm in ARMS:
            w, sc = vec(arm, r)
            sp, sn, sv = sc.transform(tp_e) @ w, sc.transform(tn_e) @ w, sc.transform(V) @ w
            y = np.r_[np.ones(len(sv)), np.zeros(len(sn))]
            s = np.r_[sv, sn]
            llr = logn(sv, sp.mean(), sp.std()) - logn(sv, sn.mean(), sn.std())
            term_rows.append({"library": r.library, "ontology": r.ontology, "term": r.term, "arm": arm,
                              "n_val": len(sv), "n_test_neg": len(sn),
                              "auroc": roc_auc_score(y, s), "aupr": average_precision_score(y, s),
                              "neg_mean": sn.mean(), "neg_std": sn.std(), "val_median": np.median(sv),
                              "pct_llr_pos": 100 * (llr > 0).mean()})
            for p, x, l in zip(prots, sv, llr):
                pair_rows.append({"library": r.library, "ontology": r.ontology, "term": r.term,
                                  "protein_id": p, "arm": arm, "score": x, "llr": l})
    per_term = pd.DataFrame(term_rows)
    per_pair = pd.DataFrame(pair_rows)

    # ---- ranks, true-term CAV swapped ----
    rank_rows = []
    # GO: rerun the Fig 2C rank script with the pilot terms' CAVs substituted
    stage = HERE / "rank_stage_cscan"
    stage.mkdir(exist_ok=True)
    go_terms = set(terms[terms.library == "GO"].term)
    for ont in ["mf", "bp", "cc"]:
        d = pd.read_csv(PAPER / f"figures/figure_data/go_specificity_ranks_{ont}.csv")
        d[d.go_term.isin(go_terms)].to_csv(stage / f"go_specificity_ranks_{ont}.csv", index=False)
    os.environ["FIG_DATA_DIR"] = str(stage)
    sys.path.insert(0, str(PAPER / "go_eval" / "rank_all_terms"))
    import compute_ranks_all_terms as cr  # noqa: E402
    orig_load = cr.load_cavs
    base_ids, base_mat, base_frames = orig_load()
    go_rows = terms[terms.library == "GO"]
    for arm in RANK_ARMS:
        mat = base_mat.copy()
        col = {g: i for i, g in enumerate(base_ids)}
        if arm != "shipped":
            for r in go_rows.itertuples(index=False):
                mat[col[r.term]] = vec(arm, r)[0]
        cr.load_cavs = lambda m=mat: (base_ids, m, base_frames)
        cr.main()
        for ont in ["mf", "bp", "cc"]:
            o = pd.read_csv(stage / f"go_specificity_ranks_allterms_{ont}.csv")
            for x in o.itertuples(index=False):
                rank_rows.append({"library": "GO", "ontology": ont, "term": x.go_term, "protein_id": x.protein_id,
                                  "arm": arm, "rank": x.cav_rank_ont_filt_pf})
    # EC: same-level filtered rank, as compute_ec_ranks_same_level.py
    EC_SETS = Path("/xdisk/clairemcwhite/proFAB/EC_sets")
    ids, vecs = [], []
    for d in sorted(glob.glob(str(EC_SETS / "ec_dataset_part*" / "ecNo_*"))):
        c = Path(d) / "random_positive_train_max1000_cav" / "concept_v1.npy"
        if c.exists():
            ids.append(Path(d).name[len("ecNo_"):].replace("-", ".")); vecs.append(np.load(c))
    ec_mat = np.stack(vecs); level = np.array([len(e.split(".")) for e in ids]); col = {e: i for i, e in enumerate(ids)}
    gold = pd.read_csv(PAPER / "ec_eval/static/ec_annotations/EC-Bench/test_ec_long.tsv", sep="\t")
    strip = lambda e: ".".join(p for p in str(e).split(".") if p not in ("-", ""))  # noqa: E731
    true_ecs = gold.assign(ec=gold["ec_number"].map(strip)).groupby("protein_id")["ec"].apply(set).to_dict()
    eidx = {}
    for i, p in enumerate(ec_emb[1]):
        eidx.setdefault(p, i)
    ec_rows = terms[terms.library == "EC"]
    for arm in RANK_ARMS:
        mat = ec_mat.copy()
        if arm != "shipped":
            for r in ec_rows.itertuples(index=False):
                mat[col[r.term]] = vec(arm, r)[0]
        S = SHARED.transform(ec_emb[0]) @ mat.T
        for x in ec_pairs.itertuples(index=False):
            i, j = eidx[x.protein_id], col[x.ec_number]
            lvl = len(x.ec_number.split("."))
            others = {".".join(e.split(".")[:lvl]) for e in true_ecs.get(x.protein_id, ())
                      if len(e.split(".")) >= lvl} - {x.ec_number}
            keep = level == lvl
            keep[[col[e] for e in others if e in col]] = False
            rk = int((S[i][keep] > S[i][j]).sum()) + 1
            if arm == "shipped":
                assert rk == x.cav_rank_level_filt, (x.protein_id, x.ec_number, rk, x.cav_rank_level_filt)
            rank_rows.append({"library": "EC", "ontology": f"level{lvl}", "term": x.ec_number,
                              "protein_id": x.protein_id, "arm": arm, "rank": rk})
    ranks = pd.DataFrame(rank_rows)
    per_pair = per_pair.merge(ranks, on=["library", "ontology", "term", "protein_id", "arm"], how="left")
    per_term.to_csv(HERE / "per_term_cscan.tsv", sep="\t", index=False)
    per_pair.to_csv(HERE / "per_pair_cscan.tsv", sep="\t", index=False)
    print("wrote per_term_cscan.tsv, per_pair_cscan.tsv")


if __name__ == "__main__":
    main()
