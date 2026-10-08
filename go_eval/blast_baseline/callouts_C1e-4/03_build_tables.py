#!/usr/bin/env python
"""
Build the GO callout tables on the C=1e-4 library (2026-10-01) from 01_select.py's candidates
and 02_run_jack.py's jackhmmer hits. Same columns and logic as the old-library callouts
(go_eval/blast_baseline/callouts/, bp_cav_only_rank1_cases.csv), recovered from the session
that built them. Every CAV-derived column uses the C=1e-4 overlay, and every BLAST column uses
the extended hits/members (C1e-4_extra/).

  1. jackhmmer: hits with E <= 1e-3 ranked with blast_ranks.py logic (same candidates as
     Fig 2C) -> jack_ranks/; best E-value (any hit, E <= 10) to a true-term positive.
     Calibration: hits for the 60 old strict queries compared with the old run's file.
  2. Metadata: UniProt REST (gene, protein name, organism, length, Pfam xrefs) and the
     GPADs (evidence, reference, earliest date).
  3. CAV top-6 terms and the margin of the rank-1 term over the runner-up, with the Fig 2C
     candidate filter (filt_pf, same ontology, obsolete terms masked).
  4. Pfam domains shared with the true term's training positives (strict survivors without
     a jackhmmer hit, as before).
  5. Homology tier: jackhmmer hit > shared Pfam > none; callout = no detectable homology.
  6. Competitors for the cav_only set (top CAV runner-ups, DeepGo-SE and BLAST top terms,
     nearest DIAMOND neighbour, known GPAD terms in the ontology).
Outputs (this directory): strict_callouts.csv (all ontologies), mf_callouts.csv and
mf_callouts_cav_runnerups.csv (old format), {mf,bp,cc}_cav_only_rank1_cases.csv,
calibration.txt.
"""
import io, os, sys, time, urllib.request
from pathlib import Path

import joblib
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[2]
GE = PAPER / "go_eval"
BB = GE / "blast_baseline"
EXTRA = BB / "C1e-4_extra"
OVERLAY = Path("/xdisk/clairemcwhite/clairemcwhite/go_ec_cavs_C1e-4_overlay/GO_sets")
OLD_JACK = Path("/xdisk/clairemcwhite/clairemcwhite/go_train_positives_embed/jackhmmer/jack_vs_train.tsv")
FIG = PAPER / "figures" / "figure_data"
HITS = ["q", "s", "pident", "length", "qlen", "slen", "evalue", "bits"]

os.environ["FIG_DATA_DIR"] = str(FIG)
sys.path.insert(0, str(GE / "rank_all_terms")); sys.path.insert(0, str(BB))
import compute_ranks_all_terms as C  # noqa: E402
C.GO_SETS = OVERLAY  # C=1e-4 CAVs; proFAB-positive filter over the overlay's terms
import blast_ranks as Bl  # noqa: E402

strict = pd.read_csv(HERE / "strict_candidates.csv")
cav_only = pd.read_csv(HERE / "cav_only_candidates.csv")
members = pd.read_csv(EXTRA / "term_members.tsv", sep="\t", names=["go_term", "s"])
cal = []

# 1. jackhmmer ---------------------------------------------------------------------------
jh = pd.read_csv(HERE / "jack_vs_train.tsv", sep="\t", names=HITS)
old = pd.read_csv(OLD_JACK, sep="\t", names=HITS)
oq = set(old.q)
a = jh[jh.q.isin(oq)].sort_values(["q", "s"]).reset_index(drop=True)
b = old.sort_values(["q", "s"]).reset_index(drop=True)
same = len(a) == len(b) and (a[["q", "s"]].values == b[["q", "s"]].values).all() and \
    np.allclose(a.evalue, b.evalue, rtol=1e-6) and np.allclose(a.bits, b.bits, rtol=1e-6)
cal.append(f"jackhmmer hits for the {len(oq)} old queries identical to the old run: {same} "
           f"({len(a)} vs {len(b)} hits)")
(HERE / "jack_ranks").mkdir(exist_ok=True)
jh[jh.evalue <= 1e-3].to_csv(HERE / "jack_vs_train_e3.tsv", sep="\t", header=False, index=False)
members.to_csv(HERE / "jack_ranks" / "term_members.tsv", sep="\t", header=False, index=False)
os.environ["HITS_FILE"] = str(HERE / "jack_vs_train_e3.tsv")
Bl.HERE, Bl.STAGE = HERE / "jack_ranks", FIG
Bl.main()
jr = pd.concat([pd.read_csv(HERE / "jack_ranks" / f"blast_ranks_{o}.csv") for o in C.ONTS])
jr = jr.rename(columns=lambda x: x.replace("blast_", "jack_"))[
    ["protein_id", "go_term", "jack_predicted", "jack_rank_opt", "jack_rank_mid"]]
best = jh.merge(members, on="s").groupby(["q", "go_term"]).evalue.min() \
    .rename("jack_bestE_any").reset_index().rename(columns={"q": "protein_id"})

# 2. metadata ----------------------------------------------------------------------------
def uniprot(ids, fields):
    out = []
    ids = sorted(ids)
    for i in range(0, len(ids), 100):
        url = ("https://rest.uniprot.org/uniprotkb/accessions?accessions=" + ",".join(ids[i:i + 100])
               + f"&fields={fields}&format=tsv&size=100")
        for t in range(4):
            try:
                out.append(pd.read_csv(io.StringIO(urllib.request.urlopen(url, timeout=90).read().decode()), sep="\t"))
                break
            except Exception as e:  # noqa: BLE001
                print("retry", e, flush=True); time.sleep(5)
        else:
            raise RuntimeError(f"UniProt failed for batch {i}")
    return pd.concat(out)

allq = set(strict.protein_id) | set(cav_only.protein_id)
u = uniprot(allq, "accession,gene_primary,protein_name,organism_name,organism_id,length,xref_pfam")
u.columns = ["protein_id", "gene", "pname", "organism", "taxon", "length", "q_pfam"]
u["q_pfam"] = u.q_pfam.fillna("").str.strip(";")
g = pd.concat([pd.read_csv(GE / "static/go_annotations" / f"{o}.gpad", sep="\t", comment="!", header=None,
                           usecols=[0, 1, 3, 4, 8, 11], names=["db", "protein_id", "go_term", "ref", "date", "ev"],
                           dtype=str) for o in C.ONTS])
g = g[g.db == "UniProtKB"]
g["protein_id"] = g.protein_id.str.split("-").str[0]
g = g[g.protein_id.isin(allq)]
ann = g.groupby(["protein_id", "go_term"]).agg(
    ref=("ref", lambda s: ",".join(sorted(set(s)))), date=("date", "min"),
    ev=("ev", lambda s: ",".join(sorted(set(x.replace("goEvidence=", "") for x in s))))).reset_index()

# 3./6. CAV runner-ups and competitors ------------------------------------------------------
ns, parents = C.parse_obo(C.OBO)
names, cur = {}, None
for l in open(C.OBO):
    if l.startswith("id: GO:"): cur = l[4:].strip()
    elif l.startswith("name: ") and cur: names[cur] = l[6:].strip(); cur = None
ids, M, _ = C.load_cavs()
col = {t: i for i, t in enumerate(ids)}
cns = np.array([ns.get(t) for t in ids], dtype=object)
known = {}
for pid, ts in g.groupby("protein_id")["go_term"]:
    s = set(ts)
    for t in list(s): s |= C.ancestors(t, parents)
    known[pid] = s
pfpos = C.load_profab_positives(allq)
scaler = joblib.load(C.SCALER)
bh = pd.read_csv(EXTRA / "val_vs_train.tsv", sep="\t", names=HITS)
bh = bh[bh.q.isin(allq) & (bh.q != bh.s)]
bq = bh.groupby(["q", "s"], as_index=False).bits.max().merge(members, on="s") \
    .groupby(["q", "go_term"], as_index=False).bits.max()
nn = bh.sort_values("bits", ascending=False).drop_duplicates("q").set_index("q")

def fmt(scores, keep, k, nd):
    o = np.argsort(-np.where(keep, scores, -np.inf))[:k]
    return "; ".join(f"{names.get(ids[i], ids[i])} ({scores[i]:.{nd}g})" for i in o if keep[i] and scores[i] > 0) or "-"

extra = []
for ont in C.ONTS:
    sub = pd.concat([strict[strict.ont == ont], cav_only[cav_only.ont == ont]]).drop_duplicates(["protein_id", "go_term"])
    if sub.empty: continue
    embs, sids = C.load_sequence_embeddings(GE / "results" / f"temporal_eval_{ont}" / "eval_embeddings" /
                                            "val_proteins.span.fasta.pkl")
    idx = {}
    for i, s in enumerate(sids): idx.setdefault(C.normalise_pid(s), i)
    S = scaler.transform(embs) @ M.T
    tool = C.load_tool(ont); tool = tool[tool.go_term.isin(col) & tool.protein_id.isin(set(sub.protein_id))]
    for r in sub.itertuples(index=False):
        p, t0 = r.protein_id, r.go_term; i, j = idx[p], col[t0]
        drop = C.ancestors(t0, parents) | known.get(p, set())
        for t in pfpos.get(p, ()): drop |= {t} | C.ancestors(t, parents)
        drop -= {t0}
        keep = cns == ont; keep[[col[t] for t in drop if t in col]] = False
        sc = np.where(keep, S[i], -np.inf); o = np.argsort(-sc)[:6]
        assert ids[o[0]] == t0, (p, t0, ids[o[0]])  # rank 1 by construction
        kc = keep.copy(); kc[j] = False
        T = np.zeros(len(ids))
        for _, x in tool[tool.protein_id == p].iterrows(): T[col[x.go_term]] = x.tool_score
        B = np.zeros(len(ids))
        for _, x in bq[bq.q == p].iterrows():
            if x.go_term in col: B[col[x.go_term]] = x.bits
        extra.append(dict(
            protein_id=p, go_term=t0,
            cav_top6="; ".join(f"{names[ids[k]]} ({sc[k]:.1f})" for k in o), cav_margin=sc[o[0]] - sc[o[1]],
            cav_true=S[i, j], cav_next=fmt(S[i], kc, 2, 2), dg_true=T[j], dg_top=fmt(T, kc, 3, 2),
            blast_true=B[j], blast_top=fmt(B, kc, 3, 3),
            nn_hit=nn.s.get(p, ""), nn_pident=nn.pident.get(p, np.nan), nn_bits=nn.bits.get(p, np.nan),
            known_terms="; ".join(sorted(names.get(t, t) for t in set(g[g.protein_id == p].go_term)
                                         if ns.get(t) == ont))[:300]))
extra = pd.DataFrame(extra)

# 4./5. strict callouts ------------------------------------------------------------------------
st = strict.merge(jr, on=["protein_id", "go_term"], how="left").merge(best, on=["protein_id", "go_term"], how="left") \
    .merge(ann, on=["protein_id", "go_term"], how="left").merge(u, on="protein_id", how="left") \
    .merge(extra, on=["protein_id", "go_term"], how="left")
need = st[~st.jack_predicted.astype(bool)]
mem_need = members[members.go_term.isin(set(need.go_term))]
pf = uniprot(set(mem_need.s), "accession,xref_pfam"); pf.columns = ["acc", "pfam"]
P = {r.acc: set(str(r.pfam).strip(";").split(";")) - {"nan", ""} for r in pf.itertuples()}
st["n_pos"] = st.go_term.map(members.groupby("go_term").size())
def shared(r):
    if r.jack_predicted: return np.nan
    q = set(str(r.q_pfam).split(";")) - {"", "nan"}; cnt = {}
    for s in mem_need[mem_need.go_term == r.go_term].s:
        for d in q & P.get(s, set()): cnt[d] = cnt.get(d, 0) + 1
    return ";".join(f"{k}x{v}" for k, v in sorted(cnt.items()))
st["shared_pfam"] = st.apply(shared, axis=1)
def tier(r):
    if r.jack_predicted: return "remote homology: jackhmmer finds true-term positives"
    if isinstance(r.shared_pfam, str) and r.shared_pfam: return "shares Pfam domain with true-term positives"
    return "no detectable homology to true-term positives"
st["homology_status"] = st.apply(tier, axis=1)
st["callout"] = st.homology_status.str.startswith("no detectable")
st["gene"] = st.gene.fillna(st.pname.str.split(" \\(").str[0])
cols = {"ont": "ontology", "protein_id": "protein_id", "gene": "gene", "pname": "protein_name", "organism": "organism",
        "length": "length", "go_term": "go_term", "name": "go_name", "ev": "evidence", "ref": "reference",
        "date": "annotation_date", "llr": "cav_llr", "cav_rank_ont_filt_pf": "cav_rank",
        "cav_margin": "cav_margin_over_runnerup", "tool_predicted": "deepgose_predicted",
        "tool_rank_ont_filt_pf": "deepgose_rank", "blast_predicted": "diamond_hit_to_term",
        "nn_hit": "diamond_nearest_training_protein", "nn_pident": "diamond_nearest_pident",
        "jack_predicted": "jackhmmer_hit_to_term_e1e-3", "jack_rank_mid": "jackhmmer_rank_mid",
        "jack_bestE_any": "jackhmmer_best_E_to_term", "q_pfam": "query_pfam",
        "shared_pfam": "pfam_shared_with_term_positives_xcount", "n_pos": "n_term_training_positives",
        "n_terms_ont": "n_candidate_terms", "homology_status": "homology_status", "callout": "callout"}
m = st[list(cols)].rename(columns=cols)
m.loc[~m["jackhmmer_hit_to_term_e1e-3"].astype(bool), "jackhmmer_rank_mid"] = np.nan
m.loc[~m.deepgose_predicted.astype(bool), "deepgose_rank"] = np.nan
m = m.sort_values(["ontology", "callout", "cav_margin_over_runnerup"], ascending=[True, False, False])
m.to_csv(HERE / "strict_callouts.csv", index=False)
m[m.ontology == "mf"].drop(columns="ontology").to_csv(HERE / "mf_callouts.csv", index=False)
rows = []
for r in st[st.ont == "mf"].itertuples():
    for k, t in enumerate(r.cav_top6.split("; "), 1):
        n, s = t.rsplit(" (", 1)
        rows.append(dict(protein_id=r.protein_id, gene=r.gene, true_go_term=r.go_term, cav_rank=k, go_name=n,
                         cav_score=float(s.rstrip(")")), is_true_term=(n == r.name)))
pd.DataFrame(rows).to_csv(HERE / "mf_callouts_cav_runnerups.csv", index=False)

# calibration against the stored old MF callouts (jackhmmer columns depend only on hits + members)
oldm = pd.read_csv(BB / "callouts" / "mf_callouts.csv")
chk = oldm[["protein_id", "go_term", "jackhmmer_hit_to_term_e1e-3", "jackhmmer_best_E_to_term"]].merge(
    best.rename(columns={"jack_bestE_any": "newE"}), on=["protein_id", "go_term"], how="left")
hit_e3 = jh[jh.evalue <= 1e-3].merge(members, on="s").groupby(["q", "go_term"]).size().rename("n").reset_index()
chk = chk.merge(hit_e3.rename(columns={"q": "protein_id"}), on=["protein_id", "go_term"], how="left")
chk["newhit"] = chk.n.fillna(0) > 0
cal.append(f"old MF callouts: jackhmmer hit flag reproduced {int((chk.newhit == chk['jackhmmer_hit_to_term_e1e-3']).sum())}/{len(chk)}; "
           f"best-E reproduced {int(np.isclose(chk.newE, chk.jackhmmer_best_E_to_term, rtol=1e-3, equal_nan=True).sum())}/{len(chk)}")

# 6. cav_only tables ---------------------------------------------------------------------------
div = pd.concat([pd.read_csv(GE / "diversity" / f"term_diversity_{o}.csv") for o in C.ONTS])
co = cav_only.merge(ann, on=["protein_id", "go_term"], how="left").merge(u, on="protein_id", how="left") \
    .merge(extra, on=["protein_id", "go_term"], how="left") \
    .merge(div[["go_term", "families_per_100", "largest_family_frac"]], on="go_term", how="left")
for o in C.ONTS:
    co[co.ont == o].to_csv(HERE / f"{o}_cav_only_rank1_cases.csv", index=False)
(HERE / "calibration.txt").write_text("\n".join(cal) + "\n")
print("\n".join(cal))
print("strict callouts:", m.groupby("ontology").callout.agg(["size", "sum"]).to_dict())
print("cav_only:", co.ont.value_counts().to_dict())
