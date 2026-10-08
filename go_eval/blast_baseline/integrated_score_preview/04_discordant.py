#!/usr/bin/env python
"""CAV rank 1 (LLR > 0, i.e. the Fig 2 rank panel's "Rank 1" tier) but BLAST rank > 10 or
no DIAMOND hit to the true term's training positives. From integrated_ranks_{ont}.csv
(identical to live figure_data/go_rank_families_{ont}.csv ranks). Protein metadata from the
UniProt REST API (as callouts_C1e-4/03_build_tables.py); jackhmmer / Pfam homology status
joined from callouts_C1e-4/strict_callouts.csv where the pair is in that set.
Sorted: no BLAST hit first, then worst BLAST rank, then highest LLR."""
import io, time, urllib.request
from pathlib import Path
import pandas as pd
HERE = Path(__file__).resolve().parent
GE = HERE.parents[1]
names, cur = {}, None
for l in open(GE / "static" / "go_annotations" / "go.obo"):
    l = l.strip()
    if l.startswith("id: GO:"): cur = l[4:]
    elif l.startswith("name: ") and cur: names[cur] = l[6:]; cur = None
rows = []
for ont in ["mf", "bp", "cc"]:
    d = pd.read_csv(HERE / f"integrated_ranks_{ont}.csv")
    d = d[(d.cav_rank == 1) & (d.llr > 0) & (~d.blast_predicted | (d.blast_rank_mid > 10))].copy()
    d.insert(0, "ontology", ont)
    rows.append(d)
d = pd.concat(rows)
d["go_name"] = d.go_term.map(names)
def uniprot(ids, fields):
    out, ids = [], sorted(ids)
    for i in range(0, len(ids), 100):
        url = ("https://rest.uniprot.org/uniprotkb/accessions?accessions=" + ",".join(ids[i:i + 100])
               + f"&fields={fields}&format=tsv&size=100")
        for t in range(4):
            try:
                out.append(pd.read_csv(io.StringIO(urllib.request.urlopen(url, timeout=90).read().decode()), sep="\t")); break
            except Exception as e:  # noqa: BLE001
                print("retry", e, flush=True); time.sleep(5)
        else:
            raise RuntimeError(i)
    return pd.concat(out)
u = uniprot(set(d.protein_id), "accession,gene_primary,protein_name,organism_name,length")
u.columns = ["protein_id", "gene", "protein_name", "organism", "length"]
d = d.merge(u.drop_duplicates("protein_id"), on="protein_id", how="left")
s = pd.read_csv(HERE.parent / "callouts_C1e-4" / "strict_callouts.csv")[
    ["protein_id", "go_term", "jackhmmer_hit_to_term_e1e-3", "homology_status"]]
d = d.merge(s, on=["protein_id", "go_term"], how="left")
d["deepgose_rank"] = d.tool_rank_ont_filt_pf.where(d.tool_predicted.astype(bool))
d["blast_rank"] = d.blast_rank_mid.where(d.blast_predicted, other=float("nan"))
d["blast_class"] = d.blast_predicted.map({False: "no hit", True: "rank >10"})
d = d.sort_values(["ontology", "blast_predicted", "blast_rank", "llr"],
                  ascending=[True, True, False, False], na_position="first")
cols = ["ontology", "protein_id", "gene", "protein_name", "organism", "length", "go_term", "go_name",
        "n_families", "cav_rank", "llr", "blast_class", "blast_rank", "blast_true_bits", "blast_true_pident",
        "blast_n_hit_terms", "deepgose_rank", "rankprod_rank", "soft_eps0.1_rank",
        "jackhmmer_hit_to_term_e1e-3", "homology_status", "n_cand"]
d[cols].to_csv(HERE / "cav_rank1_blast_low.csv", index=False)
tot = {o: len(pd.read_csv(HERE / f"integrated_ranks_{o}.csv").query("cav_rank == 1 and llr > 0")) for o in ["mf","bp","cc"]}
print(d.groupby(["ontology", "blast_class"]).size().unstack(fill_value=0).assign(cav_rank1_total=pd.Series(tot)))
fam = pd.concat([pd.read_csv(HERE / f"integrated_ranks_{o}.csv").assign(ontology=o) for o in ["mf","bp","cc"]])
print("median n_families: discordant", d.groupby("ontology").n_families.median().to_dict(),
      "| all CAV rank-1", fam.query("cav_rank==1 and llr>0").groupby("ontology").n_families.median().to_dict(),
      "| all pairs", fam.groupby("ontology").n_families.median().to_dict())
