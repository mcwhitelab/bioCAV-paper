#!/usr/bin/env python
"""Homology tier for the untested discordant pairs, same rules as ../../callouts_C1e-4/03_build_tables.py:
jackhmmer hit (E <= 1e-3) to any training positive of the true term > Pfam domain shared with the
true term's training positives (UniProt xrefs) > no detectable homology. Merges the result into a
full copy of ../cav_rank1_blast_low.csv (cav_rank1_blast_low_jack.csv) and checks the calibration
queries against the earlier run."""
import io, time, urllib.request
from pathlib import Path
import numpy as np, pandas as pd
HERE = Path(__file__).resolve().parent; BB = HERE.parents[1]
HITS = ["q", "s", "pident", "length", "qlen", "slen", "evalue", "bits"]
members = pd.read_csv(BB / "C1e-4_extra" / "term_members.tsv", sep="\t", names=["go_term", "s"])
pairs = pd.read_csv(HERE / "untested_pairs.csv")
jh = pd.read_csv(HERE / "jack_vs_train.tsv", sep="\t", names=HITS)
log = []

# calibration: the three already-tested queries must reproduce the earlier run
old = pd.read_csv(BB / "callouts_C1e-4" / "jack_vs_train.tsv", sep="\t", names=HITS)
cal = sorted(set(jh.q) - set(pairs.protein_id))
a = jh[jh.q.isin(cal)].sort_values(["q", "s"]).reset_index(drop=True)
b = old[old.q.isin(cal)].sort_values(["q", "s"]).reset_index(drop=True)
same = len(a) == len(b) and (a[["q", "s"]].values == b[["q", "s"]].values).all() and np.allclose(a.evalue, b.evalue, rtol=1e-6)
log.append(f"calibration queries {cal}: hits identical to callouts_C1e-4 run: {same} ({len(a)} vs {len(b)})")

hm = jh.merge(members, on="s")
best = hm.groupby(["q", "go_term"]).evalue.min().rename("jackhmmer_best_E_to_term").reset_index() \
    .rename(columns={"q": "protein_id"})
p = pairs.merge(best, on=["protein_id", "go_term"], how="left")
p["jack_hit"] = p.jackhmmer_best_E_to_term <= 1e-3

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
            raise RuntimeError(f"UniProt failed for batch {i}")
    return pd.concat(out)

need = p[~p.jack_hit]
mem_need = members[members.go_term.isin(set(need.go_term))]
pf = uniprot(set(mem_need.s) | set(need.protein_id), "accession,xref_pfam"); pf.columns = ["acc", "pfam"]
P = {r.acc: set(str(r.pfam).strip(";").split(";")) - {"nan", ""} for r in pf.itertuples()}
def shared(r):
    if r.jack_hit: return np.nan
    cnt = {}
    for s in mem_need[mem_need.go_term == r.go_term].s:
        for dom in P.get(r.protein_id, set()) & P.get(s, set()): cnt[dom] = cnt.get(dom, 0) + 1
    return ";".join(f"{k}x{v}" for k, v in sorted(cnt.items()))
p["pfam_shared_with_term_positives_xcount"] = p.apply(shared, axis=1)
p["query_pfam"] = p.protein_id.map(lambda a: ";".join(sorted(P.get(a, set()))))
def tier(r):
    if r.jack_hit: return "remote homology: jackhmmer finds true-term positives"
    if isinstance(r.pfam_shared_with_term_positives_xcount, str) and r.pfam_shared_with_term_positives_xcount:
        return "shares Pfam domain with true-term positives"
    return "no detectable homology to true-term positives"
p["homology_status_new"] = p.apply(tier, axis=1)

d = pd.read_csv(HERE.parent / "cav_rank1_blast_low.csv")
d = d.merge(p[["protein_id", "go_term", "jack_hit", "jackhmmer_best_E_to_term", "query_pfam",
               "pfam_shared_with_term_positives_xcount", "homology_status_new"]], on=["protein_id", "go_term"], how="left")
d["homology_source"] = np.where(d.homology_status.notna(), "callouts_C1e-4 strict run", "this run (jackhmmer_untested)")
d["homology_status"] = d.homology_status.fillna(d.homology_status_new)
d["jackhmmer_hit_to_term_e1e-3"] = d["jackhmmer_hit_to_term_e1e-3"].where(d["jackhmmer_hit_to_term_e1e-3"].notna(), d.jack_hit)
d = d.drop(columns=["jack_hit", "homology_status_new"])
d.to_csv(HERE / "cav_rank1_blast_low_jack.csv", index=False)
clean = d[d.homology_status.str.startswith("no detectable")].sort_values(["ontology", "llr"], ascending=[True, False])
clean.to_csv(HERE / "cav_rank1_blast_low_no_homology_all.csv", index=False)
log.append("homology status, newly tested pairs:\n" + p.groupby("ontology").homology_status_new.value_counts().to_string())
log.append("homology status, all discordant pairs:\n" + d.groupby("ontology").homology_status.value_counts().to_string())
(HERE / "summary.txt").write_text("\n\n".join(log) + "\n"); print("\n\n".join(log))
