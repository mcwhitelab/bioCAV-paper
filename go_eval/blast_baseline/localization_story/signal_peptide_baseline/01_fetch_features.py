#!/usr/bin/env python
"""UniProt features for every CC-term CAV training positive and every Fig 2B CC query protein:
signal peptide present (ft_signal), number of TM helices (ft_transmem), kingdom (from lineage).
Writes features.csv. UniProt signal features are curated or SignalP-derived (Swiss-Prot / TrEMBL
automatic annotation), i.e. this is a sequence-feature baseline, not a homology one."""
import io, sys, time, urllib.request
from pathlib import Path
import pandas as pd
HERE = Path(__file__).resolve().parent; PAPER = HERE.parents[3]; BB = PAPER / "go_eval/blast_baseline"
sys.path.insert(0, str(PAPER / "go_eval/rank_all_terms"))
import compute_ranks_all_terms as C  # noqa: E402
ns, _ = C.parse_obo(C.OBO)
m = pd.read_csv(BB / "C1e-4_extra/term_members.tsv", sep="\t", names=["go_term", "s"])
m = m[m.go_term.map(ns) == "cc"]
q = pd.read_csv(PAPER / "figures/figure_data/go_rank_families_cc.csv").protein_id
ids = sorted(set(m.s) | set(q))
out_f = HERE / "features_raw.tsv"
done = set(pd.read_csv(out_f, sep="\t").Entry) if out_f.exists() else set()
todo = [i for i in ids if i not in done]
print(len(ids), "proteins;", len(todo), "to fetch", flush=True)
N = 400
for i in range(0, len(todo), N):
    url = ("https://rest.uniprot.org/uniprotkb/accessions?accessions=" + ",".join(todo[i:i + N])
           + f"&fields=accession,ft_signal,ft_transmem,lineage&format=tsv&size={N}")
    for t in range(5):
        try:
            d = pd.read_csv(io.StringIO(urllib.request.urlopen(url, timeout=120).read().decode()), sep="\t"); break
        except Exception as e:  # noqa: BLE001
            print("retry", e, flush=True); time.sleep(10)
    else:
        raise RuntimeError(f"batch {i}")
    d.to_csv(out_f, sep="\t", index=False, mode="a", header=not out_f.exists())
    print(i + N, flush=True)
d = pd.read_csv(out_f, sep="\t").drop_duplicates("Entry")
KING = ["Viridiplantae", "Fungi", "Metazoa", "Bacteria", "Archaea", "Viruses"]
def king(s):
    s = str(s)
    for k in KING:
        if k + " (" in s: return k
    return "Other eukaryote" if "Eukaryota" in s else "Other"
f = pd.DataFrame({"protein_id": d.Entry, "sp": d["Signal peptide"].notna(),
                  "n_tm": d["Transmembrane"].fillna("").str.count("TRANSMEM"), "kingdom": d["Taxonomic lineage"].map(king)})
f.to_csv(HERE / "features.csv", index=False)
print(len(f), "with features;", len(set(ids) - set(f.protein_id)), "missing (obsolete/merged accessions)")
print(f.groupby("kingdom").sp.agg(["size", "mean"]))
