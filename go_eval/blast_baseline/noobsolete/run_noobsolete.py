"""Rerun Fig 2C ranks (CAV/DeepGo-SE and BLAST) with obsolete GO terms removed from
the candidate pool. Obsolete terms get namespace None, so every ont_mask excludes them.
No true term is obsolete. Writes into this directory (the scratchpad), not staging."""
import os, sys
from pathlib import Path
OUT = Path(__file__).resolve().parent
os.environ["FIG_DATA_DIR"] = str(OUT)
GE = Path("/groups/clairemcwhite/cav_workspace/bioCAV-paper/go_eval")
sys.path.insert(0, str(GE / "rank_all_terms")); sys.path.insert(0, str(GE / "blast_baseline"))
import compute_ranks_all_terms as C

_parse = C.parse_obo
def parse_obo_noobs(path):
    ns, parents = _parse(path)
    obs, cur = set(), None
    for l in open(path):
        if l.startswith("id: GO:"): cur = l[4:].strip()
        elif l.startswith("is_obsolete: true"): obs.add(cur)
    for t in obs: ns[t] = None
    print(f"masked {len(obs)} obsolete terms", flush=True)
    return ns, parents
C.parse_obo = parse_obo_noobs

step = sys.argv[1]
if step == "cav":
    C.main()
else:
    import blast_ranks as B
    B.parse_obo = parse_obo_noobs
    B.HERE, B.STAGE = OUT, OUT
    os.environ["HITS_FILE"] = str(GE / "blast_baseline" / "val_vs_train.tsv")
    B.main()
