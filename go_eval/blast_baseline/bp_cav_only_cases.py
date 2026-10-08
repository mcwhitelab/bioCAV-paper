"""Top competing terms (filt_pf candidate set) for CAV / DeepGo-SE / BLAST on the CAV-only BP rank-1 pairs."""
import sys, numpy as np, pandas as pd, joblib
from pathlib import Path
GE=Path("/groups/clairemcwhite/cav_workspace/bioCAV-paper/go_eval")
sys.path.insert(0,str(GE/"rank_all_terms"))
import compute_ranks_all_terms as C
SP=Path(sys.argv[1]); ont="bp"
x=pd.read_csv(SP/"bp_cav_only_annot.csv")
ns,parents=C.parse_obo(C.OBO)
names={}
cur=None
for l in open(C.OBO):
    if l.startswith("id: GO:"): cur=l[4:].strip()
    elif l.startswith("name: ") and cur: names[cur]=l[6:].strip(); cur=None
gpad=pd.concat([pd.read_csv(C.GO_EVAL/"static"/"go_annotations"/f"{o}.gpad",sep="\t",comment="!",header=None,usecols=[0,1,3],names=["db","pid","go"]) for o in C.ONTS])
gpad=gpad[gpad.db=="UniProtKB"]; gpad["pid"]=gpad.pid.astype(str).str.split("-").str[0]
gpad=gpad[gpad.pid.isin(set(x.protein_id))]
known={}
for pid,ts in gpad.groupby("pid")["go"]:
    s=set(ts)
    for t in list(s): s|=C.ancestors(t,parents)
    known[pid]=s
ids,M,_=C.load_cavs(); col={g:i for i,g in enumerate(ids)}
cns=np.array([ns.get(g) for g in ids],dtype=object); om=cns==ont
pf=C.load_profab_positives(set(x.protein_id))
embs,sids=C.load_sequence_embeddings(GE/"results"/f"temporal_eval_{ont}"/"eval_embeddings"/"val_proteins.span.fasta.pkl")
idx={}
for i,s in enumerate(sids): idx.setdefault(C.normalise_pid(s),i)
S=joblib.load(C.SCALER).transform(embs)@M.T
tool=C.load_tool(ont); tool=tool[tool.go_term.isin(col)&tool.protein_id.isin(set(x.protein_id))]
# BLAST
mem=pd.read_csv(GE/"blast_baseline"/"term_members.tsv",sep="\t",names=["go_term","sseqid"])
h=pd.read_csv(GE/"blast_baseline"/"val_vs_train.tsv",sep="\t",names=["q","s","pident","length","qlen","slen","evalue","bits"])
h=h[h.q.isin(set(x.protein_id))&(h.q!=h.s)]
bq=h.groupby(["q","s"],as_index=False).bits.max().merge(mem,left_on="s",right_on="sseqid").groupby(["q","go_term"],as_index=False).bits.max()
nn=h.sort_values("bits",ascending=False).drop_duplicates("q").set_index("q")
def top(scores,keep,k=3):
    o=np.argsort(-np.where(keep,scores,-np.inf))[:k]
    return "; ".join(f"{names.get(ids[i],ids[i])} ({scores[i]:.2g})" for i in o if keep[i] and scores[i]>0) or "-"
rows=[]
for r in x.itertuples(index=False):
    p,g=r.protein_id,r.go_term; i,j=idx[p],col[g]
    drop=(C.ancestors(g,parents)|known.get(p,set()))
    for t in pf.get(p,()): drop|={t}|C.ancestors(t,parents)
    drop-={g}
    keep=om.copy(); keep[[col[t] for t in drop if t in col]]=False
    kc=keep.copy(); kc[j]=False
    T=np.zeros(len(ids)); 
    for _,rr in tool[tool.protein_id==p].iterrows(): T[col[rr.go_term]]=rr.tool_score
    B=np.zeros(len(ids))
    for _,rr in bq[bq.q==p].iterrows():
        if rr.go_term in col: B[col[rr.go_term]]=rr.bits
    rows.append(dict(protein_id=p,go_term=g,
        cav_true=S[i,j],cav_next=top(S[i],kc,2),
        dg_true=T[j],dg_top=top(T,kc,3),
        blast_true=B[j],blast_top=top(B,kc,3),
        nn_hit=nn.s.get(p,""),nn_pident=nn.pident.get(p,np.nan),nn_bits=nn.bits.get(p,np.nan),
        known_bp="; ".join(sorted(names.get(t,t) for t in set(gpad[gpad.pid==p].go) if ns.get(t)=="bp"))[:300]))
out=x.merge(pd.DataFrame(rows),on=["protein_id","go_term"])
out.to_csv(SP/"bp_cav_only_competitors.csv",index=False)
print("done",len(out))
