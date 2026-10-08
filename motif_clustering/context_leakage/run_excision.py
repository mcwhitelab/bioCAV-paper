#!/usr/bin/env python3
"""
Does whole-protein context leak into domain-window CAVs?

The 27k motif library embeds the FULL protein through ESM-C and pools only the
domain window (esm2_tcav/src/embed.py: process_dataset feeds sample['sequence'];
pool_window_embeddings slices [start+1:end+1]).  Attention therefore sees the
whole chain.  If that context leaks, CAVs for two domains that co-occur in the
same protein should look alike even when the domains are structurally unrelated.

Observation being tested: pre-ATP-grasp families (Pfam clan CL0063, Rossmann)
sit adjacent to ATP-grasp families (CL0179) in the CAV dendrogram, away from the
rest of their own clan.

Two matched regimes, identical windows and identical pooling extent.  The only
difference is what ESM saw:

  FULL     positives: full protein -> pool the domain window   [current pipeline]
           negatives: full protein -> pool a random window of matched length
  EXCISED  positives: domain substring alone -> pool all
           negatives: the SAME random substring alone -> pool all

Collapse of cos(CAV_pre, CAV_grasp) from FULL to EXCISED == context leakage.
"""

import json, argparse, logging, random
from pathlib import Path

import numpy as np
import torch
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold, cross_val_score

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)

TCAV_DATA = Path("/groups/clairemcwhite/ahmad_workspace/esm_c/tcav_data")
NEG_JSONL = Path("/groups/clairemcwhite/ahmad_workspace/esm_c/neg_data/neg_10000.jsonl")
CKPT = "/groups/clairemcwhite/models/ESMplusplus_large"
LAYER = 25
MAXLEN = 1000          # cap so the full-context forward never truncates a window

# (pre-ATP-grasp family, ATP-grasp family, label).  Cognate = the two domains are
# the N- and C-terminal halves of one protein.
COGNATE = [
    ("PF01820", "PF07478", "Dala_Dala_lig_N | Dala_Dala_lig_C"),
    ("PF18030", "PF08443", "Rimk_N | RimK"),
    ("PF00289", "PF02786", "Biotin_carb_N | CPSase_L_D2"),
    ("PF02844", "PF01071", "GARS_N | GARS_A"),
    ("PF17927", "PF05770", "Ins134_P3_kin_N | Ins134_P3_kin"),
    ("PF02078", "PF02750", "Synapsin | Synapsin_C"),
]
# Controls: both members are ATP-grasp (CL0179) but live in different proteins,
# so context leakage predicts no change between regimes.
CONTROL = [
    ("PF02222", "PF02655", "ATP-grasp | ATP-grasp_3"),
    ("PF01071", "PF02955", "GARS_A | GSH-S_ATP"),
    ("PF07478", "PF08443", "Dala_Dala_lig_C | RimK"),
]

FAMILIES = sorted({f for p in COGNATE + CONTROL for f in p[:2]})


def load_positives(pfam, n, rng):
    recs = []
    for line in open(TCAV_DATA / pfam / "pos_1000.jsonl"):
        d = json.loads(line)
        s, e = d["domain_span_1based_inclusive"]
        # union across all copies, matching batch_train_cavs._domain_span_to_halfopen
        spans = d.get("all_domain_spans_1based_inclusive") or [[s, e]]
        st, en = max(0, min(x[0] for x in spans) - 1), max(x[1] for x in spans)
        if d["sequence_length"] > MAXLEN or en - st < 20:
            continue
        recs.append({"acc": d["accession"], "seq": d["sequence"], "span": (st, en)})
    rng.shuffle(recs)
    return recs[:n]


def load_negatives(n, rng):
    recs = []
    for line in open(NEG_JSONL):
        d = json.loads(line)
        seq = d["sequence"]
        if len(seq) > MAXLEN or len(seq) < 250:
            continue
        recs.append({"acc": d["accession"], "seq": seq})
    rng.shuffle(recs)
    return recs[:n]


def negative_windows(negs, pos_lengths, seed):
    """Random window per negative, length drawn from this family's domain lengths."""
    r = random.Random(seed)
    out = []
    for rec in negs:
        L = min(r.choice(pos_lengths), len(rec["seq"]))
        st = r.randrange(0, len(rec["seq"]) - L + 1)
        out.append((st, st + L))
    return out


class Embedder:
    def __init__(self, device):
        from transformers import AutoModelForMaskedLM
        log.info("loading %s", CKPT)
        self.model = AutoModelForMaskedLM.from_pretrained(CKPT, trust_remote_code=True)
        self.model = self.model.to(device).eval()
        self.tok = self.model.tokenizer
        self.device = device

    @torch.no_grad()
    def pooled(self, seqs, spans, batch_size=16):
        """spans are 0-based half-open in sequence coords; None == pool everything."""
        out = []
        for i in range(0, len(seqs), batch_size):
            bs, bsp = seqs[i:i + batch_size], spans[i:i + batch_size]
            enc = self.tok(bs, padding="longest", truncation=True,
                           max_length=MAXLEN + 2, return_tensors="pt",
                           add_special_tokens=True)
            ids = enc["input_ids"].to(self.device)
            mask = enc["attention_mask"].to(self.device)
            hs = self.model(input_ids=ids, attention_mask=mask,
                            output_hidden_states=True).hidden_states[LAYER]
            hs = hs.float().cpu().numpy()
            am = mask.cpu().numpy()
            for j, sp in enumerate(bsp):
                if sp is None:
                    n = int(am[j].sum())
                    out.append(hs[j, 1:n - 1].mean(axis=0))   # strip BOS/EOS
                else:
                    a, b = sp[0] + 1, sp[1] + 1               # BOS offset, as in pipeline
                    out.append(hs[j, a:b].mean(axis=0))
            del hs, ids, mask
            if self.device.startswith("cuda"):
                torch.cuda.empty_cache()
            if (i // batch_size) % 25 == 0:
                log.info("  %d/%d", i, len(seqs))
        return np.stack(out)


def fit_cav(Xpos, Xneg):
    """Exactly the library's recipe: StandardScaler, no PCA, balanced L2 logistic."""
    X = np.vstack([Xpos, Xneg])
    y = np.hstack([np.ones(len(Xpos)), np.zeros(len(Xneg))])
    Xs = StandardScaler().fit_transform(X)
    clf = LogisticRegression(max_iter=2000, C=1.0, class_weight="balanced")
    clf.fit(Xs, y)
    auc = cross_val_score(
        LogisticRegression(max_iter=2000, C=1.0, class_weight="balanced"),
        Xs, y, cv=StratifiedKFold(5, shuffle=True, random_state=42), scoring="roc_auc"
    ).mean()
    v = clf.coef_[0]
    return v / (np.linalg.norm(v) + 1e-10), float(auc)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-pos", type=int, default=300)
    ap.add_argument("--n-neg", type=int, default=1000)
    ap.add_argument("--out", default="results")
    ap.add_argument("--batch-size", type=int, default=16)
    args = ap.parse_args()

    outdir = Path(args.out); outdir.mkdir(parents=True, exist_ok=True)
    device = "cuda" if torch.cuda.is_available() else "cpu"
    log.info("device=%s families=%d", device, len(FAMILIES))

    rng = random.Random(0)
    negs = load_negatives(args.n_neg, rng)
    log.info("negatives: %d", len(negs))

    emb = Embedder(device)
    cavs, aucs = {}, {}

    for k, pfam in enumerate(FAMILIES):
        pos = load_positives(pfam, args.n_pos, random.Random(1000 + k))
        lens = [sp[1] - sp[0] for sp in (p["span"] for p in pos)]
        nwin = negative_windows(negs, lens, seed=2000 + k)
        log.info("[%s] %d positives, median domain len %d", pfam, len(pos), int(np.median(lens)))

        # FULL: whole protein through the model, pool the window
        log.info("[%s] FULL positives", pfam)
        P_full = emb.pooled([p["seq"] for p in pos], [p["span"] for p in pos], args.batch_size)
        log.info("[%s] FULL negatives", pfam)
        N_full = emb.pooled([n["seq"] for n in negs], nwin, args.batch_size)

        # EXCISED: the same windows, cut out and run alone
        log.info("[%s] EXCISED positives", pfam)
        P_exc = emb.pooled([p["seq"][p["span"][0]:p["span"][1]] for p in pos],
                           [None] * len(pos), args.batch_size)
        log.info("[%s] EXCISED negatives", pfam)
        N_exc = emb.pooled([n["seq"][w[0]:w[1]] for n, w in zip(negs, nwin)],
                           [None] * len(negs), args.batch_size)

        for regime, (P, N) in {"full": (P_full, N_full), "excised": (P_exc, N_exc)}.items():
            v, auc = fit_cav(P, N)
            cavs[(regime, pfam)] = v
            aucs[(regime, pfam)] = auc
            log.info("[%s] %s CAV cv-AUROC=%.4f", pfam, regime, auc)

        np.save(outdir / f"{pfam}_cav_full.npy", cavs[("full", pfam)])
        np.save(outdir / f"{pfam}_cav_excised.npy", cavs[("excised", pfam)])

    def cos(r, a, b):
        return float(cavs[(r, a)] @ cavs[(r, b)])

    rows = []
    for group, pairs in (("cognate", COGNATE), ("control", CONTROL)):
        for a, b, label in pairs:
            rows.append({"group": group, "label": label, "a": a, "b": b,
                         "cos_full": cos("full", a, b), "cos_excised": cos("excised", a, b),
                         "auc_full_a": aucs[("full", a)], "auc_full_b": aucs[("full", b)],
                         "auc_exc_a": aucs[("excised", a)], "auc_exc_b": aucs[("excised", b)]})

    # background: all unordered family pairs, for scale
    bg = {r: [cos(r, FAMILIES[i], FAMILIES[j])
              for i in range(len(FAMILIES)) for j in range(i + 1, len(FAMILIES))]
          for r in ("full", "excised")}

    import pandas as pd
    df = pd.DataFrame(rows)
    df["delta"] = df.cos_excised - df.cos_full
    df.to_csv(outdir / "pair_cosines.csv", index=False)
    json.dump({r: {"mean": float(np.mean(v)), "sd": float(np.std(v))} for r, v in bg.items()},
              open(outdir / "background.json", "w"), indent=2)

    print("\n" + df.to_string(index=False))
    print("\nbackground cos over all %d family pairs:" % len(bg["full"]))
    for r, v in bg.items():
        print(f"  {r:8s} mean={np.mean(v):+.4f} sd={np.std(v):.4f}")
    print("\nmean delta  cognate: %+.4f   control: %+.4f"
          % (df[df.group == "cognate"].delta.mean(), df[df.group == "control"].delta.mean()))


if __name__ == "__main__":
    main()
