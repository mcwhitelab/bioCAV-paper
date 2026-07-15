#!/usr/bin/env python3
"""
Score 9 additional example proteins (arxiv 2511.21614v1 Appendix A grid, plus
the Figure-3 repeat-domain examples Q9Y7R4/SET1_SCHPO and Q54SA1/PLDZ_DICDI)
against layer-25 ("display layer 26") CAVs from the 20k-motif library at
/xdisk/clairemcwhite/shamail/tcav_outputs_esmplusplus_all/cavs (per user:
"use the new cavs from the 20k").

Two peak-finding modes per (protein, domain):
  - "per_interval" (default, used for the Appendix-A grid): one peak per
    ground-truth interval -- the local max within [start, end] (padded by
    half a window on each side).
  - "top2" (used for PF00076 on Q9Y7R4 and PF00614 on Q54SA1, the two
    Figure-3 repeat-domain showcases): find the two strongest global local
    maxima (min separation = window_length), regardless of how many
    ground-truth intervals exist, and label each as "annotated" (overlaps
    a known interval) or "candidate" (does not) -- reproducing the paper's
    "second peak not annotated in UniProt/InterPro-N" finding organically
    from the actual CAV scores, rather than hardcoding peak positions.
"""
import json
import sys
from pathlib import Path

import numpy as np
import torch
from joblib import load as joblib_load
from scipy.signal import find_peaks
from transformers import AutoConfig, AutoModelForMaskedLM

REPO_DIR = Path("/groups/clairemcwhite/ahmad_workspace/esm_c")
sys.path.insert(0, str(REPO_DIR))
from run_test_detection_full_embed import generate_windows, resolve_stride, cumsum_mean_pool  # noqa: E402

CAV_DIR = Path("/xdisk/clairemcwhite/shamail/tcav_outputs_esmplusplus_all/cavs")
MODEL_PATH = "/groups/clairemcwhite/models/ESMplusplus_large"
DATA_DIR = Path("/groups/clairemcwhite/cav_workspace/bioCAV-paper/figures/figure_data/vav_motif_repro/extra_proteins")
OUT_JSON = DATA_DIR / "extra_proteins_scores.json"

DISPLAY_LAYER = 26  # hs[25] -> display "layer 26" (embedding counted as layer 1)
HS_INDEX = 25
BOS_OFFSET = 1

WINDOW_LENGTHS = {
    "PF00041": 86, "PF00536": 61, "PF00096": 23, "PF00069": 250, "PF00130": 51,
    "PF00168": 100, "PF00002": 241, "PF07645": 40, "PF00249": 47, "PF00013": 61,
    "PF00076": 68, "PF00271": 109, "PF00098": 17, "PF00856": 111, "PF00614": 27,
}

# which domains to plot per protein (subset of ground_truth_annotations,
# matching each source panel's legend), and peak-finding mode per (protein, domain)
PANELS = {
    "P29317": {"domains": ["PF00041", "PF00536"], "modes": {}},
    "P03001": {"domains": ["PF00096"], "modes": {}},
    "Q4R4U2": {"domains": ["PF00069", "PF00130", "PF00168"], "modes": {}},
    "P48960": {"domains": ["PF00002", "PF07645"], "modes": {}},
    "Q9ZQ85": {"domains": ["PF00249"], "modes": {}},
    "Q8CGX0": {"domains": ["PF00013", "PF00076"], "modes": {}},
    "Q4JG17": {"domains": ["PF00098", "PF00271"], "modes": {}},
    "Q9Y7R4": {"domains": ["PF00076", "PF00856"], "modes": {"PF00076": "top2"}},
    "Q54SA1": {"domains": ["PF00614"], "modes": {"PF00614": "top2"}},
}


def load_model():
    print(f"Loading model from {MODEL_PATH}...", file=sys.stderr)
    device = "cuda" if torch.cuda.is_available() else "cpu"
    cfg = AutoConfig.from_pretrained(MODEL_PATH, output_hidden_states=True, trust_remote_code=True)
    model = AutoModelForMaskedLM.from_pretrained(MODEL_PATH, config=cfg, trust_remote_code=True)
    tokenizer = model.tokenizer
    model = model.to(device).eval()
    return model, tokenizer, device


def score_domain(hs_L, seq_len, window_len, cav, scaler):
    stride = resolve_stride(window_len, 1)
    spans = list(generate_windows(seq_len, window_len, stride))
    if not spans:
        return None
    pooled = cumsum_mean_pool(hs_L, spans, bos_offset=BOS_OFFSET)
    X = pooled.detach().cpu().numpy().astype(np.float32)
    Xs = scaler.transform(X)
    window_scores = Xs @ cav

    residue_sum = np.zeros(seq_len, dtype=np.float64)
    residue_cnt = np.zeros(seq_len, dtype=np.float64)
    for (s, e), sc in zip(spans, window_scores):
        residue_sum[s:e] += sc
        residue_cnt[s:e] += 1
    residue_cnt[residue_cnt == 0] = np.nan
    return residue_sum / residue_cnt


def find_domain_peaks(curve, mode, gt_intervals, window_len, seq_len):
    """Returns list of {position (1-based), score, annotated (bool)}."""
    valid = ~np.isnan(curve)
    if mode == "top2":
        peaks_idx, _ = find_peaks(np.nan_to_num(curve, nan=-1e9), distance=max(window_len, 5))
        if len(peaks_idx) == 0:
            return []
        order = np.argsort(curve[peaks_idx])[::-1]
        top = peaks_idx[order][:2]
        out = []
        for p in top:
            annotated = any(s <= p + 1 <= e for s, e in gt_intervals)
            out.append({"position": int(p) + 1, "score": float(curve[p]), "annotated": bool(annotated)})
        return sorted(out, key=lambda d: d["position"])
    else:
        out = []
        pad = window_len // 2
        for s, e in gt_intervals:
            lo = max(0, s - 1 - pad)
            hi = min(seq_len, e + pad)
            seg = curve[lo:hi]
            if np.all(np.isnan(seg)):
                continue
            local_max = np.nanargmax(seg)
            out.append({"position": int(lo + local_max) + 1, "score": float(seg[local_max]), "annotated": True})
        return out


def main():
    records = json.loads((DATA_DIR / "records.json").read_text())
    model, tokenizer, device = load_model()

    results = {}
    for acc, panel in PANELS.items():
        rec = records[acc]
        seq = rec["sequence"]
        seq_len = len(seq)
        gt_by_domain = {}
        for a in rec["ground_truth_annotations"]:
            gt_by_domain.setdefault(a["motif_id"], []).append((a["start"], a["end"], a["name"]))

        enc = tokenizer(seq, return_tensors="pt", padding=False, truncation=False, add_special_tokens=True)
        input_ids = enc["input_ids"].to(device)
        attn = enc.get("attention_mask")
        if attn is not None:
            attn = attn.to(device)
        with torch.no_grad():
            out = model(input_ids=input_ids, attention_mask=attn, output_hidden_states=True, return_dict=True)
        hs_L = out.hidden_states[HS_INDEX][0]

        protein_result = {"seq_len": seq_len, "domains": {}}
        for motif_id in panel["domains"]:
            cav_file = CAV_DIR / motif_id / "L25_concept_v1.npy"
            scaler_file = CAV_DIR / motif_id / "L25_scaler_v1.pkl"
            if not cav_file.exists():
                print(f"[WARN] no CAV for {motif_id}, skipping {acc}/{motif_id}", file=sys.stderr)
                continue
            cav = np.load(cav_file)
            scaler = joblib_load(scaler_file)
            window_len = WINDOW_LENGTHS[motif_id]

            curve = score_domain(hs_L, seq_len, window_len, cav, scaler)
            if curve is None:
                continue

            gt_intervals = [(s, e) for s, e, _ in gt_by_domain.get(motif_id, [])]
            mode = panel["modes"].get(motif_id, "per_interval")
            peaks = find_domain_peaks(curve, mode, gt_intervals, window_len, seq_len)

            protein_result["domains"][motif_id] = {
                "window_length": window_len,
                "ground_truth": [{"start": s, "end": e, "name": n} for s, e, n in gt_by_domain.get(motif_id, [])],
                "curve": [None if np.isnan(v) else round(float(v), 5) for v in curve],
                "peaks": peaks,
                "mode": mode,
            }
            print(f"{acc}/{motif_id}: {len(peaks)} peak(s), mode={mode}", file=sys.stderr)

        results[acc] = protein_result

    OUT_JSON.write_text(json.dumps(results))
    print(f"Wrote {OUT_JSON}", file=sys.stderr)


if __name__ == "__main__":
    main()
