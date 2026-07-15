#!/usr/bin/env python3
"""
Score 3 "big receptor protein" candidates (MET, LDLR, Plexin A1) for a
higher-domain-count second Figure 4 worked example, against layer-25
("display layer 26") CAVs from the 20k-motif library.
"""
import json
import sys
from pathlib import Path

import numpy as np
import torch
from joblib import load as joblib_load
from transformers import AutoConfig, AutoModelForMaskedLM

REPO_DIR = Path("/groups/clairemcwhite/ahmad_workspace/esm_c")
sys.path.insert(0, str(REPO_DIR))
from run_test_detection_full_embed import generate_windows, resolve_stride, cumsum_mean_pool  # noqa: E402

CAV_DIR = Path("/xdisk/clairemcwhite/shamail/tcav_outputs_esmplusplus_all/cavs")
MODEL_PATH = "/groups/clairemcwhite/models/ESMplusplus_large"
DATA_DIR = Path("/groups/clairemcwhite/cav_workspace/bioCAV-paper/figures/figure_data/vav_motif_repro/extra_proteins")
OUT_JSON = DATA_DIR / "candidates_bigreceptors_scores.json"

HS_INDEX = 25
BOS_OFFSET = 1

WINDOW_LENGTHS = {
    "PF01403": 173, "PF01437": 46, "PF01833": 80, "PF07714": 254,
    "PF00057": 37, "PF00058": 42, "PF07645": 40, "PF14670": 36,
    "PF17960": 90, "PF18020": 94, "PF24479": 27, "PF08337": 249, "PF20170": 113,
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


def find_peaks_per_interval(curve, gt_intervals, window_len, seq_len):
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
    records = json.loads((DATA_DIR / "candidates_bigreceptors.json").read_text())
    model, tokenizer, device = load_model()

    results = {}
    for acc, rec in records.items():
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

        protein_result = {"seq_len": seq_len, "name": rec["name"], "domains": {}}
        for motif_id in gt_by_domain:
            cav_file = CAV_DIR / motif_id / "L25_concept_v1.npy"
            scaler_file = CAV_DIR / motif_id / "L25_scaler_v1.pkl"
            if not cav_file.exists():
                print(f"[WARN] no CAV for {motif_id}", file=sys.stderr)
                continue
            cav = np.load(cav_file)
            scaler = joblib_load(scaler_file)
            window_len = WINDOW_LENGTHS[motif_id]

            curve = score_domain(hs_L, seq_len, window_len, cav, scaler)
            if curve is None:
                continue

            gt_intervals = [(s, e) for s, e, _ in gt_by_domain[motif_id]]
            peaks = find_peaks_per_interval(curve, gt_intervals, window_len, seq_len)

            protein_result["domains"][motif_id] = {
                "window_length": window_len,
                "ground_truth": [{"start": s, "end": e, "name": n} for s, e, n in gt_by_domain[motif_id]],
                "curve": [None if np.isnan(v) else round(float(v), 5) for v in curve],
                "peaks": peaks,
                "mode": "per_interval",
            }
        print(f"{acc}: {len(protein_result['domains'])} domain types scored", file=sys.stderr)
        results[acc] = protein_result

    OUT_JSON.write_text(json.dumps(results))
    print(f"Wrote {OUT_JSON}", file=sys.stderr)


if __name__ == "__main__":
    main()
