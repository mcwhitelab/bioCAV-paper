#!/usr/bin/env python3
"""
Reproduce Figure 2 (Q9NHV9/VAV_DROME layerwise CAV score profiles) from the
arxiv draft, using CAVs trained fresh on 36 layers for PF00621/PF00130/
PF00017/PF00018 (see train_all_layers.log in this directory) plus Ahmad's
existing training data (read-only) for window lengths.

For each layer with a trained CAV, embeds the query protein once, slides a
window of the motif's median positive-domain length across the sequence at
stride 1, scores each window (scaler + CAV dot product, matching
run_test_detection_full_embed.py::score_motif_from_hidden), then assigns
each residue the mean score of all windows containing it -- exactly the
continuous per-residue curve the paper describes.

Outputs a single JSON with per-layer, per-domain, per-residue score arrays
so the plotting step (R, to match figures.R conventions) doesn't need torch.
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

CAV_DIR = Path("/groups/clairemcwhite/cav_workspace/bioCAV-paper/figures/figure_data/vav_motif_repro/"
                "tcav_outputs_all_layers/cavs")
MODEL_PATH = "/groups/clairemcwhite/models/ESMplusplus_large"
OUT_JSON = Path("/groups/clairemcwhite/cav_workspace/bioCAV-paper/figures/figure_data/vav_motif_repro/"
                 "vav_layerwise_scores.json")

QUERY_SEQ = (
    "MASSSSSNSFGGVAGVNGDLWRECVAWLTRCKVIPPDHKAAQPDAEIRILAMTLRDGVLLCNLVIHLDPSSLDPREFNRKPQMAQFLCSKNIKLFLDVCHNNFGIRDADLFEPTMLYDLT"
    "NFHRVLITLSKLSQCRKVQQLHPDLIGFNLQLSPTERSHSDEAIYKDLHSTTTDNIACNGTGYDHTNTKEEEVYQDLCALHRTSRSQTASSTSFEQRDYVIRELIDTESNYLDVLTALKT"
    "KFMGPLERHLNQDELRLIFPRIRELVDIHTKFLDKLRESLTPNAKVKMAQVFLDFREPFLIYGEFCSLLLGAIDYLADVCKKNQIIDQLVQKCERDYNVGKLQLRDILSVPMQRILKYHL"
    "LLDKLVKETSPLHDDYRSLERAKEAMIDVSQYINEVKRDSDHLVIIQKVKDSICDIHLLQNGNGSDLLQYGRLLLDGELHIKAHEDQKTKLRYAFVFDKILIMVKALHIKTGDMQYTYRD"
    "SHNLADYRVEQSHSRRTIGRDTRFKYQLLLARKSGKTAFTLYLKSEHERDKWRKALTEAMESLEPPGCQSTDHKMEIYTFDAPTTCRHCSKFLKGRIHQGYRCKVCQISVHKGCISSTGR"
    "CKQNPVSVPPPVCDRQLSEFNWFAGNMDRETAAHRLENRRIGTYLLRVRPQGPSTAHETMYALSLKTDDNVIKHMKINQENSGDSMLYCLSSRRHFKTIVELVSYYERNDLGENFAGLNQ"
    "SLQWPYKEVIATALYDYEPKAGSNQLQLRTDCQVLVIGKDGDSKGWWRGKIGDTVGYFPKEYVQEQKLASEEL"
)

GROUND_TRUTH = {
    "PF00621": {"start": 220, "end": 394, "name": "RhoGEF domain"},
    "PF00130": {"start": 553, "end": 597, "name": "C1 domain"},
    "PF00017": {"start": 622, "end": 706, "name": "SH2 domain"},
    "PF00018": {"start": 732, "end": 780, "name": "SH3 domain"},
}

WINDOW_LENGTHS = {"PF00017": 76, "PF00018": 47, "PF00130": 51, "PF00621": 178}
NUM_LAYERS = 36
BOS_OFFSET = 1
DEVICE = "cuda" if torch.cuda.is_available() else "cpu"


def main():
    print(f"Loading model from {MODEL_PATH} on {DEVICE}...", file=sys.stderr)
    cfg = AutoConfig.from_pretrained(MODEL_PATH, output_hidden_states=True, trust_remote_code=True)
    model = AutoModelForMaskedLM.from_pretrained(MODEL_PATH, config=cfg, trust_remote_code=True)
    tokenizer = model.tokenizer
    model = model.to(DEVICE).eval()

    enc = tokenizer(QUERY_SEQ, return_tensors="pt", padding=False, truncation=False, add_special_tokens=True)
    input_ids = enc["input_ids"].to(DEVICE)
    attention_mask = enc.get("attention_mask")
    if attention_mask is not None:
        attention_mask = attention_mask.to(DEVICE)

    with torch.no_grad():
        out = model(input_ids=input_ids, attention_mask=attention_mask,
                     output_hidden_states=True, return_dict=True)
    hs = out.hidden_states  # tuple len NUM_LAYERS+1 (embeddings + each layer)
    print(f"Got {len(hs)} hidden-state tensors (embeddings + {len(hs) - 1} layers)", file=sys.stderr)

    seq_len = len(QUERY_SEQ)
    results = {"seq_len": seq_len, "ground_truth": GROUND_TRUTH, "domains": {}}

    for motif_id, window_len in WINDOW_LENGTHS.items():
        stride = resolve_stride(window_len, 1)
        spans = list(generate_windows(seq_len, window_len, stride))
        if not spans:
            print(f"[WARN] no windows for {motif_id}, skipping", file=sys.stderr)
            continue

        per_layer_curve = {}
        for L in range(0, NUM_LAYERS):  # matches CAV filenames L0..L35 == hs[0..35] directly
            cav_file = CAV_DIR / motif_id / f"L{L}_concept_v1.npy"
            scaler_file = CAV_DIR / motif_id / f"L{L}_scaler_v1.pkl"
            if not cav_file.exists() or not scaler_file.exists():
                continue

            # Display label counts the embedding output as "layer 1" (paper's
            # "indexed from 1" convention) -- so internal L=25 (the CAV already
            # trained before this session, used in the paper's panel B) displays
            # as "Layer 26", matching the paper caption exactly.
            display_layer = L + 1
            hs_L = hs[L][0]  # (T, D)
            pooled = cumsum_mean_pool(hs_L, spans, bos_offset=BOS_OFFSET)
            X = pooled.detach().cpu().numpy().astype(np.float32)

            cav = np.load(cav_file)
            scaler = joblib_load(scaler_file)
            Xs = scaler.transform(X)
            window_scores = Xs @ cav  # (N,)

            # assign each residue the mean score of all windows containing it
            residue_sum = np.zeros(seq_len, dtype=np.float64)
            residue_cnt = np.zeros(seq_len, dtype=np.float64)
            for (s, e), sc in zip(spans, window_scores):
                residue_sum[s:e] += sc
                residue_cnt[s:e] += 1
            residue_cnt[residue_cnt == 0] = np.nan
            residue_curve = residue_sum / residue_cnt

            per_layer_curve[str(display_layer)] = [None if np.isnan(v) else round(float(v), 5) for v in residue_curve]

        results["domains"][motif_id] = {
            "window_length": window_len,
            "n_layers_available": len(per_layer_curve),
            "per_layer_curve": per_layer_curve,
        }
        print(f"{motif_id}: {len(per_layer_curve)}/{NUM_LAYERS} layers scored", file=sys.stderr)

    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_JSON, "w") as f:
        json.dump(results, f)
    print(f"Wrote {OUT_JSON}", file=sys.stderr)


if __name__ == "__main__":
    main()
