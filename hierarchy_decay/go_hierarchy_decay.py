#!/usr/bin/env python3
"""
go_hierarchy_decay.py — Supplemental figure data: for a handful of
validation-positive proteins with a specific (leaf) true GO molecular-
function term, show how the CAV projection score for that protein changes
as you score it against progressively broader is_a ancestor CAVs (leaf ->
parent -> grandparent -> great-grandparent).

Mirrors ec_hierarchy_decay.py's design: all CAVs share one global scaler,
so raw projection scores are directly comparable across the hierarchy for
the same protein. Chains are found with find_go_chains.py.

Usage
-----
python go_hierarchy_decay.py \\
    --val-pkl  go_eval/results/temporal_eval_mf/eval_embeddings/val_proteins.span.fasta.pkl \\
    --go-obo   go_eval/static/go_annotations/go.obo \\
    --go-base-dirs "$PROFAB_GO_SETS/go_dataset_part*" \\
    --scaler-pkl shared_protein/reference_population/scaler_v1.pkl \\
    --proteins O69729:GO:0000155 A3KPE8:GO:0000175 G8BDQ3:GO:0000400 \\
               F8VPJ6:GO:0000979 O94913:GO:0000993 \\
    --out-dir hierarchy_decay/results
"""

import argparse
import glob
import sys
from pathlib import Path

import joblib
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from find_go_chains import parse_obo, find_chain


def load_sequence_embeddings(pkl_path):
    sys.path.insert(0, str(Path(__file__).parent.parent.parent / "bioCAV" / "core"))
    from src.utils.data_loader import load_sequence_embeddings as _load
    return _load(pkl_path)


def compute_projections(embeddings, cav, scaler):
    X = scaler.transform(embeddings)
    return X @ cav


def load_go_cav(go_base_dir_patterns, go_id: str):
    dir_name = go_id.replace("GO:", "GO_", 1)
    for pattern in go_base_dir_patterns:
        for base in glob.glob(pattern):
            cav_file = Path(base) / dir_name / "random_positive_train_max1000_cav" / "concept_v1.npy"
            if cav_file.exists():
                return np.load(cav_file).astype(np.float64)
    return None


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--val-pkl", required=True)
    parser.add_argument("--go-obo", required=True)
    parser.add_argument("--go-base-dirs", nargs="+", required=True)
    parser.add_argument("--scaler-pkl", required=True)
    parser.add_argument("--proteins", nargs="+", required=True,
                        help="protein_id:go_id pairs, e.g. O69729:GO:0000155")
    parser.add_argument("--out-dir", required=True)
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    parents, names = parse_obo(args.go_obo)

    has_cav = set()
    for pattern in args.go_base_dirs:
        for base in glob.glob(pattern):
            for d in Path(base).glob("GO_*"):
                if (d / "random_positive_train_max1000_cav" / "concept_v1.npy").exists():
                    has_cav.add(d.name.replace("GO_", "GO:", 1))

    val_embs, val_ids = load_sequence_embeddings(args.val_pkl)
    _SKIP = {"sp", "tr", "sw", "ref"}
    id_to_idx = {}
    for i, sid in enumerate(val_ids):
        id_to_idx[sid] = i
        for part in sid.split("|"):
            if part and part not in _SKIP:
                id_to_idx.setdefault(part, i)

    scaler = joblib.load(args.scaler_pkl)

    rows = []
    for spec in args.proteins:
        protein_id, go_id = spec.split(":", 1)
        idx = id_to_idx.get(protein_id)
        if idx is None:
            print(f"WARNING: {protein_id} not in val embeddings -- skipping")
            continue
        emb = val_embs[idx:idx + 1]

        chain = find_chain(go_id, parents, has_cav, max_depth=3)
        if chain is None:
            print(f"WARNING: {protein_id} ({go_id}) has no full 4-rung chain -- skipping")
            continue
        # chain is leaf-to-root; level 4 = leaf (most specific), level 1 = most general
        chain_levels = list(enumerate(reversed(chain), start=1))  # [(1, root), ..., (4, leaf)]

        for level, term in chain_levels:
            cav = load_go_cav(args.go_base_dirs, term)
            score = float(compute_projections(emb, cav, scaler)[0])
            rows.append({
                "protein_id": protein_id, "leaf_term": go_id,
                "leaf_name": names.get(go_id, ""),
                "level": level, "level_term": term,
                "level_name": names.get(term, ""), "score": score,
            })

    df = pd.DataFrame(rows)
    out_csv = out_dir / "go_hierarchy_decay.csv"
    df.to_csv(out_csv, index=False)
    print(f"Saved {out_csv}")
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
