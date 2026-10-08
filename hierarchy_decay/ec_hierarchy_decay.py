#!/usr/bin/env python3
"""
ec_hierarchy_decay.py — Supplemental figure data: for a handful of
validation-positive proteins with a fully-specific (level-4) true EC
number, show how the CAV projection score for that protein changes as you
score it against progressively broader ancestor CAVs (level 1 -> 2 -> 3 -> 4
of the same EC lineage, e.g. "3" -> "3.1" -> "3.1.1" -> "3.1.1.3").

All CAVs share one global scaler (shared_protein/reference_population/
scaler_v1.pkl), so raw projection scores are directly comparable across
levels for the same protein without needing a level-specific reference
(positive/negative) distribution -- this is what makes it feasible to pick
proteins outside the small subset of EC terms with cached LLR reference
scores.

Usage
-----
python ec_hierarchy_decay.py \\
    --val-pkl static/ec_annotations/EC-Bench/test_ec.fasta.pkl \\
    --gold    static/ec_annotations/EC-Bench/test_ec_long.tsv \\
    --ec-base-dirs "$PROFAB_EC_SETS/ec_dataset_part*" \\
    --scaler-pkl shared_protein/reference_population/scaler_v1.pkl \\
    --proteins A0A452E9Y6:1.11.1.7 C0HLQ6:1.2.1.12 P0DPS1:2.7.7.7 \\
               P0DUB9:3.1.1.3 A0A345BJN1:4.3.1.24 \\
    --out-dir hierarchy_decay/results
"""

import argparse
import glob
import sys
from pathlib import Path

import joblib
import numpy as np
import pandas as pd


def load_sequence_embeddings(pkl_path):
    sys.path.insert(0, str(Path(__file__).parent.parent.parent / "bioCAV" / "core"))
    from src.utils.data_loader import load_sequence_embeddings as _load
    return _load(pkl_path)


def compute_projections(embeddings, cav, scaler):
    X = scaler.transform(embeddings)
    return X @ cav


def load_ec_cav_lineage(ec_base_dir_patterns, ec4: str) -> dict:
    """Return {level: cav_vector} for levels 1-4 of the given EC number's
    lineage, or an empty dict if any level is missing a trained CAV."""
    parts = ec4.split(".")
    result = {}
    for base_pattern in ec_base_dir_patterns:
        for base in glob.glob(base_pattern):
            for level in [1, 2, 3, 4]:
                if level in result:
                    continue
                key = ".".join(parts[:level])
                d = Path(base) / f"ecNo_{key.replace('.', '-')}"
                cav_file = d / "random_positive_train_max1000_cav" / "concept_v1.npy"
                if cav_file.exists():
                    result[level] = np.load(cav_file).astype(np.float64)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--val-pkl", required=True)
    parser.add_argument("--gold", required=True)
    parser.add_argument("--ec-base-dirs", nargs="+", required=True)
    parser.add_argument("--scaler-pkl", required=True)
    parser.add_argument("--proteins", nargs="+", required=True,
                        help="protein_id:ec_number pairs, e.g. P0DUB9:3.1.1.3")
    parser.add_argument("--out-dir", required=True)
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

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
        protein_id, ec4 = spec.split(":")
        idx = id_to_idx.get(protein_id)
        if idx is None:
            print(f"WARNING: {protein_id} not in val embeddings -- skipping")
            continue
        emb = val_embs[idx:idx + 1]

        cavs = load_ec_cav_lineage(args.ec_base_dirs, ec4)
        if len(cavs) != 4:
            print(f"WARNING: {protein_id} ({ec4}) missing lineage levels "
                  f"{set([1,2,3,4]) - set(cavs)} -- skipping")
            continue

        for level, cav in cavs.items():
            score = float(compute_projections(emb, cav, scaler)[0])
            level_ec = ".".join(ec4.split(".")[:level])
            rows.append({
                "protein_id": protein_id, "true_ec": ec4,
                "level": level, "level_ec": level_ec, "score": score,
            })

    df = pd.DataFrame(rows)
    out_csv = out_dir / "ec_hierarchy_decay.csv"
    df.to_csv(out_csv, index=False)
    print(f"Saved {out_csv}")
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
