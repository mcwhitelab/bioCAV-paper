#!/usr/bin/env python3
"""
make_rank_vs_training_size_table.py — Supplemental table: for each rank
bucket in the GO<->EC cosine validation (go_ec_cosine_validation.py), the
distribution of the annotated EC partner's training-set size (n_positive).

Shows that pairs where the annotated EC partner ranks poorly are
concentrated among the most severely under-trained EC CAVs -- i.e. the
handful of genuine failures in the cross-ontology validation are explained
by training-set size, not a biological mismatch.

Usage
-----
python make_rank_vs_training_size_table.py \\
    --pairs-csv   results/go_ec_cosine_validation.csv \\
    --ec-base-dirs "$PROFAB_EC_SETS/ec_dataset_part*" \\
    --out-dir     results
"""

import argparse
import glob
import json
from pathlib import Path

import pandas as pd


def load_n_positive(ec_base_dir_patterns: list) -> dict:
    """Return {ec_number (dotted, wildcard-stripped): n_positive} from
    each EC CAV's report_v1.json."""
    n_pos = {}
    for pattern in ec_base_dir_patterns:
        for base in glob.glob(pattern):
            for d in Path(base).glob("ecNo_*"):
                report = d / "random_positive_train_max1000_cav" / "report_v1.json"
                if not report.exists():
                    continue
                try:
                    j = json.load(open(report))
                    key = d.name.replace("ecNo_", "", 1).replace("-", ".")
                    n_pos[key] = j["concept_metrics"]["n_positive"]
                except (KeyError, json.JSONDecodeError):
                    continue
    return n_pos


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--pairs-csv", required=True,
                        help="go_ec_cosine_validation.csv from go_ec_cosine_validation.py")
    parser.add_argument("--ec-base-dirs", nargs="+", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--rank-breaks", nargs="+", type=int, default=[1, 3, 10],
                        help="Upper bounds of rank buckets before the final "
                             "'>N' bucket (default: 1 3 10 -> '1','2-3','4-10','>10').")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.pairs_csv)
    n_pos = load_n_positive(args.ec_base_dirs)
    df["n_positive"] = df["ec_number"].map(n_pos)

    n_missing = df["n_positive"].isna().sum()
    if n_missing:
        print(f"Warning: {n_missing}/{len(df)} pairs had no n_positive lookup")

    breaks = [0] + args.rank_breaks + [df["rank"].max()]
    labels = []
    prev = 0
    for b in args.rank_breaks:
        lo = prev + 1
        labels.append(str(lo) if lo == b else f"{lo}-{b}")
        prev = b
    labels.append(f">{args.rank_breaks[-1]}")

    df["rank_bucket"] = pd.cut(df["rank"], bins=breaks, labels=labels, include_lowest=True)

    table = (
        df.groupby("rank_bucket", observed=True)["n_positive"]
        .agg(n_pairs="count", min="min", median="median", max="max")
        .reset_index()
        .rename(columns={"rank_bucket": "Rank bucket"})
    )
    for col in ["min", "median", "max"]:
        table[col] = table[col].round(1)

    out_csv = out_dir / "supp_table_rank_vs_training_size.csv"
    table.to_csv(out_csv, index=False)
    print(f"Saved {out_csv}")
    print(table.to_string(index=False))


if __name__ == "__main__":
    main()
