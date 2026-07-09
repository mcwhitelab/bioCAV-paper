#!/usr/bin/env python3

import argparse
import json
import pandas as pd


def extract_go_scores(record):
    """
    Convert one JSONL record into:
        protein_id -> {GO_ID: score}
    """
    protein_id = record.get("protein_id")

    if protein_id is None:
        # fallback, in case only predictions.protein_info exists
        protein_id = record.get("predictions", {}).get("protein_info")

    scores = {}

    prediction_groups = (
        record.get("predictions", {})
              .get("functions", [])
    )

    for group in prediction_groups:
        category = group.get("name", "")
        functions = group.get("functions", [])

        for item in functions:
            # Expected format:
            # ["GO:0005215", "transporter activity", 0.571]
            if len(item) < 3:
                continue

            go_id = item[0]
            score = item[2]

            scores[go_id] = score

            # Optional: if you want category-specific duplicate-safe columns:
            # scores[f"{category}|{go_id}"] = score

    return protein_id, scores


def jsonl_to_csv(input_jsonl, output_csv, fill_value=""):
    rows = []

    with open(input_jsonl, "r") as f:
        for line_number, line in enumerate(f, start=1):
            line = line.strip()
            if not line:
                continue

            try:
                record = json.loads(line)
            except json.JSONDecodeError as e:
                raise ValueError(f"Could not parse JSON on line {line_number}: {e}")

            protein_id, scores = extract_go_scores(record)

            if protein_id is None:
                raise ValueError(f"No protein_id found on line {line_number}")

            row = {"protein_id": protein_id}
            row.update(scores)
            rows.append(row)

    df = pd.DataFrame(rows)

    # Put protein_id first, then sort GO columns
    go_cols = sorted([c for c in df.columns if c != "protein_id"])
    df = df[["protein_id"] + go_cols]

    df.to_csv(output_csv, sep=",", index=False, na_rep=fill_value)


def main():
    parser = argparse.ArgumentParser(
        description="Convert protein GO prediction JSONL to wide CSV."
    )
    parser.add_argument(
        "input_jsonl",
        help="Input JSONL file"
    )
    parser.add_argument(
        "output_csv",
        help="Output CSV file"
    )
    parser.add_argument(
        "--fill-value",
        default="",
        help="Value to use for missing GO predictions. Default: empty string"
    )

    args = parser.parse_args()

    jsonl_to_csv(
        input_jsonl=args.input_jsonl,
        output_csv=args.output_csv,
        fill_value=args.fill_value,
    )


if __name__ == "__main__":
    main()
