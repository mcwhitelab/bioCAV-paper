#!/usr/bin/env python3
"""
go_ec_cosine_validation.py — Cross-ontology validation of CAV concept
vectors using the curated GO<->EC mapping.

GO and EC CAVs are trained completely independently -- different label
sets, different positive-example spans, no shared supervision. But a
subset of GO molecular_function terms are curated by the GO consortium as
referring to the *same* underlying enzymatic activity as a specific EC
number (the ec2go mapping, https://current.geneontology.org/ontology/external2go/ec2go).

If CAV training recovers real biological structure rather than
ontology-specific artifacts, a GO term's concept vector should be most
cosine-similar to its annotated EC partner's concept vector, out of all
trained EC concept vectors -- exactly the same kind of specificity check
used for the protein-level EC rank analysis (ec_eval), but here applied
directly to the concept vectors themselves as an independent, protein-free
validation.

For each GO term with a trained CAV that ec2go maps to an EC number that
also has a trained CAV:
  - cosine similarity to its annotated EC partner ("matched")
  - cosine similarity to every other trained EC CAV ("background")
  - rank of the annotated partner among all trained EC CAVs

Usage
-----
python go_ec_cosine_validation.py \\
    --ec2go       static/ec2go.txt \\
    --go-base-dirs "$PROFAB_GO_SETS/go_dataset_part*" \\
    --ec-base-dirs "$PROFAB_EC_SETS/ec_dataset_part*" \\
    --out-dir     results
"""

import argparse
import glob
import logging
import re
from pathlib import Path

import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

CAV_SUBDIR = "random_positive_train_max1000_cav"
CONCEPT_FILE = "concept_v1.npy"


def strip_wildcards(ec: str) -> str:
    parts = ec.strip().split(".")
    while parts and parts[-1].strip() in ("-", ""):
        parts.pop()
    return ".".join(parts)


def parse_ec2go(path: str) -> pd.DataFrame:
    """Parse ec2go.txt lines like:
    'EC:1.1.1.1 > GO:alcohol dehydrogenase (NAD+) activity ; GO:0004022'
    """
    rows = []
    pat = re.compile(r"^EC:(\S+)\s*>\s*GO:.*;\s*(GO:\d+)\s*$")
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("!"):
                continue
            m = pat.match(line)
            if m:
                ec_number, go_id = m.group(1), m.group(2)
                rows.append({"ec_number": strip_wildcards(ec_number), "go_id": go_id})
    df = pd.DataFrame(rows).drop_duplicates()
    logger.info(f"Parsed {len(df)} EC<->GO mapping rows from {path}")
    return df


def load_cav_vectors(base_dir_patterns: list, dir_prefix: str) -> dict:
    """Load and unit-normalise CAV vectors, keyed by their normalised id
    (e.g. 'GO:0004022' or '1.1.1.1')."""
    result = {}
    for pattern in base_dir_patterns:
        for base in sorted(glob.glob(pattern)):
            base = Path(base)
            for concept_dir in sorted(base.glob(f"{dir_prefix}*")):
                cav_file = concept_dir / CAV_SUBDIR / CONCEPT_FILE
                if not cav_file.exists():
                    continue
                raw_name = concept_dir.name  # e.g. "GO_0004022" or "ecNo_1-1-1-1"
                if dir_prefix == "GO_":
                    key = raw_name.replace("GO_", "GO:", 1)
                else:  # ecNo_
                    key = strip_wildcards(raw_name.replace("ecNo_", "", 1).replace("-", "."))
                v = np.load(cav_file).astype(np.float64)
                norm = np.linalg.norm(v)
                if norm > 1e-10:
                    result[key] = v / norm
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--ec2go", required=True)
    parser.add_argument("--go-base-dirs", nargs="+", required=True,
                        help="Glob pattern(s) for GO dataset part dirs.")
    parser.add_argument("--ec-base-dirs", nargs="+", required=True,
                        help="Glob pattern(s) for EC dataset part dirs.")
    parser.add_argument("--out-dir", required=True)
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    mapping = parse_ec2go(args.ec2go)

    logger.info("Loading EC CAV vectors...")
    ec_vecs = load_cav_vectors(args.ec_base_dirs, "ecNo_")
    logger.info(f"  {len(ec_vecs)} EC CAVs loaded")

    logger.info("Loading GO CAV vectors...")
    go_vecs = load_cav_vectors(args.go_base_dirs, "GO_")
    logger.info(f"  {len(go_vecs)} GO CAVs loaded")

    # Restrict mapping to pairs where BOTH sides have a trained CAV
    mapping = mapping[mapping["ec_number"].isin(ec_vecs) & mapping["go_id"].isin(go_vecs)]
    mapping = mapping.drop_duplicates(subset=["go_id"])  # one EC partner per GO term
    logger.info(f"{len(mapping)} GO<->EC pairs with CAVs on both sides")

    if mapping.empty:
        logger.error("No overlapping GO<->EC pairs with trained CAVs found.")
        return

    ec_names = sorted(ec_vecs.keys())
    ec_matrix = np.vstack([ec_vecs[n] for n in ec_names])   # (n_ec, dim)
    ec_index = {n: i for i, n in enumerate(ec_names)}

    rows = []
    for _, row in mapping.iterrows():
        go_id, ec_number = row["go_id"], row["ec_number"]
        go_vec = go_vecs[go_id]                       # unit vector, (dim,)
        sims = ec_matrix @ go_vec                      # cosine sim to every EC CAV
        true_idx  = ec_index[ec_number]
        true_sim  = float(sims[true_idx])
        rank      = int((sims > true_sim).sum()) + 1
        rows.append({
            "go_id": go_id, "ec_number": ec_number,
            "matched_sim": true_sim, "rank": rank, "n_ec_cavs": len(ec_names),
        })

    result_df = pd.DataFrame(rows)
    result_df.to_csv(out_dir / "go_ec_cosine_validation.csv", index=False)
    logger.info(f"Saved {out_dir / 'go_ec_cosine_validation.csv'}")

    # Background: all non-matching GO-EC cosine similarities, for the
    # GO terms actually tested (excludes each GO term's own true partner).
    matched_go_ids = list(mapping["go_id"])
    go_matrix = np.vstack([go_vecs[g] for g in matched_go_ids])   # (n_go, dim)
    full_sim = go_matrix @ ec_matrix.T                             # (n_go, n_ec)
    true_idx_per_go = np.array([ec_index[e] for e in mapping["ec_number"]])
    mask = np.ones_like(full_sim, dtype=bool)
    mask[np.arange(len(matched_go_ids)), true_idx_per_go] = False
    background_sims = full_sim[mask]

    np.save(out_dir / "go_ec_background_sims.npy", background_sims)

    n_total = len(result_df)
    pct_rank1  = (result_df["rank"] == 1).mean() * 100
    pct_top10  = (result_df["rank"] <= 10).mean() * 100
    med_matched = result_df["matched_sim"].median()
    med_bg      = float(np.median(background_sims))

    from scipy import stats as spstats
    u_stat, p_val = spstats.mannwhitneyu(result_df["matched_sim"], background_sims,
                                          alternative="greater")

    logger.info("\n" + "=" * 60)
    logger.info(f"GO<->EC cosine validation: {n_total} pairs "
                f"({len(ec_names)} EC CAVs, {len(go_vecs)} GO CAVs total)")
    logger.info(f"  Rank 1 (annotated EC partner is most similar):  {pct_rank1:.1f}%")
    logger.info(f"  Rank <= 10:                                     {pct_top10:.1f}%")
    logger.info(f"  Median matched-pair cosine sim:                 {med_matched:.3f}")
    logger.info(f"  Median background (non-matching) cosine sim:    {med_bg:.3f}")
    logger.info(f"  Mann-Whitney U (matched > background):          "
                f"U={u_stat:.0f}, p={p_val:.2e}")
    logger.info("=" * 60)

    summary = {
        "n_pairs": n_total, "n_ec_cavs": len(ec_names), "n_go_cavs": len(go_vecs),
        "pct_rank1": pct_rank1, "pct_rank_le10": pct_top10,
        "median_matched_sim": med_matched, "median_background_sim": med_bg,
        "mannwhitney_u": u_stat, "mannwhitney_p": p_val,
    }
    pd.DataFrame([summary]).to_csv(out_dir / "go_ec_cosine_validation_summary.csv", index=False)
    logger.info(f"Saved {out_dir / 'go_ec_cosine_validation_summary.csv'}")


if __name__ == "__main__":
    main()
