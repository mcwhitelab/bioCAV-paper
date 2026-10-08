#!/usr/bin/env python3
"""
find_go_chains.py — Find (protein, GO term) pairs from a positive-pair
validation table where the GO term has an is_a ancestor chain of length 4
(the term itself + 3 ancestors) and every rung of the chain has a trained
CAV.

Usage
-----
python find_go_chains.py \\
    --go-obo go_eval/static/go_annotations/go.obo \\
    --pairs  figures/figure_data/go_specificity_ranks_mf.csv \\
    --go-base-dirs "$PROFAB_GO_SETS/go_dataset_part*"
"""

import argparse
import glob
from pathlib import Path

import pandas as pd


def parse_obo(go_obo):
    parents = {}
    names = {}
    cur_id = None
    with open(go_obo) as f:
        for line in f:
            line = line.rstrip("\n")
            if line == "[Term]":
                cur_id = None
            elif line.startswith("id: "):
                cur_id = line[4:].strip()
                parents.setdefault(cur_id, [])
            elif line.startswith("name: ") and cur_id:
                names[cur_id] = line[6:].strip()
            elif line.startswith("is_a: ") and cur_id:
                parent = line[6:].split("!")[0].strip()
                parents[cur_id].append(parent)
    return parents, names


def find_chain(term, parents, has_cav, max_depth=3):
    """DFS for an ancestor chain [term, parent, grandparent, ...] of length
    max_depth+1 where every node has a trained CAV. Returns the chain
    (leaf-to-root order) or None."""
    def dfs(node, depth):
        if depth == max_depth:
            return [node]
        for p in parents.get(node, []):
            if p in has_cav:
                sub = dfs(p, depth + 1)
                if sub is not None:
                    return [node] + sub
        return None
    if term not in has_cav:
        return None
    return dfs(term, 0)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--go-obo", required=True)
    parser.add_argument("--pairs", required=True)
    parser.add_argument("--go-base-dirs", nargs="+", required=True)
    parser.add_argument("--n-show", type=int, default=30)
    args = parser.parse_args()

    parents, names = parse_obo(args.go_obo)

    has_cav = set()
    for pattern in args.go_base_dirs:
        for base in glob.glob(pattern):
            for d in Path(base).glob("GO_*"):
                if (d / "random_positive_train_max1000_cav" / "concept_v1.npy").exists():
                    has_cav.add(d.name.replace("GO_", "GO:", 1))
    print(f"{len(has_cav)} GO terms have a trained CAV")

    pairs = pd.read_csv(args.pairs)
    found = []
    seen_terms = {}
    for _, row in pairs.iterrows():
        term = row["go_term"]
        if term not in seen_terms:
            seen_terms[term] = find_chain(term, parents, has_cav)
        chain = seen_terms[term]
        if chain is not None:
            found.append({
                "protein_id": row["protein_id"], "go_term": term,
                "go_name": names.get(term, ""), "chain": " -> ".join(chain),
            })

    out = pd.DataFrame(found)
    print(f"\n{len(out)} (protein, term) pairs with a full 4-rung CAV'd ancestor chain")
    print(out.head(args.n_show).to_string(index=False))


if __name__ == "__main__":
    main()
