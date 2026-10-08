#!/usr/bin/env python3
"""
build_go_namespace_lists.py — Write MF/BP/CC include-list text files for
cav_viz_pattern.py's --include-list flag.

GO terms aren't split into MF/BP/CC subdirectories on disk -- all GO CAVs
live under a single flat glob (go_dataset_part*/GO_*/.../concept_v1.npy).
This script parses go.obo for each term's `namespace:` field and writes
one plain-text file per namespace (one concept name per line, matching
the "GO_0003824"-style names cav_viz_pattern.py derives from the path),
so a single glob can be restricted to one namespace via --include-list
without any directory restructuring or symlinks.

Usage
-----
python build_go_namespace_lists.py --out-dir go_namespace_lists
"""

import argparse
import glob
from pathlib import Path

GO_OBO_DEFAULT    = "/groups/clairemcwhite/cav_workspace/bioCAV-paper/go_eval/static/go_annotations/go.obo"
PROFAB_GO_DEFAULT = "/xdisk/clairemcwhite/proFAB/GO_sets"

NS_TO_SHORT = {
    "molecular_function": "mf",
    "biological_process": "bp",
    "cellular_component": "cc",
}


def parse_go_namespaces(go_obo: str) -> dict:
    """Return {'GO:0000001': 'biological_process', ...} parsed from go.obo."""
    ns_map = {}
    cur_id = None
    with open(go_obo) as f:
        for line in f:
            line = line.rstrip("\n")
            if line == "[Term]":
                cur_id = None
            elif line.startswith("id: "):
                cur_id = line[4:].strip()
            elif line.startswith("namespace: ") and cur_id:
                ns_map[cur_id] = line[len("namespace: "):].strip()
    return ns_map


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--profab-go", default=PROFAB_GO_DEFAULT)
    parser.add_argument("--go-obo", default=GO_OBO_DEFAULT)
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    ns_map = parse_go_namespaces(args.go_obo)
    print(f"Parsed {len(ns_map)} GO id->namespace entries from {args.go_obo}")

    paths = glob.glob(f"{args.profab_go}/go_dataset_part*/GO_*/random_positive_train_max1000_cav/concept_v1.npy")
    print(f"Found {len(paths)} GO concept files under {args.profab_go}")

    names_by_ns = {short: [] for short in NS_TO_SHORT.values()}
    skipped = 0
    for p in paths:
        go_dir_name = Path(p).parent.parent.name        # e.g. "GO_0003824"
        go_id = go_dir_name.replace("GO_", "GO:", 1)
        ns = ns_map.get(go_id)
        if ns not in NS_TO_SHORT:
            skipped += 1
            continue
        names_by_ns[NS_TO_SHORT[ns]].append(go_dir_name)

    for short, names in names_by_ns.items():
        out_file = out_dir / f"{short}_ids.txt"
        out_file.write_text("\n".join(sorted(names)) + "\n")
        print(f"  {short}: {len(names)} concepts -> {out_file}")
    print("Skipped (no namespace match, e.g. obsolete terms):", skipped)


if __name__ == "__main__":
    main()
