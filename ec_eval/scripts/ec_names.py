#!/usr/bin/env python3
"""
ec_names.py — Parse ExPASy's ENZYME nomenclature files into an
EC-number -> plain-English-name lookup, covering both fully-specific
numbers (enzyme.dat) and broader class/subclass/sub-subclass categories
(enzclass.txt).

Source: https://ftp.expasy.org/databases/enzyme/{enzyme.dat,enzclass.txt}
"""

import re
from pathlib import Path


def load_ec_names(enzyme_dat: str, enzclass_txt: str) -> dict:
    names = {}

    # enzclass.txt: broad class/subclass/sub-subclass names, e.g.
    # "1. 1. 1.-    With NAD(+) or NADP(+) as acceptor."
    pat = re.compile(r"^(\S+)\.\s*(\S+)\.\s*(\S+)\.-\s+(.*)$")
    with open(enzclass_txt, encoding="latin-1") as f:
        for line in f:
            m = pat.match(line.rstrip("\n"))
            if not m:
                continue
            a, b, c, desc = m.groups()
            parts = [p for p in (a, b, c) if p != "-"]
            key = ".".join(parts)
            names[key] = desc.strip().rstrip(".")

    # enzyme.dat: fully-specific EC numbers, e.g.
    # ID   1.1.1.1
    # DE   alcohol dehydrogenase.
    cur_id = None
    with open(enzyme_dat, encoding="latin-1") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("ID   "):
                cur_id = line[5:].strip()
            elif line.startswith("DE   ") and cur_id:
                desc = line[5:].strip().rstrip(".")
                names[cur_id] = desc
                cur_id = None

    return names


if __name__ == "__main__":
    import sys
    names = load_ec_names(sys.argv[1], sys.argv[2])
    print(f"Loaded {len(names)} EC names")
    for k in ["1", "1.1", "1.1.1", "1.1.1.1", "3.4.14.11"]:
        print(f"  {k}: {names.get(k, '<missing>')}")
