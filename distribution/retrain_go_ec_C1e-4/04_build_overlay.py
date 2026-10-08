#!/usr/bin/env python
"""
Symlink overlay of proFAB with the C=1e-4 CAVs (2026-10-01), so the unchanged Fig 2/3
pipelines (which take --go-base-dirs / --ec-base-dirs) score with the new library.

  /xdisk/clairemcwhite/clairemcwhite/go_ec_cavs_C1e-4_overlay/{GO_sets,EC_sets}/<part>/<term>/
    every proFAB file of the term dir    -> symlink to proFAB (spans, fastas, test pkls)
    random_positive_train_max1000_cav/
      concept_v1.npy, report_v1.json     -> symlink to the C=1e-4 library
      scaler_v1.pkl, manifest_v1.json    -> symlink to the shipped CAV dir (loaded by
                                            load_cav_artifacts but never used to score:
                                            every pipeline passes the shared scaler_v1.pkl)
Every term with a C=1e-4 CAV gets a CAV dir: the shipped set (6,677 GO + 952 EC) plus the
two new-only CAVs, EC 3.4.22.28 (positives truncated to 2046 aa) and GO:0032436 (added
2026-10-01 at the user's request). Those two have no shipped dir, so their scaler_v1.pkl
links to the shared reference scaler (again loaded but never used to score).
Idempotent; never writes into proFAB.
"""
import os
from pathlib import Path

PROFAB = Path("/xdisk/clairemcwhite/proFAB")
NEW = Path("/xdisk/clairemcwhite/clairemcwhite/go_ec_cavs_C1e-4_sharedscaler")
OV = Path("/xdisk/clairemcwhite/clairemcwhite/go_ec_cavs_C1e-4_overlay")
SUB = "random_positive_train_max1000_cav"
SHARED_SCALER = Path("/groups/clairemcwhite/cav_workspace/bioCAV-paper/shared_protein/reference_population/scaler_v1.pkl")


def link(src, dst):
    if not dst.is_symlink():
        os.symlink(src, dst)


n = {"terms": 0, "cavs": 0}
for sets, pat in [("GO_sets", "go_dataset_part*/GO_*"), ("EC_sets", "ec_dataset_part*/ecNo_*")]:
    for d in sorted((PROFAB / sets).glob(pat)):
        if not d.is_dir():
            continue
        rel = d.relative_to(PROFAB)
        o = OV / rel
        o.mkdir(parents=True, exist_ok=True)
        n["terms"] += 1
        for f in d.iterdir():
            if f.name != SUB:
                link(f, o / f.name)
        new = NEW / rel / SUB
        shipped = (d / SUB / "concept_v1.npy").exists()
        if shipped:
            assert (new / "concept_v1.npy").exists(), new
        if (new / "concept_v1.npy").exists():
            oc = o / SUB
            oc.mkdir(exist_ok=True)
            for f in ["concept_v1.npy", "report_v1.json"]:
                link(new / f, oc / f)
            if shipped:
                for f in ["scaler_v1.pkl", "manifest_v1.json"]:
                    if (d / SUB / f).exists():
                        link(d / SUB / f, oc / f)
            else:
                link(SHARED_SCALER, oc / "scaler_v1.pkl")
                n["new_only"] = n.get("new_only", []) + [str(rel)]
            n["cavs"] += 1
print(n)
