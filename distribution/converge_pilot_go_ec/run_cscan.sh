#!/bin/bash
# Wait for the CPU embed array (24048937), then fit (02b), evaluate (03b) and summarise (04).
cd "$(dirname "$0")"
PY=/groups/clairemcwhite/envs/core_pkgs4/bin/python
while squeue -h -j 24048937 2>/dev/null | grep -q .; do sleep 120; done
echo "pooled positives: $(ls /xdisk/clairemcwhite/clairemcwhite/converge_pilot_go_ec/pos/*.pos.npy 2>/dev/null | wc -l)/45"
OMP_NUM_THREADS=8 $PY 02b_fit_cscan.py > fit_cscan.log 2>&1; echo "fit exit $?"; tail -1 fit_cscan.log
OMP_NUM_THREADS=8 $PY 03b_evaluate_cscan.py > eval_cscan.log 2>&1; echo "eval exit $?"; tail -2 eval_cscan.log
$PY 04_summarise_cscan.py
