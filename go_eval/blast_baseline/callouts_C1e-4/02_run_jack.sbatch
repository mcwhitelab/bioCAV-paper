#!/bin/bash
#SBATCH --job-name=callout_jackhmmer
#SBATCH --output=jack_%j.out
#SBATCH --account=mcwhite
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=03:00:00
cd /groups/clairemcwhite/cav_workspace/bioCAV-paper/go_eval/blast_baseline/callouts_C1e-4
/groups/clairemcwhite/envs/core_pkgs4/bin/python 02_run_jack.py
