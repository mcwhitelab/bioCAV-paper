#! /bin/bash

# -------------
### Directives
# -------------
#SBATCH --job-name=peps12_cofold
#SBATCH --account=mcwhite
#SBATCH --partition=gpu_standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G
#SBATCH --time=24:00:00
#SBATCH --gres=gpu:1
#SBATCH --array=0-5
#SBATCH --mail-type=END,FAIL

# -------------
### Purpose
# -------------
# Does peps12 occupy the RpoZ/omega binding site? Sequence CAVs could not
# answer this (they encode omega-FAMILY identity, and peps12 is hypothesised
# to be a non-homologous convergent mimic — outside their validated reach).
# Co-folding brings in the information the CAVs lacked: the partner.
#
# The omega site on beta-prime is discontinuous (E. coli contacts at 362-364,
# 414-419, 438, 473-488, 614-618, 905-913, 1359-1364, measured at 5A from
# PDB 6ALF chains J/K). Residues 340-640 capture 30 of 33 contacts (91%),
# so that trimmed domain is the docking scaffold — the full plastid
# equivalent (rpoC1A+rpoC1B+rpoC2 ~4200 aa) is far too large to co-fold.
#
# B1 is the load-bearing control: if AF2 cannot reproduce the KNOWN
# E. coli beta'/omega interaction on this trimmed domain, the scaffold has
# lost the pocket and every other run here is uninterpretable. Check B1
# before reading anything into B3/P1.
#
#   B1  posctrl  bp(340-640) + E.coli omega      <- MUST succeed
#   B2  posctrl  bp(340-640) + Synechococcus omega (cross-species)
#   B3  TEST     bp(340-640) + peps12 mature (transit peptide removed)
#   B4  negctrl  bp(340-640) + size-matched unrelated chlamy protein
#   P1  TEST     plastid rpoC1B + rpoC2(1-250) + peps12 mature
#   P2  ctrl     plastid rpoC1B + rpoC2(1-250) + Synechococcus omega
#
# Read ipTM (interface confidence) from ranking_debug.json. Rough guide:
# ipTM > 0.8 confident, 0.6-0.8 plausible, < 0.5 likely no interface.
# What matters is B3/P1 RELATIVE to B1/B2 (true positives) and B4 (decoy).

RUNDIR=/groups/clairemcwhite/cav_workspace/bioCAV-paper/interface_search/fastas/pairs/rpoz_chlamy/cofold
NAMES=(B1_posctrl_ecoliZ B2_posctrl_synZ B3_TEST_peps12 B4_negctrl_decoy P1_TEST_peps12 P2_ctrl_synZ)
NAME=${NAMES[$SLURM_ARRAY_TASK_ID]}

DB=/contrib/datasets/alphafold/2.3.0
mkdir -p "$RUNDIR/out"
cd "$RUNDIR" || exit 1

echo "=== $NAME ==="; date

apptainer run --nv \
  -B "$RUNDIR":/work \
  -B $DB:/data \
  --pwd /app/alphafold \
  /contrib/singularity/alphafold/alphafold-2.3.0.sif \
  --data_dir=/data \
  --fasta_paths=/work/inputs/${NAME}.fasta \
  --output_dir=/work/out \
  --model_preset=multimer \
  --db_preset=full_dbs \
  --max_template_date=2024-01-01 \
  --uniref90_database_path=/data/uniref90/uniref90.fasta \
  --mgnify_database_path=/data/mgnify/mgy_clusters_2022_05.fa \
  --template_mmcif_dir=/data/pdb_mmcif/mmcif_files \
  --obsolete_pdbs_path=/data/pdb_mmcif/obsolete.dat \
  --bfd_database_path=/data/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt \
  --uniref30_database_path=/data/uniref30/UniRef30_2021_03 \
  --pdb_seqres_database_path=/data/pdb_seqres/pdb_seqres.txt \
  --uniprot_database_path=/data/uniprot/uniprot.fasta \
  --norun_relax \
  --use_gpu_relax=false

echo "=== done $NAME ==="; date
