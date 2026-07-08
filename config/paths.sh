#!/usr/bin/env bash
# Central path/config source for bioCAV-paper driver scripts.
# Usage: source "$(dirname "$0")/../config/paths.sh"   (adjust relative depth per caller)
#
# BIOCAV_REPO / BIOCAV_PAPER_ROOT are resolved relative to THIS file's own
# location (not the caller's cwd), so this works the same whether invoked
# interactively, from an sbatch script, or from an array-job task with an
# unpredictable working directory.

_PATHS_SH_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIOCAV_PAPER_ROOT="$(cd "${_PATHS_SH_DIR}/.." && pwd)"
BIOCAV_REPO="$(cd "${_PATHS_SH_DIR}/../../bioCAV" && pwd)"
export BIOCAV_PAPER_ROOT BIOCAV_REPO

# --- Tools external to both repos ---
ESM_MODEL_DIR=/groups/clairemcwhite/models/ESMplusplus_large
HF_EMBED_SCRIPT=/groups/clairemcwhite/claire_workspace/github/mcwlab_utils/hf_embed_new.py
CONDA_ENV=core_pkgs4
export ESM_MODEL_DIR HF_EMBED_SCRIPT CONDA_ENV

# --- Other-lab-member / shared scratch resources ---
# These are genuinely external (owned by other lab members or shared cluster
# scratch); they stay real absolute paths. Centralizing them here just means
# editing one place instead of a dozen scripts.
REFERENCE_NEG_FASTA=/groups/clairemcwhite/ahmad_workspace/esm_c/neg_data/neg_10000.fasta
AHMAD_CAV_OUTPUTS_DIR=/groups/clairemcwhite/ahmad_workspace/esm_c/tcav_outputs_esmplusplus_all
AHMAD_GO_RESULTS_DIR=/groups/clairemcwhite/ahmad_workspace/go_results
GENEFORMER_MODEL_DIR=/groups/clairemcwhite/rshaw_workspace/stiffness_prediction/geneformer_model
PROFAB_GO_SETS=/xdisk/clairemcwhite/proFAB/GO_sets
PROFAB_EC_SETS=/xdisk/clairemcwhite/proFAB/EC_sets
UNIPROT_HUMAN_FASTA=/xdisk/clairemcwhite/clairemcwhite/uniprot_human_all.fasta
NEXTPROT_MF_FASTA=/xdisk/clairemcwhite/nextprot/data/mf/nextprot_mf.fasta
MOTIF_AGENT_FASTAS=/groups/clairemcwhite/claire_workspace/github/motif_agent/fastas
export REFERENCE_NEG_FASTA AHMAD_CAV_OUTPUTS_DIR AHMAD_GO_RESULTS_DIR GENEFORMER_MODEL_DIR
export PROFAB_GO_SETS PROFAB_EC_SETS UNIPROT_HUMAN_FASTA NEXTPROT_MF_FASTA MOTIF_AGENT_FASTAS
