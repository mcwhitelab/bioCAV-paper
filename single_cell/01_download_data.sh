#ds=3f7c572c-cd73-4b51-a313-207c7f20f188
ds=b617ee1b-f8c8-4de9-b82b-e803ab93550d
apikey=$GOOGLE_API_KEY

source "$(dirname "$0")/../config/paths.sh"

python $BIOCAV_REPO/specific_scripts/build_cav_library_agent.py \
    --dataset-id $ds \
    --prompt "Inspect the cell_type, tissue, and disease metadata columns. For every combination of (cell_type × tissue × disease_state) where disease_state is normal/healthy or tumor/cancer (use whatever values exist in the disease column): create one CAV named '<cell_type>__<tissue>__<disease_state>' using cells matching that combination. Replace spaces with underscores in the name. Use max 500 cells per group, skip if fewer than 50. Do not include blood, PBMC, or lymph node samples. Do not pool across tissues or disease states." \
    --api-key $apikey  


#ds=9d8e5dca-03a3-457d-b7fb-844c75735c83
#python $BIOCAV_REPO/specific_scripts/build_cav_library_agent.py \
#    --dataset-id $ds \
#    --prompt "Create one-vs-rest CAVs for all cell types in this dataset. \
#              Skip any cell type with fewer than 50 cells. \
#              Use max 500 cells per group." \
#    --api-key $api_key 
    

