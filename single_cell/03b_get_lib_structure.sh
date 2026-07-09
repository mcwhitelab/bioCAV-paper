LIB=cav_library/b617ee1b-f8c8-4de9-b82b-e803ab93550d
#LIB=cav_library/9d8e5dca-03a3-457d-b7fb-844c75735c83
source "$(dirname "$0")/../config/paths.sh"
GEMINI_API_KEY=$GOOGLE_API_KEY



# Step 1: infer library structure (run once, needs API key)
python $BIOCAV_REPO/specific_scripts/analyze_cav_library.py \
    --lib-dir $LIB \
    --api-key $GEMINI_API_KEY
