source "$(dirname "$0")/../../config/paths.sh"

python - <<'EOF'
import os, sys
biocav_repo = os.environ["BIOCAV_REPO"]
sys.path.insert(0, os.path.join(biocav_repo, "specific_scripts"))
sys.path.insert(0, os.path.join(biocav_repo, "core"))

from build_reference_population import fit_global_pca
from src.utils.data_loader import load_sequence_embeddings

embs, ids = load_sequence_embeddings("reference_population/neg_10000.pkl")
fit_global_pca(embs, "reference_population/global_pca_v1.pkl", pca_dim=128)
EOF


