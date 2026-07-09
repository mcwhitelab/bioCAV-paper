# Frozen/archived: not called by any live pipeline. retrieve_GO.sh is a
# sibling in this same archive/ dir, and retrieve_GO.sh itself still has
# pre-reorg broken internal paths (not fixed, since this cluster is not
# maintained).
bash "$(dirname "$0")/retrieve_GO.sh" GO_0005887 # Integral component of plasma membrane
#bash "$(dirname "$0")/retrieve_GO.sh" GO_0005886 # cc plasma membrane
#bash "$(dirname "$0")/retrieve_GO.sh" GO_0072657 # protein localization to membrane

