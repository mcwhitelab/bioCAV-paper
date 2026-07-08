# retrieve_GO.sh lives in archive/ — this call was broken (looked for it in
# the current dir). Pointed at its actual location; move it out of archive/
# if it's meant to be current.
bash "$(dirname "$0")/archive/retrieve_GO.sh" GO_0005887 # Integral component of plasma membrane
#bash "$(dirname "$0")/archive/retrieve_GO.sh" GO_0005886 # cc plasma membrane
#bash "$(dirname "$0")/archive/retrieve_GO.sh" GO_0072657 # protein localization to membrane

