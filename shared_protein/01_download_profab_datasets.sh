source "$(dirname "$0")/../config/paths.sh"

cd "$PROFAB_GO_SETS"
wget https://profab.kansil.org/profab_data/go_dataset_part1.tar.gz
wget https://profab.kansil.org/profab_data/go_dataset_part2.tar.gz
wget https://profab.kansil.org/profab_data/go_dataset_part3.tar.gz
wget https://profab.kansil.org/profab_data/go_dataset_part4.tar.gz
wget https://profab.kansil.org/profab_data/go_dataset_part5.tar.gz
wget https://profab.kansil.org/profab_data/go_dataset_part6.tar.gz
wget https://profab.kansil.org/profab_data/go_dataset_part7.tar.gz
wget https://profab.kansil.org/profab_data/go_dataset_part8.tar.gz
wget https://profab.kansil.org/profab_data/go_dataset_part9.tar.gz
wget https://profab.kansil.org/profab_data/go_dataset_part10.tar.gz

cd "$PROFAB_EC_SETS"
wget https://profab.kansil.org/profab_data/ec_dataset_part1.tar.gz
wget https://profab.kansil.org/profab_data/ec_dataset_part2.tar.gz
wget https://profab.kansil.org/profab_data/ec_dataset_part3.tar.gz
wget https://profab.kansil.org/profab_data/ec_dataset_part4.tar.gz

# Example extraction (repeat per part/dir as needed)
cd "$PROFAB_GO_SETS"
tar -xzvf go_dataset_part1.tar.gz
cd go_dataset_part1
for f in *zip; do unzip "$f"; rm "$f"; done
cd ../

# remove_profab_files.sh ships inside the proFAB tarballs themselves (not part
# of this repo) - run it from GO_sets/ after extraction to strip files proFAB
# doesn't need.
cd "$PROFAB_GO_SETS"
bash remove_profab_files.sh
