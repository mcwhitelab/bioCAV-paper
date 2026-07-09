

source "$(dirname "$0")/../config/paths.sh"
conda activate /groups/clairemcwhite/envs/$CONDA_ENV/


for dir in $PROFAB_EC_SETS/ec_dataset_part*/ec*;
do

        echo $dir
	python $BIOCAV_REPO/specific_scripts/format_proFAB.py --data-dir $dir --max-n 1000 --seed 42

done

for dir in $PROFAB_GO_SETS/go_dataset_part*/G*;
do


	python $BIOCAV_REPO/specific_scripts/format_proFAB.py --data-dir $dir --max-n 1000 --seed 42

done
