

source "$(dirname "$0")/../../config/paths.sh"
conda activate /groups/clairemcwhite/envs/$CONDA_ENV/


for dir in static/EC_sets/ec_dataset_part*/ec*;
do

        echo $dir
	python $BIOCAV_REPO/specific_scripts/format_proFAB.py --data-dir $dir --max-n 1000 --seed 42

done

for dir in static/GO_sets/go_dataset_part*/G*;
do


	python $BIOCAV_REPO/specific_scripts/format_proFAB.py --data-dir $dir --max-n 1000 --seed 42

done
