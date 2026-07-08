source "$(dirname "$0")/../../config/paths.sh"

rm cav_EC_COMMANDS.sh
touch cav_EC_COMMANDS.sh
for f in $PROFAB_EC_SETS/ec*/e*/random_positive_train_max1000.span;

do
      echo "bash $BIOCAV_REPO/specific_scripts/train_cav_from_span.sh $f" >> cav_EC_COMMANDS.sh

done

rm cav_GO_COMMANDS.sh
touch cav_GO_COMMANDS.sh
for f in $PROFAB_GO_SETS/go*/G*/random_positive_train_max1000.span;

do
      echo "bash $BIOCAV_REPO/specific_scripts/train_cav_from_span.sh $f" >> cav_GO_COMMANDS.sh

done 



