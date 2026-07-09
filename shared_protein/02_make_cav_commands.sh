source "$(dirname "$0")/../config/paths.sh"

EC_COMMANDS="$BIOCAV_PAPER_ROOT/ec_eval/COMMANDS_files/cav_EC_COMMANDS.sh"
GO_COMMANDS="$BIOCAV_PAPER_ROOT/go_eval/COMMANDS_files/cav_GO_COMMANDS.sh"

rm -f "$EC_COMMANDS"
touch "$EC_COMMANDS"
for f in $PROFAB_EC_SETS/ec*/e*/random_positive_train_max1000.span;

do
      echo "bash $BIOCAV_REPO/specific_scripts/train_cav_from_span.sh $f" >> "$EC_COMMANDS"

done

rm -f "$GO_COMMANDS"
touch "$GO_COMMANDS"
for f in $PROFAB_GO_SETS/go*/G*/random_positive_train_max1000.span;

do
      echo "bash $BIOCAV_REPO/specific_scripts/train_cav_from_span.sh $f" >> "$GO_COMMANDS"

done



