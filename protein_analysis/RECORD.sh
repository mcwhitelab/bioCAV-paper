source "$(dirname "$0")/../config/paths.sh"


cd $PROFAB_GO_SETS

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

cd $PROFAB_GO_SETS
wget https://profab.kansil.org/profab_data/ec_dataset_part1.tar.gz
wget https://profab.kansil.org/profab_data/ec_dataset_part2.tar.gz
wget https://profab.kansil.org/profab_data/ec_dataset_part3.tar.gz
wget https://profab.kansil.org/profab_data/ec_dataset_part4.tar.gz


# example processing
tar -xzvf go_dataset_part1.tar.gz
cd go_dataset_part1
for f in *zip; do unzip $f; rm $f; done
cd ../

From GO_sets/
bash remove_profab_files.sh 


# Make span files. select max 1000 proteins 
bash format_proFAB.sh


# Create commands for array job
bash COMMANDS_files/make_cav_commands.sh 


# Get validation proteins with temporal split
Downloaded from QuickGO
$BIOCAV_PAPER_ROOT/protein_analysis/results/temporal_eval_mf/eval_embeddings/val_proteins.span.fasta
$BIOCAV_PAPER_ROOT/protein_analysis/results/temporal_eval_bp/eval_embeddings/val_proteins.span.fasta
$BIOCAV_PAPER_ROOT/protein_analysis/results/temporal_eval_cc/eval_embeddings/val_proteins.span.fasta



# evaluate cavs by GO temporal split
#bash scripts/eval_cavs.sbatch

for ont in bp cc mf;
do

python $BIOCAV_REPO/specific_scripts/summarize_temporal_eval.py \
    --results  results/temporal_eval_${ont}/eval_temporal_results.tsv \
    --go-obo static/go_annotations/go.obo \
    --out-dir  results/temporal_eval_${ont}/ \
    --figure-data-dir figure_data \
    --label $ont

done





# Format API jsonl
dgw_dir=$AHMAD_GO_RESULTS_DIR
for go in bp cc mf;
do
python scripts/format_deepgoweb_api.py $dgw_dir/val_proteins_${go}.deepgo.jsonl outside_tools/val_proteins_${go}.deepgo.csv
done



# Format DeepGOSE

# Run on the model
# Val proteins
sbatch /xdisk/clairemcwhite/test_deepgose_val.sbatch

# Negative test proteins
sbatch /xdisk/clairemcwhite/test_deepgose.sbatch

for f in $PROFAB_GO_SETS/go_dataset_part*/*/*gz; do gunzip $f; done

# These are all the deepgose preds on negative test protes
andom_negative_test_max1000.span_preds_mf.tsv
cat $PROFAB_GO_SETS/go_dataset_part*/*/random_negative*.span_preds_bp.tsv | sort -u > outside_tools/neg_test_deepgose_preds_bp.tsv
cat $PROFAB_GO_SETS/go_dataset_part*/*/random_negative*.span_preds_mf.tsv | sort -u > outside_tools/neg_test_deepgose_preds_mf.tsv
cat $PROFAB_GO_SETS/go_dataset_part*/*/random_negative*.span_preds_cc.tsv | sort -u > outside_tools/neg_test_deepgose_preds_cc.tsv


# Combine positives and negatives
cat outside_tools/neg_test_deepgose_preds_bp.tsv outside_tools/val_proteins.span_preds_bp.tsv > outside_tools/deepgose_preds_bp.tsv 
cat outside_tools/neg_test_deepgose_preds_mf.tsv outside_tools/val_proteins.span_preds_mf.tsv > outside_tools/deepgose_preds_mf.tsv 
cat outside_tools/neg_test_deepgose_preds_cc.tsv outside_tools/val_proteins.span_preds_cc.tsv > outside_tools/deepgose_preds_cc.tsv 


for go in mf bp cc;
do
python $BIOCAV_REPO/specific_scripts/compare_tool_temporal.py \
    --tool-predictions  outside_tools/deepgose_preds_${go}.tsv \
    --results           results/temporal_eval_${go}/eval_temporal_results.tsv \
    --per-term-summary  results/temporal_eval_${go}/eval_temporal_per_term_summary.tsv \
    --out-dir           results/temporal_eval_${go}/ \
    --output            results/temporal_eval_${go}/tool_comparison_deepgose.tsv \
    --tool-format long_tsv \
    --val-pkl results/temporal_eval_${go}/eval_embeddings/val_proteins.span.fasta.pkl \
    --go-base-dirs $PROFAB_GO_SETS/go_dataset_part* \
    --scaler-pkl reference_population/scaler_v1.pkl \
    --figure-data-dir figure_data \
    --label $go \
    --llr-threshold 2.3

done



python $BIOCAV_REPO/specific_scripts/inspect_go_predictions.py \
    --go-term          GO:0070403 \
    --results          results/temporal_eval_mf/eval_temporal_results.tsv \
    --out-dir          results/temporal_eval_mf/ \
    --tool-predictions outside_tools/deepgose_preds_mf.tsv \
    --tool-format      long_tsv


for go in mf bp cc;
do
python specific_scripts/inspect_go_predictions.py \
    --all-terms \
    --results results/temporal_eval_${go}/eval_temporal_results.tsv \
    --out-dir results/temporal_eval_${go}/ \
    --tool-predictions outside_tools/deepgose_preds_${go}.tsv \
    --output results/temporal_eval_${go}/all_term_predictions.tsv
done



# For deepgoweb (unused)
for go in bp cc mf;
do
python $BIOCAV_REPO/specific_scripts/compare_tool_temporal.py \
    --tool-predictions  outside_tools/val_proteins_${go}.deepgo.csv \
    --results           results/temporal_eval_${go}/eval_temporal_results.tsv \
    --per-term-summary  results/temporal_eval_${go}/eval_temporal_per_term_summary.tsv \
    --out-dir           results/temporal_eval_${go}/ \
    --output            results/temporal_eval_${go}/tool_comparison.tsv
done





##### Format EC from EC-Bench
# Not using Price, too many gold standard labels seem wrong
python $BIOCAV_REPO/specific_scripts/reformat_ec_goldstandard.py \
    --input  static/ec_annotations/EC-Bench/test_ec.csv \
    --output static/ec_annotations/EC-Bench/test_ec_long.tsv

python $BIOCAV_REPO/specific_scripts/reformat_ec_goldstandard.py \
    --input  static/ec_annotations/EC-Bench/price-149.csv \
    --output static/ec_annotations/EC-Bench/price-149_long.tsv



# Evaluate cav for each EC gold standard
python $BIOCAV_REPO/specific_scripts/eval_ec.py \
    --gold-standard  static/ec_annotations/EC-Bench/test_ec_long.tsv\
    --ec-base-dirs   $PROFAB_EC_SETS/ec_dataset_part* \
    --out-dir        results/ec_eval/ecbenchtest/ \
    --scaler-pkl     reference_population/scaler_v1.pkl \
    --min-n          1


python $BIOCAV_REPO/specific_scripts/eval_ec.py \
    --gold-standard  static/ec_annotations/EC-Bench/price-149_long.tsv\
    --ec-base-dirs   $PROFAB_EC_SETS/ec_dataset_part* \
    --out-dir        results/ec_eval/ecbenchprice/ \
    --scaler-pkl     reference_population/scaler_v1.pkl \
    --val-fasta      static/ec_annotations/EC-Bench/price-149.fasta \
    --min-n          1


python $BIOCAV_REPO/specific_scripts/summarize_ec_eval.py \
    --results  results/ec_eval/ecbenchtest/eval_ec_results.tsv \
    --out-dir  results/ec_eval/ecbenchtest/ \
        --exclude           static/ec_annotations/EC-Bench/ec_train_overlap_pairs.tsv

python $BIOCAV_REPO/specific_scripts/summarize_ec_eval.py \
    --results  results/ec_eval/ecbenchprice/eval_ec_results.tsv \
    --out-dir  results/ec_eval/ecbenchprice/


# Compare to external tools 




python $BIOCAV_REPO/specific_scripts/compare_ec_tools.py \
    --tool-predictions  static/ec_annotations/EC-Bench/results_test_100_go_3.csv \
    --out-dir           results/ec_eval/ecbenchtest \
    --cav-results       results/ec_eval/ecbenchtest/eval_ec_results.tsv \
    --exclude           static/ec_annotations/EC-Bench/ec_train_overlap_pairs.tsv \
    --llr-threshold     2.3 \
    --val-pkl           results/ec_eval/ecbenchtest/eval_embeddings/val_proteins.span.fasta.pkl \
    --ec-base-dirs      $PROFAB_EC_SETS/ec_dataset_part* \
    --scaler-pkl        reference_population/scaler_v1.pkl \
    --figure-data-dir figure_data 









# Without CAV results yet

## With CAV results
#python $BIOCAV_REPO/specific_scripts/compare_ec_tools.py \
#    --tool-predictions   static/ec_annotations/EC-Bench/results_test_30_go_3.csv  \
#    --out-dir           results/ec_eval/ecbenchtest \
#    --cav-results       results/ec_eval/ecbenchtest/eval_ec_results.tsv \
#        --exclude           static/ec_annotations/EC-Bench/ec_train_overlap_pairs.tsv
#
#
## This one is most fair 
#python $BIOCAV_REPO/specific_scripts/compare_ec_tools.py \
#    --tool-predictions   static/ec_annotations/EC-Bench/results_test_100_go_3.csv  \
#    --out-dir           results/ec_eval/ecbenchtest \
#    --cav-results       results/ec_eval/ecbenchtest/eval_ec_results.tsv \
#    --exclude           static/ec_annotations/EC-Bench/ec_train_overlap_pairs.tsv \
#    --llr-threshold 1
#
#
## Don't trust this gold standard, the true labels don't match the structural matches...
#python $BIOCAV_REPO/specific_scripts/compare_ec_tools.py \
#    --tool-predictions   static/ec_annotations/EC-Bench/results_price_149_30_go_3.csv  \
#    --out-dir           results/ec_eval/ecbenchprice \
#    --cav-results       results/ec_eval/ecbenchprice/eval_ec_results.tsv
#
#python $BIOCAV_REPO/specific_scripts/compare_ec_tools.py \
#    --tool-predictions   static/ec_annotations/EC-Bench/results_price_149_100_go_3.csv  \
#    --out-dir           results/ec_eval/ecbenchprice \
#    --cav-results       results/ec_eval/ecbenchprice/eval_ec_results.tsv
#
#
#
#python $BIOCAV_REPO/specific_scripts/specificity_ec_rank.py \
#    --val-pkl       results/ec_eval/ecbenchtest/eval_embeddings/val_proteins.span.fasta.pkl \
#    --gold-standard static/ec_annotations/EC-Bench/test_ec_long.tsv \
#    --ec-base-dirs  $PROFAB_EC_SETS/ec_dataset_part* \
#    --scaler-pkl    reference_population/scaler_v1.pkl \
#    --out-dir       results/ec_eval/ecbenchtest/ \
#    --exclude       static/ec_annotations/EC-Bench/ec_train_overlap_pairs.tsv



Rscript figures.R



