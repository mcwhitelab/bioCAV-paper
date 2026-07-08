

EC_BASE="/xdisk/clairemcwhite/proFAB/EC_sets"
RESULTS="results/ec_eval/ecbenchprice/eval_ec_results.tsv"   # adjust if needed

n_overlap=0
n_total=0

while IFS=$'\t' read -r ec_number ec_cav_id protein_id rest; do
    [[ "$ec_number" == "ec_number" ]] && continue  # skip header
    n_total=$((n_total + 1))

    train_fasta=$(find "$EC_BASE" -path "*/${ec_cav_id}/random_positive_train_max1000.span.fasta" \
                  2>/dev/null | head -1)

    if [[ -n "$train_fasta" ]] && grep -q "$protein_id" "$train_fasta"; then
        echo "OVERLAP  $ec_number  $protein_id"
        n_overlap=$((n_overlap + 1))
    fi
done < "$RESULTS"

echo ""
echo "Total (protein, EC) pairs evaluated : $n_total"
echo "Pairs also in EC-matched training   : $n_overlap"

# Now for the Test proteins

EC_BASE="/xdisk/clairemcwhite/proFAB/EC_sets"
RESULTS="results/ec_eval/ecbenchtest/eval_ec_results.tsv"   # adjust if needed

n_overlap=0
n_total=0

while IFS=$'\t' read -r ec_number ec_cav_id protein_id rest; do
    [[ "$ec_number" == "ec_number" ]] && continue  # skip header
    n_total=$((n_total + 1))

    train_fasta=$(find "$EC_BASE" -path "*/${ec_cav_id}/random_positive_train_max1000.span.fasta" \
                  2>/dev/null | head -1)

    if [[ -n "$train_fasta" ]] && grep -q "$protein_id" "$train_fasta"; then
        echo "OVERLAP  $ec_number  $protein_id"
        n_overlap=$((n_overlap + 1))
    fi
done < "$RESULTS"

echo ""
echo "Total (protein, EC) pairs evaluated : $n_total"
echo "Pairs also in EC-matched training   : $n_overlap"
