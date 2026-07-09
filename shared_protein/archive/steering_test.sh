# 300 aa protein, append 23 new TM residues
python ../github/TCAV_ss/specific_scripts/cav_guided_generation.py \
    --sequence gfp.fasta \
    --model /groups/clairemcwhite/models/ESMplusplus_large \
    --cav-dir ./cavs/GO_0005887/ \
    --mask-start 200 --mask-end 250 \
    --cav-weight 0.0 \
    --out generated_cterm_tm.fasta

# Baseline (no steering)
#python specific_scripts/cav_guided_generation.py ... --cav-weight 0.0 --out baseline.fasta

