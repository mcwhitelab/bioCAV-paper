# interface_search cleanup, 2026-09-28

Moved here (paths relative to interface_search/, preserved under this folder) because they are not used by the live
SERK/CIK (Fig 5 D/E) or MYC/MAD1 analyses. Live = used by export_fig5E_data.py, export_fig5E_noparalog.py,
homology_filtered_eval.py, compare_no_serk_positives.py, run_spanfile_{myc,mad1}_contacts.sh,
run_myc_query_sharedscaler.sbatch, or cited in lrr_serkcik/SESSION_FINDINGS.md / BASELINE_FRAMING.md.
Raw inputs behind live results (per-gene FASTAs, BLAST dbs, genome-wide query tables, homolog JSONs) were kept.

Restore one:   mv archive/2026-09-28_cleanup/<path> <path>        (run from interface_search/)
Restore all:   see restore.sh

| path | reason |
|---|---|
| fastas/pairs/rpoz_chlamy | case study not in paper |
| fastas/pairs/lrr | case study not in paper |
| fastas/pairs/bzip | case study not in paper |
| fastas/pairs/colicin | case study not in paper |
| run_cofold_peps12.sh | run script for an archived case study |
| run_pdb_bfactor_Im9.sh | run script for an archived case study |
| run_pdb_bfactor_RLK5.sh | run script for an archived case study |
| run_pdb_bfactor_SERKmulti.sh | run script for an archived case study |
| run_query_Im9mutants.sh | run script for an archived case study |
| run_query_ORYSJ_serk.sh | run script for an archived case study |
| run_query_taxonomy399_serk.sh | run script for an archived case study |
| run_spanfile_FLS2serkinterface.sh | run script for an archived case study |
| run_spanfile_Im9E9interface.sh | run script for an archived case study |
| run_spanfile_PSKR1serkinterface.sh | run script for an archived case study |
| run_spanfile_RLK5serkinterface.sh | run script for an archived case study |
| run_spanfile_SERKinterface_multi.sh | run script for an archived case study |
| run_window_score_chlamy_interface.sh | run script for an archived case study |
| run_window_score_chlamy_mixed.sh | run script for an archived case study |
| run_window_score_chlamy.sh | run script for an archived case study |
| slurm-22828656.out | old slurm log |
| slurm-22828720.out | old slurm log |
| slurm-22847336.out | old slurm log |
| slurm-22850716.out | old slurm log |
| slurm-22850717.out | old slurm log |
| slurm-23344302.out | old slurm log |
| slurm-23344362.out | old slurm log |
| slurm-23344363.out | old slurm log |
| slurm-23344504.out | old slurm log |
| slurm-23344570.out | old slurm log |
| slurm-23344601.out | old slurm log |
| slurm-23344679.out | old slurm log |
| slurm-23344680.out | old slurm log |
| slurm-23348088.out | old slurm log |
| slurm-23416093.out | old slurm log |
| slurm-23419352.out | old slurm log |
| slurm-23428436.out | old slurm log |
| slurm-23428499_0.out | old slurm log |
| slurm-23428612_0.out | old slurm log |
| slurm-23428620_0.out | old slurm log |
| fastas/ORYSJ_LRR.fasta | search space used only by archived lrr/ORYSJ/taxonomy399 scripts |
| fastas/ORYSJ_LRR.fasta.pkl | search space used only by archived lrr/ORYSJ/taxonomy399 scripts |
| fastas/ORYSJ_LRR.fasta.pkl.description | search space used only by archived lrr/ORYSJ/taxonomy399 scripts |
| fastas/ORYSJ_LRR.fasta.pkl.seqnames | search space used only by archived lrr/ORYSJ/taxonomy399 scripts |
| fastas/uniprotkb_IPR001611_AND_model_organism_2026_07_10.fasta | search space used only by archived lrr/ORYSJ/taxonomy399 scripts |
| fastas/uniprotkb_IPR001611_AND_model_organism_2026_07_10.fasta.pkl | search space used only by archived lrr/ORYSJ/taxonomy399 scripts |
| fastas/uniprotkb_IPR001611_AND_model_organism_2026_07_10.fasta.pkl.description | search space used only by archived lrr/ORYSJ/taxonomy399 scripts |
| fastas/uniprotkb_IPR001611_AND_model_organism_2026_07_10.fasta.pkl.seqnames | search space used only by archived lrr/ORYSJ/taxonomy399 scripts |
| fastas/uniprotkb_IPR001611_AND_model_organism_db.pdb | search space used only by archived lrr/ORYSJ/taxonomy399 scripts |
| fastas/uniprotkb_IPR001611_AND_model_organism_db.phr | search space used only by archived lrr/ORYSJ/taxonomy399 scripts |
| fastas/uniprotkb_IPR001611_AND_model_organism_db.pin | search space used only by archived lrr/ORYSJ/taxonomy399 scripts |
| fastas/uniprotkb_IPR001611_AND_model_organism_db.pjs | search space used only by archived lrr/ORYSJ/taxonomy399 scripts |
| fastas/uniprotkb_IPR001611_AND_model_organism_db.pot | search space used only by archived lrr/ORYSJ/taxonomy399 scripts |
| fastas/uniprotkb_IPR001611_AND_model_organism_db.psq | search space used only by archived lrr/ORYSJ/taxonomy399 scripts |
| fastas/uniprotkb_IPR001611_AND_model_organism_db.ptf | search space used only by archived lrr/ORYSJ/taxonomy399 scripts |
| fastas/uniprotkb_IPR001611_AND_model_organism_db.pto | search space used only by archived lrr/ORYSJ/taxonomy399 scripts |
| fastas/uniprotkb_IPR001611_AND_taxonomy_id_399_2026_07_10.fasta | search space used only by archived lrr/ORYSJ/taxonomy399 scripts |
| fastas/uniprotkb_IPR001611_AND_taxonomy_id_399_2026_07_10.fasta.pkl | search space used only by archived lrr/ORYSJ/taxonomy399 scripts |
| fastas/uniprotkb_IPR001611_AND_taxonomy_id_399_2026_07_10.fasta.pkl.description | search space used only by archived lrr/ORYSJ/taxonomy399 scripts |
| fastas/uniprotkb_IPR001611_AND_taxonomy_id_399_2026_07_10.fasta.pkl.seqnames | search space used only by archived lrr/ORYSJ/taxonomy399 scripts |
| fastas/pairs/lrr_serkcik/make_full_figure.py | scaffold/panel builder superseded by figures.R + export_fig5E_data.py |
| fastas/pairs/lrr_serkcik/figure_full_scaffold.png | scaffold/panel builder superseded by figures.R + export_fig5E_data.py |
| fastas/pairs/lrr_serkcik/make_figure_EF.py | scaffold/panel builder superseded by figures.R + export_fig5E_data.py |
| fastas/pairs/lrr_serkcik/figure_panels_EF.png | scaffold/panel builder superseded by figures.R + export_fig5E_data.py |
| fastas/pairs/lrr_serkcik/make_panelE_density.py | scaffold/panel builder superseded by figures.R + export_fig5E_data.py |
| fastas/pairs/lrr_serkcik/panelF_density.png | scaffold/panel builder superseded by figures.R + export_fig5E_data.py |
| fastas/pairs/lrr_serkcik/make_crossspecies_final4_figure.py | pre-correction baseline (isolated-window embeddings); superseded by serkcik_nn_results_corrected.pkl |
| fastas/pairs/lrr_serkcik/serkcik_crossspecies_final4_density.png | pre-correction baseline (isolated-window embeddings); superseded by serkcik_nn_results_corrected.pkl |
| fastas/pairs/lrr_serkcik/serkcik_crossspecies_final5_density.png | pre-correction baseline (isolated-window embeddings); superseded by serkcik_nn_results_corrected.pkl |
| fastas/pairs/lrr_serkcik/make_crossspecies_nn_figure.py | pre-correction baseline (isolated-window embeddings); superseded by serkcik_nn_results_corrected.pkl |
| fastas/pairs/lrr_serkcik/serkcik_crossspecies_nn_density.png | pre-correction baseline (isolated-window embeddings); superseded by serkcik_nn_results_corrected.pkl |
| fastas/pairs/lrr_serkcik/nn_blast_results.pkl | pre-correction baseline (isolated-window embeddings); superseded by serkcik_nn_results_corrected.pkl |
| fastas/pairs/lrr_serkcik/nn_cosine_results.pkl | pre-correction baseline (isolated-window embeddings); superseded by serkcik_nn_results_corrected.pkl |
| fastas/pairs/lrr_serkcik/make_serkcik_figure.py | v1 27-gene LOO figure, superseded |
| fastas/pairs/lrr_serkcik/serkcik_loo_density.png | v1 27-gene LOO figure, superseded |
| fastas/pairs/lrr_serkcik/serkcik_loo_results.pkl | v1 27-gene LOO figure, superseded |
| fastas/pairs/lrr_serkcik/make_hardneg_slope_figure.py | exploratory panel not in current figure |
| fastas/pairs/lrr_serkcik/serkcik_hardneg_slope.png | exploratory panel not in current figure |
| fastas/pairs/lrr_serkcik/serkcik_hardneg_slope.pkl | exploratory panel not in current figure |
| fastas/pairs/lrr_serkcik/make_2class_blast_heatmap.py | exploratory panel not in current figure |
| fastas/pairs/lrr_serkcik/serkcik_2class_blast_heatmap.png | exploratory panel not in current figure |
| fastas/pairs/lrr_serkcik/serkcik_2class_blast_matrix.pkl | exploratory panel not in current figure |
| fastas/pairs/lrr_serkcik/make_family_and_clustered_heatmaps.py | exploratory panel not in current figure |
| fastas/pairs/lrr_serkcik/serkcik_blast_clustered_heatmap.png | exploratory panel not in current figure |
| fastas/pairs/lrr_serkcik/serkcik_clefamily_blast_heatmap.png | exploratory panel not in current figure |
| fastas/pairs/lrr_serkcik/make_singleton_heatmap.py | exploratory panel not in current figure |
| fastas/pairs/lrr_serkcik/serkcik_singleton_vs_all_heatmap.png | exploratory panel not in current figure |
| fastas/pairs/lrr_serkcik/serkcik_singleton_vs_all_matrix.pkl | exploratory panel not in current figure |
| fastas/pairs/lrr_serkcik/serkcik_singleton_meta.pkl | exploratory panel not in current figure |
| fastas/pairs/lrr_serkcik/serkcik_all24.fasta | 24-gene set superseded by all27/expanded42 |
| fastas/pairs/lrr_serkcik/serkcik_all24.fasta.pkl | 24-gene set superseded by all27/expanded42 |
| fastas/pairs/lrr_serkcik/serkcik_all24.fasta.pkl.description | 24-gene set superseded by all27/expanded42 |
| fastas/pairs/lrr_serkcik/serkcik_all24.fasta.pkl.seqnames | 24-gene set superseded by all27/expanded42 |
| fastas/pairs/lrr_serkcik/train35_vs_train35_blast.tsv | different BLAST run from panel D; using it silently flips PSKR1 (see memory gotcha) |
| fastas/pairs/lrr_serkcik/export_fig5E_data.17pos.py.bak | 17-positive export, superseded 2026-09-28 (data archived in figures/figure_data/archive_fig5E_17pos_20260928) |
| fastas/pairs/myc_mad/MYC_concept | April-May 2026 PCA-128 own-scaler generation (366-407 span); superseded by shared-scaler MYC_HUMAN.fasta_concept |
| fastas/pairs/myc_mad/MAD1_concept | April-May 2026 PCA-128 own-scaler generation (366-407 span); superseded by shared-scaler MYC_HUMAN.fasta_concept |
| fastas/pairs/myc_mad/whole_bHLH_results | April-May 2026 PCA-128 own-scaler generation (366-407 span); superseded by shared-scaler MYC_HUMAN.fasta_concept |
| fastas/pairs/myc_mad/top_proteins_MYC.tsv | April-May 2026 PCA-128 own-scaler generation (366-407 span); superseded by shared-scaler MYC_HUMAN.fasta_concept |
| fastas/pairs/myc_mad/top_proteins_MAD1.tsv | April-May 2026 PCA-128 own-scaler generation (366-407 span); superseded by shared-scaler MYC_HUMAN.fasta_concept |
| fastas/pairs/myc_mad/validation_MYC_to_MAD1.tsv | April-May 2026 PCA-128 own-scaler generation (366-407 span); superseded by shared-scaler MYC_HUMAN.fasta_concept |
| fastas/pairs/myc_mad/validation_MAD1_to_MYC.tsv | April-May 2026 PCA-128 own-scaler generation (366-407 span); superseded by shared-scaler MYC_HUMAN.fasta_concept |
| fastas/pairs/myc_mad/MYC_HUMAN.fasta.spans.npy | April-May 2026 PCA-128 own-scaler generation (366-407 span); superseded by shared-scaler MYC_HUMAN.fasta_concept |
| fastas/pairs/myc_mad/MAD1_HUMAN.fasta.spans.npy | April-May 2026 PCA-128 own-scaler generation (366-407 span); superseded by shared-scaler MYC_HUMAN.fasta_concept |
| fastas/pairs/myc_mad/MYC_HUMAN.fasta.spans.prev | April-May 2026 PCA-128 own-scaler generation (366-407 span); superseded by shared-scaler MYC_HUMAN.fasta_concept |
| fastas/pairs/myc_mad/neg.npy | April-May 2026 PCA-128 own-scaler generation (366-407 span); superseded by shared-scaler MYC_HUMAN.fasta_concept |
| fastas/pairs/myc_mad/MYC_HUMAN.fasta.bHLH.spans | alternative MYC/MAD1 span definition (bHLH 369-421 / DNA-contact residues); paper uses 366-420 |
| fastas/pairs/myc_mad/MYC_HUMAN.fasta.bHLH_concept | alternative MYC/MAD1 span definition (bHLH 369-421 / DNA-contact residues); paper uses 366-420 |
| fastas/pairs/myc_mad/MYC_HUMAN.fasta.bHLH_pos.npy | alternative MYC/MAD1 span definition (bHLH 369-421 / DNA-contact residues); paper uses 366-420 |
| fastas/pairs/myc_mad/MYC_HUMAN.fasta.bHLH_query_results.tsv | alternative MYC/MAD1 span definition (bHLH 369-421 / DNA-contact residues); paper uses 366-420 |
| fastas/pairs/myc_mad/MYC_HUMAN.fasta.DNAcontact.span | alternative MYC/MAD1 span definition (bHLH 369-421 / DNA-contact residues); paper uses 366-420 |
| fastas/pairs/myc_mad/MYC_HUMAN.fasta.DNAcontact.span_concept | alternative MYC/MAD1 span definition (bHLH 369-421 / DNA-contact residues); paper uses 366-420 |
| fastas/pairs/myc_mad/MYC_HUMAN.fasta.DNAcontact.span_pos.npy | alternative MYC/MAD1 span definition (bHLH 369-421 / DNA-contact residues); paper uses 366-420 |
| fastas/pairs/myc_mad/MYC_HUMAN.fasta.DNAcontact.span_query_results.tsv | alternative MYC/MAD1 span definition (bHLH 369-421 / DNA-contact residues); paper uses 366-420 |
| fastas/pairs/myc_mad/MAD1_HUMAN.fasta.bHLH.spans | alternative MYC/MAD1 span definition (bHLH 369-421 / DNA-contact residues); paper uses 366-420 |
| fastas/pairs/myc_mad/MAD1_HUMAN.fasta.DNAcontact.span | alternative MYC/MAD1 span definition (bHLH 369-421 / DNA-contact residues); paper uses 366-420 |
