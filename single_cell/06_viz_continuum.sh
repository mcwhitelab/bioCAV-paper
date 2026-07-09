




export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
LIB=cav_library/b617ee1b-f8c8-4de9-b82b-e803ab93550d

python "$(dirname "$0")/scripts/cav_continuum_viz.py" \
    --coords        $LIB/results/hierarchy/cell_coordinates.tsv \
    --h5ad          $LIB/data/cells.h5ad \
    --gene-corr-dir $LIB/results/gene_correlation_celltype_tissue/ \
    --group-col     cell_type \
    --context-col   tissue \
    --condition-col disease \
    --level         L2 \
    --pairs \
        endothelial_cell__skin_epidermis__normal_vs_melanoma \
        epithelial_cell__skin_epidermis__normal_vs_melanoma \
        fibroblast__skin_epidermis__normal_vs_melanoma \
        T_cell__skin_epidermis__normal_vs_melanoma \
    --title "Melanoma vs. normal — transcriptional continuum (L2)" \
    --out   $LIB/results/figures/cav_continuum_melanoma.png \
    --top-n-each 5 \
--style strips 


python "$(dirname "$0")/scripts/cav_continuum_viz.py" \
    --coords        $LIB/results/hierarchy/cell_coordinates.tsv \
    --h5ad          $LIB/data/cells.h5ad \
    --gene-corr-dir $LIB/results/gene_correlation_celltype_tissue/ \
    --group-col     cell_type \
    --context-col   tissue \
    --condition-col disease \
    --level         L2 \
    --pairs \
       B_cell__lung__normal_vs_lung_cancer \
       endothelial_cell__lung__normal_vs_lung_cancer \
       epithelial_cell__lung__normal_vs_lung_cancer \
       fibroblast__lung__normal_vs_lung_cancer \
       mast_cell__lung__normal_vs_lung_cancer \
       mononuclear_phagocyte__lung__normal_vs_lung_cancer \
       neutrophil__lung__normal_vs_lung_cancer \
       plasmacytoid_dendritic_cell__lung__normal_vs_lung_cancer \
       T_cell__lung__normal_vs_lung_cancer \
    --title "Lung cancer vs. normal — transcriptional continuum (L2)" \
    --out   $LIB/results/figures/cav_continuum_lungcancer.png \
--top-n-each 4 \
--style strips 


python "$(dirname "$0")/scripts/cav_continuum_viz.py" \
    --coords        $LIB/results/hierarchy/cell_coordinates.tsv \
    --h5ad          $LIB/data/cells.h5ad \
    --gene-corr-dir $LIB/results/gene_correlation_celltype_tissue/ \
    --group-col     cell_type \
    --context-col   tissue \
    --condition-col disease \
    --level         L2 \
    --pairs \
B_cell__breast__normal_vs_breast_cancer \
mononuclear_phagocyte__breast__normal_vs_breast_cancer \
neutrophil__breast__normal_vs_breast_cancer \
T_cell__breast__normal_vs_breast_cancer \
    --title "Breast cancer vs. normal — transcriptional continuum (L2)" \
    --out   $LIB/results/figures/cav_continuum_breastcancer.png \
--top-n-each 4 \
--style strips 



python "$(dirname "$0")/scripts/cav_continuum_viz.py" \
    --coords        $LIB/results/hierarchy/cell_coordinates.tsv \
    --h5ad          $LIB/data/cells.h5ad \
    --gene-corr-dir $LIB/results/gene_correlation_celltype_tissue/ \
    --group-col     cell_type \
    --context-col   tissue \
    --condition-col disease \
    --level         L2 \
    --pairs \
BB_cell__colorectum__normal_vs_colorectal_cancer \
endothelial_cell__colorectum__normal_vs_colorectal_cancer \
epithelial_cell__colorectum__normal_vs_colorectal_cancer \
fibroblast__colorectum__normal_vs_colorectal_cancer \
mast_cell__colorectum__normal_vs_colorectal_cancer \
mononuclear_phagocyte__colorectum__normal_vs_colorectal_cancer \
T_cell__colorectum__normal_vs_colorectal_cancer \
    --title "Colorectal cancer vs. normal — transcriptional continuum (L2)" \
    --out   $LIB/results/figures/cav_continuum_colorectalcancer.png \
--top-n-each 4 \
--style strips 






