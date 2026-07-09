#Prepare pooled embedding vectors from a .pkl file for CAV training.

#Reads a .pkl + .pkl.info embedding file and a spans file, pools each sample's
#embeddings according to its span, and saves the result as a .npy file.

#Run this separately for positive and negative sets, then pass the resulting
#.npy files to train_cav_from_embeddings.py.

#Spans file format (tab-separated, no header, '#' lines are comments):
#    accession                       # whole-sequence
#    accession  TAB  pos             # single position (e.g. PTM)
#    accession  TAB  start  TAB  end # window [start, end)

#Examples
#--------
# Motif positives (window spans):
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
source "$(dirname "$0")/../../config/paths.sh"
workdir=$BIOCAV_PAPER_ROOT/protein_analysis
for f in $workdir/fastas/pairs/myc_mad/*fasta
do
echo $f
python $BIOCAV_REPO/scripts/prepare_embeddings.py \
    --pkl ${f}.pkl \
    --info ${f}.pkl.seqnames \
    --spans ${f}.spans \
    --out ${f}.spans.npy
done


