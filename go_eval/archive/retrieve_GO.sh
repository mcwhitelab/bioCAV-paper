GO=$1
echo $GO
#GO_0000018

cd static/GO_sets

# Download and extract only if directory doesn't exist
if [ ! -d "$GO" ]; then
    if [ ! -f "${GO}.zip" ]; then
        wget https://profab.kansil.org/go_dataset/${GO}.zip
    fi
    unzip ${GO}.zip
    rm ${GO}.zip
else
    echo "Directory $GO already exists, skipping download and extraction"
fi

cd $GO

# Generate .ids files only if they don't exist
for f in temporal*paac*; do
    if [ ! -f "${f}.ids" ]; then
        awk -F'\t' '{print $1}' $f > ${f}.ids
        echo "Created ${f}.ids"
    else
        echo "Skipping ${f}.ids, already exists"
    fi
done

# Run create_files.py only if *_ids.list files don't already exist
if ! ls *_ids.list 1> /dev/null 2>&1; then
    python ../../../scripts/create_files.py
else
    echo "Skipping create_files.py, *_ids.list files already exist"
fi

# Retrieve FASTAs only if they don't exist
for f in *_ids.list; do
    out="${f%_ids.list}.fasta"
    if [ ! -f "$out" ]; then
        python ../../../scripts/retrieve_fasta.py $f $out --batch_size 100
    else
        echo "Skipping $out, already exists"
    fi
done

cd ../../../

for f in static/GO_sets/$GO/*fasta; do
    echo "Exists: $f"
done

for f in static/GO_sets/$GO/*t[re][as]*fasta; do
    python ../github/mcwlab_utils/hf_embed_new.py -f $f -o ${f}.pkl -ss mean -s -l  -9 -10 -11 -12   -m /groups/clairemcwhite/models/ESMplusplus_large -b 1 --max_length 2048 --subsample 200
done

for f in static/GO_sets/$GO/*valid*fasta; do
    python ../github/mcwlab_utils/hf_embed_new.py -f $f -o ${f}.pkl -a -l  -9 -10 -11 -12   -m /groups/clairemcwhite/models/ESMplusplus_large -b 1 --max_length 2048 --subsample 200
done




