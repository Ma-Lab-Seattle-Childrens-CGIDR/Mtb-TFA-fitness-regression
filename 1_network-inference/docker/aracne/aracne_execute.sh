#!/bin/bash

# See https://github.com/califano-lab/ARACNe-AP?tab=readme-ov-file#aracne-ap for usage
# Running time is ~3m for 200 samples, 20531 genes, and 1813 transcription factors according to
# https://github.com/califano-lab/GPU-ARACNe (presumably for one bootstrap)

#
# Set parameters
#

expression_file="$1"
regulators_file="$2"
out_file="${3:-/root/mount/aracne_network.txt}"
transpose="${4:-n}"
p_value="${5:-1E-8}"
threads="${6:-4}"
bootstraps="${7:-10}"
memory_gbs="${8:-12}"

out_dir="$(dirname $out_file)/aracne_temp-$(date +%s)"
mkdir -p "$out_dir"

if [ "$transpose" = "y" ] || [ "$transpose" = "yes" ]
then
    datamash transpose < "$expression_file" > "$out_dir/expression_file.tsv"
    expression_file="$out_dir/expression_file.tsv"
fi

# ARACNe runtime can't handle carriage returns -- messes everything up silently
tr -d '\r' < "$expression_file" > temp_expression_file
cp temp_expression_file "$out_dir/expression_file.tsv" # mv often fails here in container; just copy
expression_file="$out_dir/expression_file.tsv"

tr -d '\r' < "$regulators_file" > "$out_dir/regulators_file.tsv"
regulators_file="$out_dir/regulators_file.tsv"

#
# Step 1: Calculate threshold (result is written to a file)
#

# seed must be provided with --calculateThreshold to avoid a NPE
java -Xmx${memory_gbs}G -jar /root/dist/aracne.jar \
    -e "$expression_file" \
    -o "$out_dir" \
    --tfs "$regulators_file" \
    --pvalue $p_value \
    --seed 1 \
    --calculateThreshold

#
# Step 2: Run bootstraps of the network (results are written to file, one per bootstrap)
#

for (( i = 0; i < $bootstraps; i++ ))
do
    java -Xmx${memory_gbs}G -jar /root/dist/aracne.jar \
        -e "$expression_file" \
        -o "$out_dir" \
        --tfs "$regulators_file" \
        --pvalue $p_value \
        --threads $threads
done

#
# Step 3: Consolidate bootstraps into a single result
#

java -Xmx${memory_gbs}G -jar /root/dist/aracne.jar \
    -o "$out_dir" \
    --consolidate

#
# Clean up
#

# sort by 4th field (p-val), with secondary sort / tie-break using 3rd field (mutual information);
# then select fields 1, 2, and 4 of results (tf, gene, p-val), remove header line, and save
sort -k 4,4g -k 3,3gr "$out_dir/network.txt" | \
    cut -f 1,2,4 | \
    tail -n +2 > "$out_file"
cp "$out_dir/network.txt" "${out_file}_details.txt"
rm -rf "$out_dir"
