#!/bin/bash
for sample in *_aa_scaffolds.fas
do
    grep -c ">" ${sample} > ${sample}_gene_count.txt
done

for file in *_gene_count.txt;
do
    echo "$file"$'\t'"$(cat -- "$file")" > "$file"
done

cat *_gene_count.txt > mespa_gene_count_all.txt

rm *_gene_count.txt