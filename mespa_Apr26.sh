#!/bin/bash
#SBATCH --time=5:00:00
#SBATCH --mem=50GB
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --job-name=mespa
#SBATCH --mail-user=living.li@csiro.au

for sample in $(cat sample_names.txt)
do
    python /scratch1/li266/WGS/mespa/mespa_Feb28.py \
    -a /scratch1/li266/WGS/assembly/${sample}_spades.fasta \
    -p /scratch1/li266/WGS/mespa/TCAST_R_aa.fas \
    -s /scratch1/li266/WGS/mespa/mespa_wgs_Apr08.cfg \
    -c 4
    mv mespa_results/onlygenemodels.fa mespa_results/${sample}_onlygenemodels.fa
    mv mespa_results/order_scaff.txt mespa_results/${sample}_order_scaff.txt
    mv mespa_results/output_summary.txt mespa_results/${sample}_output_summary.txt
    mv mespa_results/protein_scaff_relationship.tab mespa_results/${sample}_protein_scaff_relationship.tab
    mv mespa_results/cov_table.tab mespa_results/${sample}_cov_table.tab
    rm mespa_results/assembly_edited.fa
    rm mespa_results/list_after_scaff.txt
    rm mespa_results/list_before_scaff.txt
    rm mespa_results/potential_duplicates.tab
    rm mespa_results/remaining_contigs.fa
    rm mespa_results/scaffolds_for_manual_curation.txt
    mv mespa_results/scaffolds.gff mespa_results/${sample}_scaffolds.gff
    mv mespa_results/scaffolds.mfa mespa_results/${sample}_scaffolds.mfa
    rm mespa_results/spaln_alignment_summary.tab
    rm mespa_results/spaln_assembly_edited.aln
    mv mespa_results/spaln_assembly_edited.gff mespa_results/${sample}_assembly.gff
    rm mespa_results/spaln_assembly_edited.mb
    rm mespa_results/spaln_scaffolds.aln
    rm mespa_results/spaln_scaffolds.mb
    rm mespa_results/spaln_translated_proteins_assembly_edited.fa
    mv mespa_results/spaln_translated_proteins_assembly_edited_reformatted.fa mespa_results/${sample}_aa_assembly.fas
    rm mespa_results/spaln_translated_proteins_scaffolds.fa
    mv mespa_results/spaln_translated_proteins_scaffolds_reformatted.fa mespa_results/${sample}_aa_scaffolds.fas
    rm mespa_results/timelog.txt
    rm mespa_results/validation_statistics.tab
    rm mespa_results/xeno_filtered.txt
done
