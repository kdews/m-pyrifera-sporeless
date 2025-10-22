#!/bin/bash
#SBATCH -p gpu
#SBATCH -J metadata_analysis
#SBATCH --mem=180gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=24
#SBATCH -o %x_%j.log

# Print date and time
date
echo

# Load R and bgzip modules
module purge
module load apptainer/1.4.1 rstats/4.5.1
Rscript --version
echo

# R script to run
rscript="metadata_analysis.R"

# Input
# Scripts directory
scripts_dir="m-pyrifera-sporeless/"
rscript="${scripts_dir}${rscript}"
# Table of individual metadata
# sample_meta_file="070721_metadata.csv"
sample_meta_file="$1"
# snpEff-annotated VCF (metadata only)
# ann_vcf_file="raw_haploid_559_indv_on_CI_03_genome_final.info_only.ann.vcf.gz"
ann_vcf_file="$2"
# Gene feature file
# gff_file="genes.gff"
gff_file="$3"
# CSV file of meiosis-related genes from JGI "Search" GUI
# meiotic_gene_list_file="jgi_gui_meiotic_genes.csv"
meiotic_gene_list_file="$4"
# Prefix for gene annotation table filenames
# gene_annot_tab_base="Macpyr2_GeneCatalog_proteins_20220914"
gene_annot_tab_base="$5"

# Run Rscript file
cmd=(
    "Rscript"
    "$rscript"
    "$sample_meta_file"
    "$ann_vcf_file"
    "$gff_file"
    "$meiotic_gene_list_file"
    "$gene_annot_tab_base"
    "$scripts_dir"
)
echo "${cmd[*]}"
"${cmd[@]}"

# Print time at end
echo
date
