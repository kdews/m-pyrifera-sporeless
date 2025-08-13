#!/bin/bash
#SBATCH -p gpu
#SBATCH --mem=180gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=24
#SBATCH -J sporeless_gk
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
rscript="sporeless_gk.R"

# Input
# Scripts directory
scripts_dir="m-pyrifera-sporeless/"
rscript="${scripts_dir}${rscript}"
# Table of individual metadata
meta_file="070721_metadata.csv"
# VCF ID file
vcf_id_file="raw_vcf_ids.txt"
# Annotated VCF
ann_vcf_file="raw_haploid_559_indv_on_CI_03_genome_final.info_only.ann.vcf.gz"
# Gene feature file
gff_file="genes.gff"
# CSV file of meiosis-related genes from JGI "Search" GUI
meiotic_gene_list_file="jgi_gui_meiotic_genes.csv"

# Run Rscript file
cmd=(
    "Rscript"
    "$rscript"
    "$meta_file"
    "$vcf_id_file"
    "$ann_vcf_file"
    "$gff_file"
    "$meiotic_gene_list_file"
    "$scripts_dir"
)
echo "${cmd[*]}"
"${cmd[@]}"

# Print time at end
echo
date
