#!/bin/bash
#SBATCH -p gpu
#SBATCH --mem=200gb
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=32
#SBATCH -J sporeless_gk
#SBATCH -o %x_%j.log

# Print date and time
date
echo

# Load R and bgzip modules
module purge
module --latest load gcc htslib
module load r/4.4.3
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

# Run Rscript file
cmd=(
    "Rscript"
    "$rscript"
    "$meta_file"
    "$vcf_id_file"
    "$ann_vcf_file"
    "$scripts_dir"
)
echo "${cmd[*]}"
"${cmd[@]}"

# Print time at end
echo
date
