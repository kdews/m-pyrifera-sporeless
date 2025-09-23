#!/bin/bash
#SBATCH -p gpu
#SBATCH --mem=180gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=24
#SBATCH -J metadata_analysis
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
# High impact variants (subset of annotated VCF)
# sub_vcf_file="high_eff_subset_raw_haploid_559_indv_on_CI_03_genome_final.ann.vcf.gz"
sub_vcf_file="$1"

# Run Rscript file
cmd=(
    "Rscript"
    "$rscript"
    "$sub_vcf_file"
    "$scripts_dir"
)
echo "${cmd[*]}"
"${cmd[@]}"

# Print time at end
echo
date
