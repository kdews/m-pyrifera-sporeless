#!/bin/bash
#SBATCH --mem-per-cpu=4gb
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=12
#SBATCH -J vcf_prep
#SBATCH -o %x.log

# Print date and time
date
echo

# Load conda
cond=~/.conda_for_sbatch.sh
if [[ -a "$cond" ]]
then
  source "$cond"
else
  echo "Error on source of $cond"
  exit 1
fi

# Raw input VCF
raw_vcf="raw_haploid_559_indv_on_CI_03_genome_final.vcf.gz"
vcf_header="raw_vcf_header.txt"
vcf_ids="raw_vcf_ids.txt"

# Load BCFtools env
conda activate bcftools
bcftools --version

# Save VCF header and list of IDs to file
cmd="bcftools head $raw_vcf"
echo "$cmd > $vcf_header"
$cmd > $vcf_header
cmd="bcftools query --list-samples $raw_vcf"
echo "$cmd > $vcf_ids"
$cmd > $vcf_ids

# Index raw VCF file
cmd="bcftools index --threads $SLURM_CPUS_PER_TASK $raw_vcf"
echo "$cmd"
$cmd

