#!/bin/bash
#SBATCH -J subset_mult
#SBATCH --mem-per-cpu=4gb
#SBATCH --cpus-per-task=12
#SBATCH --time=10:00:00
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
# Output subset VCFs
mult_snps="mult_snps.vcf.gz"
mult_indels="mult_indels.vcf.gz"

# Load BCFtools env
conda activate bcftools
bcftools --version

# Calculate variant stats before filtering
cmd1="bcftools +counts $raw_vcf"
echo "$cmd1"
$cmd1

# Keep only multiallelic SNPs
cmd2=(
  bcftools view
  --threads "$SLURM_CPUS_PER_TASK"
  -i 'TYPE="snp" && N_ALT>1'
  -Oz -o "$mult_snps"
  "$raw_vcf"
)
echo "${cmd2[*]}"
"${cmd2[@]}"

# Calculate variant stats after filtering
cmd3="bcftools +counts $mult_snps"
echo "$cmd3"
$cmd3

# Keep only multiallelic indels
cmd4=(
  bcftools view
  --threads "$SLURM_CPUS_PER_TASK"
  -i 'TYPE="indel" && N_ALT>1'
  -Oz -o "$mult_indels"
  "$raw_vcf"
)
echo "${cmd4[*]}"
"${cmd4[@]}"

# Calculate variant stats after filtering
cmd5="bcftools +counts $mult_indels"
echo "$cmd5"
$cmd5
