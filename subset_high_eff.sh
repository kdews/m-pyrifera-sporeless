#!/bin/bash
#SBATCH -J subset_high_eff
#SBATCH --mem=60gb
#SBATCH --cpus-per-task=8
#SBATCH --time=01:00:00
#SBATCH -o %x_%j.log

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

# Input
subset_tab="high_eff_subset.tsv"
in_vcf="raw_haploid_559_indv_on_CI_03_genome_final.ann.vcf.gz"

# Output
# Name from input
subset_basename="$(basename "$subset_tab")"
subset_base="${subset_basename%%.*}"
vcf_basename="$(basename "$in_vcf")"
vcf_base="${vcf_basename%%.*}"
vcf_ext="${vcf_basename#*.}"
subset_vcf="${subset_base}_${vcf_base}.$vcf_ext"

# Multithread support
if [[ -n "$SLURM_CPUS_PER_TASK" ]]
then
    nthrd="$SLURM_CPUS_PER_TASK"
else
    nthrd="1"
fi

# Load BCFtools env
conda activate bcftools
bcftools --version
echo

# Filter for sites in table
cmd=(
  bcftools view
  --threads "$nthrd"
  -R "$subset_tab"
  -Oz -o "$subset_vcf"
  "$vcf_basename"
)
echo "${cmd[*]}"
"${cmd[@]}"

# Done
echo
echo "Finished at:"
date
