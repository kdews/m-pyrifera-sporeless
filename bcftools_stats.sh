#!/bin/bash
#SBATCH -J bcftools_stats
#SBATCH --mem=20gb
#SBATCH --time=10:00:00
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
in_vcf="$1"
if [[ -f "$in_vcf" ]]
then
  echo "Input VCF: $in_vcf"
else
  echo "Error: Input VCF file ($in_vcf) not found."
  exit 1
fi
echo

# Output
# Name from input
vcf_basename="$(basename "$in_vcf")"
vcf_base="${vcf_basename%%.vcf.*}"
out_stats="${vcf_base}.stats"

# Load BCFtools env
conda activate bcftools
bcftools --version
echo

# Run bcftools stats on input VCF
cmd=(
  bcftools stats "$in_vcf"
)
echo "${cmd[*]} > $out_stats"
"${cmd[@]}" > "$out_stats"

# Done
echo
echo "Finished at:"
date
