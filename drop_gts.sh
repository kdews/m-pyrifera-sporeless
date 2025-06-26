#!/bin/bash
#SBATCH -p gpu
#SBATCH --cpus-per-task=32
#SBATCH --mem=100gb
#SBATCH --time=05:00:00
#SBATCH -J drop_gts
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

# Input VCF file
in_vcf="$1"
vcf_basename="$(basename "$in_vcf")"
vcf_base="${vcf_basename%%.*}"
vcf_ext="${vcf_basename#*.}"
# Output VCF named from input
out_vcf="${vcf_base}.info_only.${vcf_ext}"

# Save current directory path
CURDIR="$(pwd)"
echo "Current directory: $CURDIR"
# Create temporary working directory in /tmp filesystem
# (Avoids mmap errors with BLAST in /scratch1 and /project)
WORKDIR="/tmp/${USER}/${SLURM_JOB_NAME}"
mkdir -p "$WORKDIR"
cd "$WORKDIR" || exit 1
echo "Temporary working directory: $WORKDIR"
echo

# Copy input to working directory
if [[ -f "${CURDIR}/${in_vcf}" ]]
then
    echo "Copying input files to $WORKDIR"
    # VCF
    cmd="rsync -htu --progress ${CURDIR}/${in_vcf} ."
    echo "$cmd"
    $cmd
else
    echo "Error: Input VCF file ($in_vcf) not found."
    exit 1
fi
echo

# Load BCFtools env
conda activate bcftools
bcftools --version
echo

# Cut first 9 columns from tab-delimited VCF file
echo "Dropping genotypes from $in_vcf..."
cmd="bcftools view --threads $SLURM_CPUS_PER_TASK -G -o $out_vcf $vcf_basename"
echo "$cmd"
$cmd

# Copy output from /tmp to current directory
# Preserve timestamps (-t)
# Only copy/update files if newer (-u)
echo "Copying output in $WORKDIR to $CURDIR"
rsync -htu --progress "$WORKDIR/"* "$CURDIR/"
