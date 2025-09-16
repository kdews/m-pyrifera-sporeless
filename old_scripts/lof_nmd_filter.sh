#!/bin/bash
#SBATCH --mem=100gb
#SBATCH --time=05:00:00
#SBATCH --cpus-per-task=32
#SBATCH -J lof_nmd_filter
#SBATCH -o %x_%j.log

# Print date and time
date
echo

# Load open source Java
module purge
module --latest load gcc openjdk

# Test Java install (must be >= v1.8)
if java --version
then 
    echo "Java successfully loaded."
else
    echo "Error arose while testing Java install. Exiting..."
    exit 1
fi

# Default Java memory=8g, unless memory is specified by SLURM
if [[ -n "$SLURM_MEM_PER_NODE" ]]
then
    mem="$(( 85 * SLURM_MEM_PER_NODE / 100000 ))g"
else
    mem="8g"
fi
echo "Java memory set to $mem"
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

# Activate bcftools conda env
conda activate bcftools
bcftools --version
echo

# Multithread support for bcftools
if [[ -n "$SLURM_CPUS_PER_TASK" ]]
then
    nthrd="$SLURM_CPUS_PER_TASK"
else
    nthrd="1"
fi

# SnpSift configuration
snpEff_path="snpEff/"
snpEff_path="$(realpath -e "$snpEff_path")"
SnpSift="$snpEff_path/SnpSift.jar"
# Input
# VCF file
in_vcf="$1"
vcf_basename="$(basename "$in_vcf")"
vcf_base="${vcf_basename%%.*}"
out_vcf="$vcf_base.lof_nmd.ann.vcf.gz"

# Save current directory path
CURDIR="$(pwd)"
echo "Current directory: $CURDIR"
# Create temporary working directory in /tmp filesystem
# (Avoids mmap errors with BLAST in /scratch1 and /project)
WORKDIR="/tmp/$USER"
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
    if [[ -f "${CURDIR}/${in_vcf}.csi" ]]
    then
        # VCF index
        cmd="rsync -htu --progress ${CURDIR}/${in_vcf}.csi ."
        echo "$cmd"
        $cmd
    else
        echo "Error: Input VCF file index (${in_vcf}.csi) not found."
        exit 1
    fi
else
    echo "Error: Input VCF file ($in_vcf) not found."
    exit 1
fi
echo

# Run SnpSift filtration
if [[ -f "$vcf_basename" ]]
then
    echo "Filtering $vcf_basename with SnpSift..."
    cmd=(
        "java"
        "-Xmx${mem}"
        "-jar"
        "$SnpSift"
        "filter"
        "( LOF != '' ) | ( NMD != '' )"
        "$vcf_basename"
    )
    echo "${cmd[*]} | bgzip -@ $nthrd > $out_vcf"
    "${cmd[@]}" | bgzip -@ "$nthrd" > "$out_vcf"
    echo
else
    echo "Error: $vcf_basename not found in $WORKDIR"
fi
echo

# Index bgzipped VCF with bcftools
cmd=(
      "bcftools"
      "index"
      "--threads"
      "$nthrd"
      "$out_vcf"
)
echo "${cmd[*]}"
"${cmd[@]}"
echo

# Copy output from /tmp to current directory
# Preserve timestamps (-t)
# Only copy/update files if newer (-u)
echo "Copying output in $WORKDIR to $CURDIR"
rsync -htu --progress "$WORKDIR/"* "$CURDIR/"
# Clean up /tmp
if [[ -d "$WORKDIR" ]]; then
    echo "Cleaning up $WORKDIR"
    rm -rfv "$WORKDIR"
fi
