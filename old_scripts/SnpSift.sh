#!/bin/bash
#SBATCH --mem=150gb
#SBATCH --time=01:00:00
#SBATCH -J SnpSift
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
echo

# SnpSift configuration
snpEff_path="/scratch1/kdeweese/giant_kelp_scratch/snpEff"
SnpSift="$(realpath -e "$snpEff_path/SnpSift.jar")"

# Inputs
in_vcf="$1"
gene_id_file="$2"

in_vcf="$(realpath -e "$in_vcf")"
gene_id_file="$(realpath -e "$gene_id_file")"
vcf_basename="$(basename "$in_vcf")"
vcf_base="${vcf_basename%%.vcf.*}"
out_vcf="${vcf_base}.gene_filt.vcf"

# Output filename (in CURDIR)
CURDIR="$(pwd)"

# Default Java memory=8g, unless memory is specified by SLURM
if [[ -n "$SLURM_MEM_PER_NODE" ]]
then
    mem="$(( 85 * SLURM_MEM_PER_NODE / 100000 ))g"
else
    mem="8g"
fi
echo "Java memory set to $mem"
echo

# Make /tmp workdir
WORKDIR="/tmp/$USER"
mkdir -p "$WORKDIR"
cd "$WORKDIR" || exit 1
echo "Running in temporary directory: $WORKDIR"
echo

# Copy input VCF and index to /tmp
echo "Copying input VCF and tabix index (.tbi) to $WORKDIR..."
rsync -htu --progress "$in_vcf" .
rsync -htu --progress "${in_vcf}.tbi" .

# Run SnpSift to filter for gene IDs in file
if [[ -f "$vcf_basename" && -f "$gene_id_file" ]]
then
    echo "Filtering for $gene_id_file in $vcf_basename using SnpSift..."
    cmd=(
        java "-Xmx${mem}" -jar
        "$SnpSift" filter -v
        # --set "$gene_id_file"
        # "ANN[*].GENEID in SET[0]"
        "( ANN[*].GENEID = 'gene_15764' )"
        "$vcf_basename"
    )
    echo "${cmd[*]} > $out_vcf"
    "${cmd[@]}" > "$out_vcf"
else
    echo "Error: Input VCF ($vcf_basename) or gene ID list ($gene_id_file) not found."
    exit 1
fi

# Copy output from /tmp to current directory, preserving timestamps (-t)
echo "Copying output in $WORKDIR/ to $CURDIR/"
rsync -ht --progress "$out_vcf" "$CURDIR/"
# Clean up /tmp
if [[ -d "$WORKDIR" ]]; then
    echo "Cleaning up $WORKDIR"
    rm -rfv "$WORKDIR"
fi

# Done
echo
echo "Finished at:"
date
