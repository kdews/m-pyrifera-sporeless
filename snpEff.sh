#!/bin/bash
#SBATCH --mem=150gb
#SBATCH --time=1-0
#SBATCH -J snpEff
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

# snpEff configuration
snpEff_path="/scratch1/kdeweese/giant_kelp_scratch/snpEff"
snpEff_path="$(realpath -e "$snpEff_path")"
snpEff="$snpEff_path/snpEff.jar"
cpath="/project/noujdine_61/mkovalev/snpEff/snpEff.config"
dpath="/project/noujdine_61/mkovalev/snpEff/data"
# Input
# Genome build
genome_base="Macpyr2"
# VCF file
in_vcf="$1"
in_vcf="$(realpath -e "$in_vcf")"
vcf_basename="$(basename "$in_vcf")"
vcf_base="${vcf_basename%%.*}"
out_csv="${vcf_base}.csv"
out_html="${vcf_base}.html"
out_vcf="${vcf_base}.ann.vcf"

# Default Java memory=8g, unless memory is specified by SLURM
if [[ -n "$SLURM_MEM_PER_NODE" ]]
then
    mem="$(( 85 * SLURM_MEM_PER_NODE / 100000 ))g"
else
    mem="8g"
fi
echo "Java memory set to $mem"

# Run snpEff annotation
if [[ -f "$in_vcf" ]]
then
	echo "Running snpEff on $in_vcf using $genome_base database build..."
	cmd=(
        java "-Xmx${mem}" -jar
        "$snpEff" -v
        -config "$cpath"
        -dataDir "$dpath"
        -csvStats "$out_csv"
        -s "$out_html"
        -no-upstream
        -no-downstream
        "$genome_base"
        "$in_vcf"
    )
    echo "${cmd[*]} > $out_vcf"
    "${cmd[@]}" > "$out_vcf"
else
	echo "Error: Input VCF file ($in_vcf) not found."
	exit 1
fi
