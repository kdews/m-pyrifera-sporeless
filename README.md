# Mutation annotation of the brown macroalgae *Macrocystis pyrifera* (giant kelp)
A pipeline to detect and annotate the effects of deleterious mutations in *M. pyrifera* whole genome sequencing data.

## Pipeline
### 1. Index VCF file and extract basic metadata
```
sbatch vcf_prep.sh <in_vcf>
```
Output files: `raw_vcf_header.txt`, `raw_vcf_ids.txt`, and `<in_vcf>.csi`.
##### Optional: Calculate variant statistics
```
sbatch bcftools_stats.sh <in_vcf>
```

### 2. Run SnpEff to annotate predicted variant impacts
```
sbatch snpEff.sh <in_vcf>
```
Output file: `<in_vcf>.ann.vcf`

### 3. Compress and index resulting annotated VCF file
```
sbatch compr_idx.sh <in_vcf>
```
Output files: `<in_vcf>.ann.vcf.gz` and `<in_vcf>.ann.vcf.gz.csi`.

### 4. Drop genotype information from snpEff-annotated VCF file
#### Keep only metdata for R analysis
```
sbatch drop_gts.sh <in_vcf>.ann.vcf.gz
```

### 5. Run 
```
sbatch sporeless_gk.sh <metadata> <in_vcf>.ann.vcf.gz <gff> <gene_list_csv>
```
Output files:
- Sample ID index table: `simple_metadata.tsv`
- 
