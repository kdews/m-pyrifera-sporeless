# Mutation annotation of the brown macroalgae *Macrocystis pyrifera* (giant kelp)
A pipeline to detect and annotate the effects of deleterious mutations in *M. pyrifera* whole genome sequencing data.

This pipeline was created using the *Macrocystis pyrifera* CI_03 v1.0 genome[<sup>1</sup>](https://doi.org/10.1186/s12864-023-09658-x) and the asscoiated Joint Genome Institute (JGI) [PhycoCosm entry for *M. pyrifera*](https://mycocosm.jgi.doe.gov/Macpyr2/Macpyr2.home.html).

1. Diesel J, Molano G, Montecinos GJ, DeWeese K, Calhoun S, Kuo A, Lipzen A, Salamov A, Grigoriev IV, Reed DC, Miller RJ, Nuzhdin SV, Alberto F. A scaffolded and annotated reference genome of giant kelp (*Macrocystis pyrifera*). BMC Genomics. 2023 Sep 13;24(1):543. doi: [10.1186/s12864-023-09658-x](https://doi.org/10.1186/s12864-023-09658-x)

## Pipeline
### 1. Index VCF file and extract basic metadata
```
sbatch vcf_prep.sh <input>.vcf(.gz)
```
##### Output files
- `raw_vcf_header.txt`: header of input VCF
- `raw_vcf_ids.txt`: list of raw sample IDs as they appear in input VCF
- `<input>.csi`: CSI-format index of input VCF

#### Optional: Calculate variant statistics
```
sbatch bcftools_stats.sh <input>.vcf(.gz)
```

### 2. Run SnpEff to annotate predicted variant impacts
```
sbatch snpEff.sh <input>.vcf(.gz)
```
##### Output files
- `<input>.ann.vcf`: snpEff-annotated VCF file

### 3. Compress and index resulting annotated VCF file
```
sbatch compr_idx.sh <input>.ann.vcf
```
##### Output files
- `<input>.ann.vcf.gz`: compressed snpEff-annotated VCF
- `<input>.ann.vcf.gz.csi`: CSI-format index of snpEff-annotated VCF

### 4. Drop genotype information from snpEff-annotated VCF file
#### Keep only metdata for R analysis
```
sbatch drop_gts.sh <input>.ann.vcf.gz
```
##### Output files
- `<input>.info_only.ann.vcf.gz`: snpEff-annotated VCF (metadata only)

### 5. Analyze and filter snpEff variant annotations
```
sbatch metadata_analysis.sh <sample_metadata>.csv <input>.ann.vcf.gz <genes>.gff <gene_list>.csv <annotation_prefix>
```
##### Input files
- `<sample_metadata>.csv`: CSV of sample metadata containing at least the 2 titled columns "SampleID" and "GametophyteCode"
- `<input>.ann.vcf.gz`: snpEff-annotated VCF (metadata only)
- `<genes>.gff`: gene feature file (GFF) for *M. pyrifera* from [JGI PhycoCosm](https://mycocosm.jgi.doe.gov/Macpyr2/Macpyr2.home.html)
- `<gene_list>.csv`: CSV of meiosis genes identified on the [*M. pyrifera* JGI Genome Portal GUI](https://mycocosm.jgi.doe.gov/Macpyr2/Macpyr2.home.html)
- `<annotation_prefix>`: prefix for [JGI's functional gene annotation tables](https://genome.jgi.doe.gov/portal/Macpyr2/Macpyr2.download.html) for *M. pyrifera* (this pipeline used "Macpyr2_GeneCatalog_proteins_20220914")
##### Output files
###### Tables
- `sample_idx.tsv`: index of sample IDs
- `<input>.split.tsv.gz`: split VCF metadata table
- `<input>.split.col_types.tsv`: column types of split VCF metadata (improves R file reading speed)
- `debug_<input>.split.tsv.gz`: table of snpEff debugging messages
- `high_eff_subset.tsv`: 2-column (CHR POS) table to subset VCF
- `all_meiotic_proteins.tab`: table of meiotic protein GFF + GO/IPR/KEG/KOGG annotations
- `all_meiotic_protein_IDs.txt`: list of protein IDs implicated in meiosis
###### Plots
- `all_var_qc_boxplots.png`: variant quality control boxplots
![alt text](https://github.com/kdews/m-pyrifera-sporeless/blob/main/all_var_qc_boxplots.png)
- `all_eff_qc_boxplots.png`: effect quality control boxplots
![alt text](https://github.com/kdews/m-pyrifera-sporeless/blob/main/all_eff_qc_boxplots.png)
- `all_eff_bar.png`: barplot summary of all observed variant effects
![alt text](https://github.com/kdews/m-pyrifera-sporeless/blob/main/all_eff_bar.png)
- `all_impact_bar.png`: basic distribution plot of impacts
![alt text](https://github.com/kdews/m-pyrifera-sporeless/blob/main/all_impact_bar.png)
- `all_eff_impact_type_bar.png`: barplot of Variant Type vs Impact/Effect
![alt text](https://github.com/kdews/m-pyrifera-sporeless/blob/main/all_eff_impact_type_bar.png)
- `debug_snpeff.png`: barplot of snpEff debugging messages
![alt text](https://github.com/kdews/m-pyrifera-sporeless/blob/main/debug_snpeff.png)

### 6. Subset high impact variants in meiosis-associated genes
```
sbatch subset_high_eff.sh <input>.ann.vcf.gz
```
##### Output files
- `high_eff_subset_<input>.ann.vcf.gz`: subset annotated VCF

### 7. Generate crossing plan based on variant impact predictions
```
sbatch predict_sporeless.sh high_eff_subset_<input>.ann.vcf.gz
```
##### Output files
- `high_eff_subset_<input>.ann.split.tsv.gz`: table of putative high impact meiotic variants (split annotated VCF) with genotypes
- `high_eff_subset_<input>.ann.prot.split.tsv.gz`: same table as above with protein annotation columns and rows ranked by variant severity and relevance to meiosis
- `cross_table.tsv`: detailed table of predicted non-reproductive crosses
- `cross_list.tsv`: highly simplified table (Rank, Cross_ID, F_genotype, M_genotype) of predicted non-reproductive crosses
