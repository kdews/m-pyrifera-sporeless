# Clear environment
rm(list = ls())
# Required packages
library(tidyverse)
library(vroom)
library(ggpubr)
library(ggrepel)
library(tidygraph)
library(ggraph)
if (require(showtext)) {
  showtext_auto()
  if (interactive())
    showtext_opts(dpi = 100)
  else
    showtext_opts(dpi = 300)
}

# Input
# Only take command line input if not running interactively
if (interactive()) {
  wd <- "/scratch1/kdeweese/giant_kelp_scratch"
  setwd(wd)
  # Table of individual metadata
  meta_file <- "070721_metadata.csv"
  # VCF ID file
  vcf_id_file <- "raw_vcf_ids.txt"
  # Annotated VCF (metadata only)
  ann_vcf_meta_file <-
    "raw_haploid_559_indv_on_CI_03_genome_final.info_only.ann.vcf.gz"
  # Gene feature file
  gff_file <- "genes.gff"
  # CSV of meiotic genes from JGI GUI
  meiotic_gene_list_file <- "jgi_gui_meiotic_genes.csv"
  # Output directory
  outdir <- "m-pyrifera-sporeless/"
} else {
  line_args <- commandArgs(trailingOnly = T)
  meta_file <- line_args[1]
  vcf_id_file <- line_args[2]
  ann_vcf_meta_file <- line_args[3]
  gff_file <- line_args[4]
  meiotic_gene_list_file <- line_args[5]
  outdir <- line_args[6]
}
go_tab_file <- "Macpyr2_GeneCatalog_proteins_20220914_GO.tab.gz"
ipr_tab_file <- "Macpyr2_GeneCatalog_proteins_20220914_IPR.tab.gz"
kegg_tab_file <- "Macpyr2_GeneCatalog_proteins_20220914_KEGG.tab.gz"
kog_tab_file <- "Macpyr2_GeneCatalog_proteins_20220914_KOG.tab.gz"
# Output
simple_meta_file <- "simple_metadata.tsv"
# Output table
ann_vcf_meta_base <-
  tools::file_path_sans_ext(basename(ann_vcf_meta_file), compression = T)
split_ann_vcf_meta_file <- paste0(ann_vcf_meta_base, ".split.tsv.gz")
split_ann_col_types_file <- paste0(ann_vcf_meta_base, ".split.col_types.tsv")
debug_df_file <- paste0("debug", "_", ann_vcf_meta_base, ".split.tsv.gz")
debug_plot <- "debug_snpeff.png"
var_qc_boxplot_plot <- "all_var_qc_boxplots.png"
eff_qc_boxplot_plot <- "all_eff_qc_boxplots.png"
all_eff_plot <- "all_eff_bar.png"
all_impact_plot <- "all_impact_bar.png"
top10_eff_plot <- "top10_eff_bar.png"
all_eff_impact_type_plot <- "all_eff_impact_type_bar.png"
meiotic_gene_list_base <- tools::file_path_sans_ext(
  basename(meiotic_gene_list_file)
)
meiotic_gene_bed_file <- paste0(meiotic_gene_list_base, ".bed")
all_meiotic_prot_ids_file <- "all_meiotic_protein_IDs.txt"
all_meiotic_prot_annot_file <- "all_meiotic_proteins.tab"
ann_vcf_base <- gsub("\\..*", "", basename(ann_vcf_meta_file))
subset_tab_base <- "high_eff_subset"
subset_tab_file <- paste0(subset_tab_base, ".tsv")
sub_vcf_file <- paste0(subset_tab_base, "_", ann_vcf_base, ".ann.vcf.gz")
split_sub_vcf_file <- paste0(
  subset_tab_base, "_", ann_vcf_base, ".split.tsv.gz"
)
# If exists, prepend output directory to output filenames of plots
if (dir.exists(outdir)) {
  debug_plot <- paste0(outdir, debug_plot)
  var_qc_boxplot_plot <- paste0(outdir, var_qc_boxplot_plot)
  eff_qc_boxplot_plot <- paste0(outdir, eff_qc_boxplot_plot)
  all_eff_plot <- paste0(outdir, all_eff_plot)
  all_impact_plot <- paste0(outdir, all_impact_plot)
  top10_eff_plot <- paste0(outdir, top10_eff_plot)
  all_eff_impact_type_plot <- paste0(outdir, all_eff_impact_type_plot)
}

# Standard vectors of column names
# GFF3
gff_names <-
  c("CHROM",
    "Source",
    "Type",
    "Start",
    "End",
    "Score",
    "Strand",
    "Phase",
    "Attributes")
# FORMAT column names
format_names <-
  c(
    # Genotype
    "GT",
    # Allelic Depth
    "AD",
    # Read Depth
    "DP",
    # Genotype Quality
    "GQ",
    # Phred-scaled Genotype Likelihood
    "PL"
  )
# INFO column patterns
info_names <-
  c(
    # Allele count in genotypes, for each ALT allele, in same order as listed
    "AC",
    # Allele Frequency, for each ALT allele, in the same order as listed
    "AF",
    # Total number of alleles in called genotypes
    "AN",
    # Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities
    "BaseQRankSum",
    # Approximate read depth; some reads may have been filtered
    "DP",
    # Were any of the samples downsampled?
    "DS",
    # Stop position of the interval
    "END",
    # Phred-scaled p-value for exact test of excess heterozygosity
    "ExcessHet",
    # Phred-scaled p-value using Fisher's exact test to detect strand bias
    "FS",
    # Inbreeding coefficient as estimated from the genotype likelihoods
    # per-sample when compared against the Hardy-Weinberg expectation
    "InbreedingCoeff",
    # Maximum likelihood expectation (MLE) for the allele counts
    # (not necessarily the same as the AC), for each ALT allele,
    # in the same order as listed
    "MLEAC",
    # Maximum likelihood expectation (MLE) for the allele frequency
    # (not necessarily the same as the AF), for each ALT allele,
    # in the same order as listed
    "MLEAF",
    # RMS Mapping Quality
    "MQ",
    # Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities
    "MQRankSum",
    # Variant Confidence/Quality by Depth
    "QD",
    # Raw data (sum of squared MQ and total depth) for improved RMS Mapping
    # Quality calculation. Incompatible with deprecated RAW_MQ formulation.
    "RAW_MQandDP",
    # Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias
    "ReadPosRankSum",
    # Symmetric Odds Ratio of 2x2 contingency table to detect strand bias
    "SOR",
    # Functional annotations: 'Allele | Annotation | Annotation_Impact |
    # Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType |
    # Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length |
    # AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO'
    "ANN",
    # Predicted loss of function effects for this variant.
    # Format: 'Gene_Name | Gene_ID | Number_transcripts | Percent_affected'
    "LOF",
    # Predicted nonsense mediated decay effects for this variant.
    # Format: 'Gene_Name | Gene_ID | Number_transcripts | Percent_affected'
    "NMD"
  )
# ANN column
ann_names <-
  c(
    "Allele",
    "Effect",
    "Impact",
    "Protein_ID",
    "Gene_ID",
    "Feature_Type",
    "Feature_ID",
    "Transcript_BioType",
    "Exon_or_Intron_Rank/Total",
    "HGVS.c",
    "HGVS.p",
    "cDNA_position/cDNA_len",
    "CDS_position/CDS_len",
    "Protein_position/Protein_len",
    "Dist_to_Feature",
    "Debug"
  )
# LOF column
lof_names <-
  c("Protein_ID",
    "Gene_ID",
    "Total_Transcripts_LOF",
    "Percent_Transcripts_LOF")
# NMD column
nmd_names <-
  c("Protein_ID",
    "Gene_ID",
    "Total_Transcripts_NMD",
    "Percent_Transcripts_NMD")
# LOF/NMD Effects
lof_nmd_effs <-
  c(
    "frameshift_variant",
    "stop_gained",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "start_lost"
  )
# Meiosis-associated search terms
meiotic_terms <- c(
  "meio",
  "reproduct",
  "double-strand break",
  " homologous recomb",
  "Spo11",
  "Rad51",
  "RecA",
  "DMC1",
  "Red1",
  "Hop1",
  "Smc5"
)

# Functions
# Selecting function for columns containing commas in any values
contains_commas <- function(df) {
  nms <- df %>%
    select(where(is.character)) %>%
    keep(~ any(replace_na(str_detect(.x, ","), FALSE))) %>%
    names()
  # Include ALT
  nms <- nms[nms %in% c(info_names, "ALT")]
  # Don't include snpEff colnames
  nms <- nms[!(nms %in% c("ANN", "LOF", "NMD"))]
}
# Fix long labels for plotting
pretty_labs <- function(my_vector) {
  my_vector %>%
    # Remove odd characters
    str_replace_all(c(
      "_variant" = "",
      "_" = " ",
      "&" = " & "
    )) %>%
    str_replace_all(regex("prime", ignore_case = T), "'")
  # # Wrap long labels
  # str_wrap(width = 20)
}

# Data analysis
# Save simplified table of sample IDs and gametophyte codes
if (!file.exists(simple_meta_file)) {
  # Read individual metadata
  metadata <- read_csv(meta_file)
  vcf_ids <- read_table(vcf_id_file, col_names = "vcfSampleID")
  vcf_ids <- vcf_ids %>%
    arrange(vcfSampleID) %>%
    pull(vcfSampleID)
  simple_meta <- metadata %>%
    select(SampleID, GametophyteCode) %>%
    unique() %>%
    arrange(SampleID)
  # Check metadata table IDs versus VCF IDs
  all(simple_meta$SampleID == vcf_ids)
  write_tsv(simple_meta, simple_meta_file)
}

# Parse meiotic gene annotations
if (!file.exists(all_meiotic_prot_annot_file)) {
  # Parse gene feature file
  gff <- vroom(gff_file, delim = "\t", comment = "#", col_names = gff_names)
  print(paste("GFF colnames:", paste(colnames(gff), collapse = " ")))
  # Separate INFO col into columns named by regular expression (e.g., AC=##)
  gff_split <- gff %>%
    mutate(Attributes = str_split(Attributes, ";") %>%
             map(~ {
               key_vals <- str_split_fixed(.x, "=", 2)
               set_names(key_vals[, 2], key_vals[, 1])
             })) %>%
    unnest_wider(Attributes, names_repair = "unique")
  gff_split <- gff_split %>%
    mutate(across(.cols = where(is.character), .fns = parse_guess)) %>%
    select(
      CHROM,
      Type,
      Start,
      End,
      Strand,
      Phase,
      ID,
      Protein_ID = proteinId,
      Transcript_ID = transcriptId,
      Parent,
      Protein_Product = product_name
    )
  print(paste(
    "GFF colnames (Attributes parsed):",
    paste(colnames(gff_split), collapse = " ")
  ))
  # Correlate gene IDs across all GFF lines
  gff_split <- gff_split %>%
    mutate(Parent = case_when(is.na(Parent) ~ ID,
                              .default = Parent)) %>%
    separate_wider_delim(
      Parent,
      delim = "_",
      names = c(NA, "Gene_ID", NA),
      too_few = "align_start",
      too_many = "drop",
      cols_remove = F
    ) %>%
    mutate(Gene_ID=paste0("gene_", Gene_ID)) %>%
    # Fill missing protein and transcript IDs using gene ID grouping
    group_by(Gene_ID) %>%
    fill(Protein_ID, Transcript_ID, Protein_Product, .direction = "downup")
  
  # Select genes of interest
  # Use CSV from JGI GUI search results
  if (!file.exists(meiotic_gene_bed_file)) {
    # Write BED file to subset VCF
    meiotic_gene_list <- read_csv(meiotic_gene_list_file)
    meiotic_gene_bed <- meiotic_gene_list %>%
      mutate(Strand = gsub("\\(|\\)", "", Strand)) %>%
      select(
        CHROM = Scaffold,
        Gene_Start = Start,
        Gene_End = End,
        Strand,
        Protein_ID = `Protein Id`,
        Transcript_ID = `Transcript Id`,
        Model_Name = Name
      )
    meiotic_gene_bed_export <- meiotic_gene_bed %>%
      rename("#CHROM" = CHROM)
    write_tsv(meiotic_gene_bed_export, meiotic_gene_bed_file, col_names = T)
  } else {
    meiotic_gene_bed <- read_tsv(meiotic_gene_bed_file, col_names = T)
    meiotic_gene_bed <- meiotic_gene_bed %>%
      rename(CHROM = "#CHROM")
  }
  # Filter JGI annotation tables with search terms
  meiotic_key <- paste(meiotic_terms, collapse = "|")
  go_tab <- read_tsv(go_tab_file) %>%
    # Remove "#" from colnames
    rename_with(~ gsub("#", "", .)) %>%
    # Collapse to one line per proteinId
    summarize(GO=paste(goName, collapse = ";"), .by = proteinId)
  ipr_tab <- read_tsv(ipr_tab_file, na = c("\\N")) %>%
    # Remove "#" from colnames
    rename_with(~ gsub("#", "", .)) %>%
    # Collapse to one line per proteinId
    summarize(IPR=paste(iprDesc, collapse = ";"), .by = proteinId)
  kegg_tab <- read_tsv(kegg_tab_file, na = c("\\N")) %>%
    # Remove "#" from colnames
    rename_with(~ gsub("#", "", .)) %>%
    # Collapse to one line per proteinId
    summarize(KEGG=paste(definition, collapse = ";"), .by = proteinId)
  kog_tab <- read_tsv(kog_tab_file, na = c("\\N")) %>%
    # Remove "#" from colnames
    rename_with(~ gsub("#", "", .)) %>%
    # Collapse to one line per proteinId
    summarize(KOG=paste(kogdefline, collapse = ";"), .by = proteinId)
  # Find protein IDs that match search terms
  all_meiotic_prot_ids <- rbind(
    go_tab %>%
      filter(grepl(meiotic_key, GO, ignore.case = T)) %>%
      distinct(proteinId),
    ipr_tab %>%
      filter(grepl(meiotic_key, IPR, ignore.case = T)) %>%
      distinct(proteinId),
    kegg_tab %>%
      filter(grepl(meiotic_key, KEGG, ignore.case = T)) %>%
      distinct(proteinId),
    kog_tab %>%
      filter(grepl(meiotic_key, KOG, ignore.case = T)) %>%
      distinct(proteinId),
    # Add previously identified protein IDs
    meiotic_gene_bed %>%
      select(proteinId = Protein_ID)
  ) %>%
    distinct() %>%
    pull(proteinId)
  # Filter annotation tables with meiotic protein IDs
  go_tab_filt <- go_tab %>%
    filter(proteinId %in% all_meiotic_prot_ids)
  ipr_tab_filt <- ipr_tab %>%
    filter(proteinId %in% all_meiotic_prot_ids)
  kegg_tab_filt <- kegg_tab %>%
    filter(proteinId %in% all_meiotic_prot_ids)
  kog_tab_filt <- kog_tab %>%
    filter(proteinId %in% all_meiotic_prot_ids)
  all_meiotic_prot <- go_tab_filt %>%
    full_join(ipr_tab_filt) %>%
    full_join(kegg_tab_filt) %>%
    full_join(kog_tab_filt) %>%
    rename(Protein_ID = proteinId)
  # Join annotations to GFF
  all_meiotic_prot_annot <- all_meiotic_prot %>%
    left_join(gff_split) %>%
    # Remove redundant lines (gene == mRNA lines)
    filter(Type != "mRNA") %>%
    mutate(Protein_ID=as.character(Protein_ID))
  # Write list of all protein IDs implicated in meiosis
  write.table(
    all_meiotic_prot_ids,
    all_meiotic_prot_ids_file,
    sep = "\t",
    row.names = F,
    col.names = F,
    quote = F
  )
  # Write table of meiotic protein GFF + GO/IPR/KEG/KOGG annotations
  write.table(
    all_meiotic_prot_annot,
    all_meiotic_prot_annot_file,
    sep = "\t",
    row.names = F
  )
} else {
  all_meiotic_prot_annot <- read_tsv(all_meiotic_prot_annot_file)
}

# Write 2-column table (CHR POS) to subset VCF
if (!file.exists(subset_tab_file)) {
  # Parse and QC annotated VCF metadata (without genotypes)
  if (!file.exists(split_ann_vcf_meta_file)) {
    # # Ensure bgzip is accessible
    # if (system("bgzip --version") != 0) {
    #   stop("bgzip is not available in this environment.")
    # }
    # # Use pipe("bgzip ...") to unzip large VCF file
    # ann_vcf_meta_raw <-
    #   vroom(pipe(paste(
    #     "bgzip -@ $SLURM_CPUS_PER_TASK -d -c", ann_vcf_meta_file
    #   )), delim = "\t", comment = "##")
    ann_vcf_meta_raw <- vroom(ann_vcf_meta_file, delim = "\t", comment = "##")
    # If interactive, keep random subset of annotated VCF + LOF/NMD lines
    # Otherwise, keep entire dataset in df
    if (interactive()) {
      set.seed(42)
      ann_vcf_meta <- ann_vcf_meta_raw %>% slice_sample(n = 50000)
    } else {
      ann_vcf_meta <- ann_vcf_meta_raw
    }
    ann_vcf_meta <- ann_vcf_meta %>%
      # Remove remaining "#" in colnames
      rename_with(~ gsub("#", "", .)) %>%
      # Drop irrelevant cols
      select(-c(ID, FILTER))
    print(paste("Colnames:", paste(colnames(ann_vcf_meta), collapse = " ")))
    # Separate INFO col into columns named by regular expression (e.g., AC=##)
    ann_vcf_meta <- ann_vcf_meta %>%
      mutate(INFO = str_split(INFO, ";") %>%
               map(~ {
                 key_vals <- str_split_fixed(.x, "=", 2)
                 set_names(key_vals[, 2], key_vals[, 1])
               })) %>%
      unnest_wider(INFO, names_repair = "unique")
    print(
      paste(
        "Colnames (INFO parsed):",
        paste(colnames(ann_vcf_meta), collapse = " ")
      )
    )
    # Strip special characters
    ann_vcf_meta <- ann_vcf_meta %>%
      # Remove parentheses from LOF and NMD cols
      mutate(across(c(LOF, NMD), ~ gsub("\\(|\\)", "", .))) %>%
      # Keep only protein ID numbers from "jgi.p|Macpyr2|####" for downstream
      mutate(
        ANN = gsub("jgi\\.p\\|Macpyr2\\|", "", ANN),
        LOF = gsub("jgi\\.p_Macpyr2_", "", LOF),
        NMD = gsub("jgi\\.p_Macpyr2_", "", NMD)
      )
    # Second reformat: Decompose VCF
    ann_vcf_meta <- ann_vcf_meta %>%
      # Separate into rows when multiple ALT per variant
      separate_longer_delim(cols = all_of(contains_commas(.)), delim = ",")
    # Third reformat: Split ANN col
    split_ann_vcf_meta <- ann_vcf_meta %>%
      # Separate into rows when multiple ANN per variant
      separate_longer_delim(ANN, delim = ",") %>%
      # Split ANN into named columns
      separate_wider_delim(ANN, delim = "|", names = ann_names) %>%
      # Keep only matching ANN annotations
      filter(ALT == Allele)
    # Fourth reformat: Split LOF & NMD cols
    lof_df <- ann_vcf_meta %>%
      filter(!is.na(LOF)) %>%
      select(CHROM, POS, ALT, LOF) %>%
      # Separate into additional rows when multiple per variant
      separate_longer_delim(LOF, delim = ",") %>%
      # Split LOF into named columns
      separate_wider_delim(LOF, delim = "|", names = lof_names)
    nmd_df <- ann_vcf_meta %>%
      filter(!is.na(NMD)) %>%
      select(CHROM, POS, ALT, NMD) %>%
      # Separate into additional rows when multiple per variant
      separate_longer_delim(NMD, delim = ",") %>%
      # Split NMD into named columns
      separate_wider_delim(NMD, delim = "|", names = nmd_names)
    # Join LOF and NMD back to ann_vcf_meta using Gene_ID and Protein_ID
    split_ann_vcf_meta <- split_ann_vcf_meta %>%
      left_join(
        lof_df, by = c("CHROM", "POS", "ALT", "Gene_ID", "Protein_ID")
      ) %>%
      left_join(
        nmd_df, by = c("CHROM", "POS", "ALT", "Gene_ID", "Protein_ID")
      )
    # Filter LOF & NMD cols
    split_ann_vcf_meta <- split_ann_vcf_meta %>%
      mutate(
        isLOF = case_when(
          !is.na(Total_Transcripts_LOF) &
            !is.na(HGVS.c) &
            # Filter out 3' UTR LOF hits!grepl("c\\.\\*", HGVS.c) &
            # Keep LOF baesd on Effect criteria
            grepl(paste(lof_nmd_effs, collapse = "|"), Effect) ~ T,
          .default = F
        ),
        isNMD = case_when(
          !is.na(Total_Transcripts_NMD) &
            !is.na(HGVS.c) &
            # Filter out 3' UTR NMD hits!grepl("c\\.\\*", HGVS.c) &
            # Keep NMD baesd on Effect criteria
            grepl(paste(lof_nmd_effs, collapse = "|"), Effect) ~ T,
          .default = F
        )
      ) %>%
      mutate(
        LOF = ifelse(isLOF, LOF, NA),
        Total_Transcripts_LOF = ifelse(isLOF, Total_Transcripts_LOF, NA),
        Percent_Transcripts_LOF = ifelse(isLOF, Percent_Transcripts_LOF, NA),
        NMD = ifelse(isNMD, NMD, NA),
        Total_Transcripts_NMD = ifelse(isNMD, Total_Transcripts_NMD, NA),
        Percent_Transcripts_NMD = ifelse(isNMD, Percent_Transcripts_NMD, NA)
      ) %>%
      # Guess types for newly split character columns with readr::parse_guess
      mutate(across(.cols = where(is.character), .fns = parse_guess))
    # Annotate variant types (e.g., SNP, INDEL)
    split_ann_vcf_meta <- split_ann_vcf_meta %>%
      mutate(
        Variant_Type = case_when(
          ALT == "*" ~ "INDEL",
          nchar(REF) != nchar(ALT) ~ "INDEL",
          .default = "SNP"
        ),
        .after = ALT
      )
    # Parse Error/Warning/Info into separate columns
    split_ann_vcf_meta <- split_ann_vcf_meta %>%
      mutate(Debug =
               case_when(!is.na(Debug) ~ str_split(Debug, "&") %>%
                           map(~ {
                             key_vals <- str_split_fixed(.x, "_", 2)
                             set_names(key_vals[, 2], key_vals[, 1])
                           }), .default = NA)) %>%
      unnest_longer(
        Debug,
        values_to = "Debug_Message",
        indices_to = "Debug_Code",
        keep_empty = T
      )
    debug_df <- split_ann_vcf_meta %>%
      select(CHROM,
             POS,
             REF,
             ALT,
             Effect,
             Protein_ID,
             Debug_Message,
             Debug_Code)
    # Save debugging messages
    vroom_write(debug_df, debug_df_file, delim = "\t")
    split_ann_vcf_meta <- split_ann_vcf_meta %>%
      select(!starts_with("Debug")) %>%
      distinct()
    
    # Save split table column types
    split_ann_col_types <- sapply(split_ann_vcf_meta, function(x) class(x)[1])
    split_ann_col_types <- data.frame(column = names(split_ann_col_types),
                                      type = unname(split_ann_col_types))
    write.table(
      split_ann_col_types,
      split_ann_col_types_file,
      sep = "\t",
      row.names = F,
      quote = F
    )
    # QC plots of all variant annotations
    # Save total number of variants x effects
    n_tot <- dim(split_ann_vcf_meta)[1]
    n_tot_lab <- paste0("n = ", scales::comma(n_tot))
    if (
      !file.exists(debug_plot) &
      !file.exists(var_qc_boxplot_plot) &
      !file.exists(eff_qc_boxplot_plot) &
      !file.exists(all_eff_plot) &
      !file.exists(all_impact_plot) &
      !file.exists(top10_eff_plot) &
      !file.exists(all_eff_impact_type_plot)
    ) {
      print("Plotting Effect QC boxplots...")
      # Effect vs AF
      p1 <- split_ann_vcf_meta %>%
        separate_longer_delim(Effect, delim = "&") %>%
        ggplot(aes(x = Effect, y = AF, fill = Effect)) +
        geom_boxplot(width = 0.5, outliers = F, show.legend = F) +
        scale_x_discrete(labels = as_labeller(pretty_labs)) +
        scale_fill_discrete(name = "Effect") +
        labs(x = "", y = "Allele Frequency") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      
      # Effect vs QUAL
      p2 <- split_ann_vcf_meta %>%
        separate_longer_delim(Effect, delim = "&") %>%
        ggplot(aes(x = Effect, y = QUAL, fill = Effect)) +
        geom_boxplot(width = 0.5, outliers = F, show.legend = F) +
        scale_x_discrete(labels = as_labeller(pretty_labs)) +
        scale_fill_discrete(name = "Effect") +
        labs(x = "", y = "Phred-Scaled Quality Score (QUAL)") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      
      # Effect vs Depth (DP)
      p3 <- split_ann_vcf_meta %>%
        separate_longer_delim(Effect, delim = "&") %>%
        ggplot(aes(x = Effect, y = DP, fill = Effect)) +
        geom_boxplot(width = 0.5, outliers = F, show.legend = F) +
        scale_x_discrete(labels = as_labeller(pretty_labs)) +
        scale_fill_discrete(name = "Effect") +
        labs(x = "", y = "Read Depth (DP)") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      
      # Effect vs QD
      p4 <- split_ann_vcf_meta %>%
        separate_longer_delim(Effect, delim = "&") %>%
        filter(!is.na(QD)) %>%
        ggplot(aes(x = Effect, y = QD, fill = Effect)) +
        geom_boxplot(width = 0.5, outliers = F, show.legend = F) +
        scale_x_discrete(labels = as_labeller(pretty_labs)) +
        scale_fill_discrete(name = "Effect") +
        labs(x = "", y = "Quality by Depth (QD)") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      
      # Effect vs AC
      p5 <- split_ann_vcf_meta %>%
        separate_longer_delim(Effect, delim = "&") %>%
        ggplot(aes(x = Effect, y = AC, fill = Effect)) +
        geom_boxplot(width = 0.5, outliers = F, show.legend = F) +
        scale_x_discrete(labels = as_labeller(pretty_labs)) +
        scale_fill_discrete(name = "Effect") +
        labs(x = "", y = "Allele Count in Genotypes (AC)") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      
      # Save Effect QC boxplot plots
      myboxplots <- ggarrange(p1, p2, p3, p4, p5, common.legend = T)
      myboxplots <- annotate_figure(
        myboxplots,
        top = paste0("QC by Effect (", n_tot_lab, ")"),
      )
      ggsave(
        plot = myboxplots,
        filename = eff_qc_boxplot_plot,
        bg = "white",
        width = 15,
        height = 15
      )
      rm(p1, p2, p3, p4, p5, myboxplots)
      
      print("Plotting Variant QC boxplots...")
      # Variant_Type vs AF
      p1 <- split_ann_vcf_meta %>%
        ggplot(aes(x = Variant_Type, y = AF, fill = Variant_Type)) +
        geom_boxplot(width = 0.5, outliers = F) +
        scale_fill_discrete(name = "Variant Type") +
        labs(x = "", y = "Allele Frequency") +
        theme_minimal()
      
      # Variant_Type vs QUAL
      p2 <- split_ann_vcf_meta %>%
        ggplot(aes(x = Variant_Type, y = QUAL, fill = Variant_Type)) +
        geom_boxplot(width = 0.5, outliers = F) +
        scale_fill_discrete(name = "Variant Type") +
        labs(x = "", y = "Phred-Scaled Quality Score (QUAL)") +
        theme_minimal()
      
      # Variant_Type vs Depth (DP)
      p3 <- split_ann_vcf_meta %>%
        ggplot(aes(x = Variant_Type, y = DP, fill = Variant_Type)) +
        geom_boxplot(width = 0.5, outliers = F) +
        scale_fill_discrete(name = "Variant Type") +
        labs(x = "", y = "Read Depth (DP)") +
        theme_minimal()
      
      # Variant_Type vs QD
      p4 <- split_ann_vcf_meta %>%
        filter(!is.na(QD)) %>%
        ggplot(aes(x = Variant_Type, y = QD, fill = Variant_Type)) +
        geom_boxplot(width = 0.5, outliers = F) +
        scale_fill_discrete(name = "Variant Type") +
        labs(x = "", y = "Quality by Depth (QD)") +
        theme_minimal()
      
      # Variant_Type vs AC
      p5 <- split_ann_vcf_meta %>%
        ggplot(aes(x = Variant_Type, y = AC, fill = Variant_Type)) +
        geom_boxplot(width = 0.5, outliers = F) +
        scale_fill_discrete(name = "Variant Type") +
        labs(x = "", y = "Allele Count in Genotypes (AC)") +
        theme_minimal()
      
      # Save variant QC boxplot plots
      myboxplots <- ggarrange(p1, p2, p3, p4, p5, common.legend = T)
      myboxplots <- annotate_figure(
        myboxplots,
        top = paste0("QC by Effect (", n_tot_lab, ")"),
      )
      ggsave(
        plot = myboxplots,
        filename = var_qc_boxplot_plot,
        bg = "white",
        width = 10,
        height = 6
      )
      rm(p1, p2, p3, p4, p5, myboxplots)
      
      print("Plotting effect bar graphs...")
      # Summary of all observed variant effects
      p <- split_ann_vcf_meta %>%
        separate_longer_delim(Effect, delim = "&") %>%
        count(Effect, sort = T) %>%
        ggplot(aes(y = reorder(Effect, n), x = n)) +
        geom_col(fill = "steelblue") +
        geom_text(aes(label = scales::comma(n)), hjust = 0, size = 3) +
        scale_x_log10(labels = scales::label_log(), limits = c(1, n_tot)) +
        scale_y_discrete(labels = as_labeller(pretty_labs)) +
        labs(subtitle = n_tot_lab, y = "Effect", x = "Count") +
        coord_cartesian(clip = "off") +
        theme_minimal()
      ggsave(
        plot = p,
        filename = all_eff_plot,
        bg = "white",
        width = 10,
        height = 6
      )
      rm(p)
      
      # Basic distribution of impacts
      p <- split_ann_vcf_meta %>%
        count(Impact) %>%
        ggplot(aes(x = Impact, y = n, fill = Impact)) +
        geom_col(show.legend = F) +
        geom_text(
          aes(
            label = paste0(
              scales::comma(n),
              " (",
              round((n / n_tot) * 100, digits = 2),
              "%)"
            )
          ),
          vjust = -1,
          size = 3
        ) +
        scale_y_continuous(labels = scales::label_comma()) +
        labs(subtitle = n_tot_lab, x = "Impact", y = "Count") +
        theme_minimal()
      ggsave(
        plot = p,
        filename = all_impact_plot,
        bg = "white",
        width = 7,
        height = 5
      )
      rm(p)
      
      # Variant_Type vs Impact/Effect
      p <- split_ann_vcf_meta %>%
        separate_longer_delim(Effect, delim = "&") %>%
        count(Variant_Type, Effect, Impact) %>%
        ggplot(aes(
          x = reorder(Effect, n),
          y = n,
          fill = Variant_Type
        )) +
        geom_col(position = position_dodge2(width = 0.9, preserve = "single")) +
        geom_text_repel(
          aes(label = ifelse(n < quantile(n, probs = 0.40), n, "")),
          # for dodged bars
          position = position_dodge2(width = 0.9, preserve = "single"),
          vjust = -0.3,
          size = 3
        ) +
        scale_x_discrete(labels = as_labeller(pretty_labs)) +
        scale_y_log10(labels = scales::label_log()) +
        scale_fill_discrete(name = "Type") +
        facet_grid(cols = vars(Impact),
                   scales = "free_x",
                   space = "free_x") +
        labs(
          title = "Variant Type and Effect by Impact",
          subtitle = n_tot_lab,
          x = "Effect",
          y = "Count"
        ) +
        theme_light() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      ggsave(
        plot = p,
        filename = all_eff_impact_type_plot,
        bg = "white",
        width = 10,
        height = 6
      )
      rm(p)
      
      # Plot of debugging messages
      print("Plotting debug messages by effect...")
      rm(split_ann_vcf_meta)
      debug_df <- vroom(debug_df_file,
                        delim = "\t",
                        col_types = cols(.default = col_character()))
      debug_df <- debug_df %>%
        mutate(across(.cols = where(is.character), .fns = parse_guess))
      p <- debug_df %>%
        filter(!if_any(starts_with("Debug_"), is.na)) %>%
        separate_longer_delim(Effect, delim = "&") %>%
        count(Debug_Message, Effect) %>%
        ggplot(aes(
          x = reorder(Effect, n),
          y = n,
          fill = Debug_Message
        )) +
        geom_col(position = position_dodge2(width = 0.9, preserve = "single")) +
        geom_text_repel(
          aes(label = ifelse(n < quantile(n, probs = 0.40), n, "")),
          # for dodged bars
          position = position_dodge2(width = 0.9, preserve = "single"),
          vjust = -0.3,
          size = 3
        ) +
        scale_x_discrete(labels = as_labeller(pretty_labs)) +
        scale_y_log10(labels = scales::label_log()) +
        scale_fill_discrete(labels = as_labeller(pretty_labs)) +
        labs(title = "Debug Messages by Variant Effect", x = "Effect") +
        theme_minimal() +
        theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top"
        )
      ggsave(
        plot = p,
        filename = debug_plot,
        bg = "white",
        width = 10,
        height = 6
      )
    }
    
    # Save split table
    vroom_write(split_ann_vcf_meta, split_ann_vcf_meta_file, delim = "\t")
  } else {
    print(paste("Reading existing formatted table:", split_ann_vcf_meta_file))
    # Load col type metadata
    split_ann_col_types <- read.delim(
      split_ann_col_types_file,
      stringsAsFactors = F
    )
    # Create a `col_types` specification for vroom
    vroom_col_types <- do.call(
      cols, setNames(lapply(split_ann_col_types$type, function(t)
        switch(
          t,
          character = col_character(),
          numeric   = col_double(),
          integer   = col_integer(),
          logical   = col_logical(),
          factor    = col_factor(),
          Date      = col_date(),
          col_guess()  # fallback
        )), split_ann_col_types$column)
    )
    split_ann_vcf_meta <- vroom(
      split_ann_vcf_meta_file,
      delim = "\t",
      col_types = vroom_col_types
    )
  }
  # Filter snpEff-annotated variants for genes of interest
  filt_ann_vcf_meta <- split_ann_vcf_meta %>%
    filter(Protein_ID %in% as.character(all_meiotic_prot_ids)) %>%
    left_join(all_meiotic_prot_annot, relationship = "many-to-many") %>%
    filter(POS >= Start & POS <= End)
  # Filter for high impact variants
  high_ann_vcf_meta <- filt_ann_vcf_meta %>%
    filter(Impact == "HIGH")
  # Write tab-delimited text file of CHROM POS to subset from full VCF
  subset_tab <- high_ann_vcf_meta %>%
    distinct(CHROM, POS)
  write.table(
    subset_tab,
    subset_tab_file,
    sep = "\t",
    quote = F,
    row.names = F,
    col.names = F
  )
}

# Parse subset of full annotated VCF (including genotypes)
# Verify subsetted VCF exists
if (!file.exists(sub_vcf_file)) {
  stop(paste("Subset", sub_vcf_file, "with subset_high_eff.sh before running."))
}

# # Ensure bgzip is accessible
# if (system("bgzip --version") != 0) {
#   stop("bgzip is not available in this environment.")
# }
# # Use pipe("bgzip ...") to decompress bgzipped VCF
# sub_vcf_raw <-
#   vroom(pipe(paste(
#     "bgzip -@ $SLURM_CPUS_PER_TASK -d -c", sub_vcf_file
#   )), delim = "\t", comment = "##")
sub_vcf_raw <- vroom(sub_vcf_file, delim = "\t", comment = "##")
sub_vcf <- sub_vcf_raw %>%
  # Remove remaining "#" in colnames
  rename_with(~ gsub("#", "", .)) %>%
  # Mark original row IDs
  mutate(orig_row = row_number(), .before = "CHROM") %>%
  # Drop irrelevant cols
  select(-c(ID, FILTER))
print(paste("Colnames:", paste(colnames(sub_vcf), collapse = " ")))
# Separate INFO col into columns named by regular expression (e.g., AC=##)
sub_vcf <- sub_vcf %>%
  mutate(INFO = str_split(INFO, ";") %>%
           map(~ {
             key_vals <- str_split_fixed(.x, "=", 2)
             set_names(key_vals[, 2], key_vals[, 1])
           })) %>%
  unnest_wider(INFO, names_repair = "unique")
print(
  paste("Colnames (INFO parsed):", paste(colnames(sub_vcf), collapse = " "))
)
# Strip special characters
sub_vcf <- sub_vcf %>%
  # Remove parentheses from LOF and NMD cols
  mutate(across(c(LOF, NMD), ~ gsub("\\(|\\)", "", .))) %>%
  # Keep only protein ID numbers from "jgi.p|Macpyr2|####" for downstream
  mutate(
    ANN = gsub("jgi\\.p\\|Macpyr2\\|", "", ANN),
    LOF = gsub("jgi\\.p_Macpyr2_", "", LOF),
    NMD = gsub("jgi\\.p_Macpyr2_", "", NMD)
  )
# Second reformat: Decompose VCF
sub_vcf <- sub_vcf %>%
  # Separate into rows when multiple ALT per variant
  separate_longer_delim(
    cols = all_of(contains_commas(.)),
    delim = ","
  ) %>%
  # Mark index of each allele split by ","
  group_by(orig_row) %>%
  mutate(decomp_index = row_number(), .after = orig_row) %>%
  ungroup()
# Third reformat: Split ANN col
split_sub_vcf <- sub_vcf %>%
  # Separate into rows when multiple ANN per variant
  separate_longer_delim(ANN, delim = ",") %>%
  # Split ANN into named columns
  separate_wider_delim(ANN, delim = "|", names = ann_names) %>%
  # Keep only matching ANN annotations
  filter(ALT == Allele)
# Fourth reformat: Split LOF & NMD cols
lof_df <- sub_vcf %>%
  filter(!is.na(LOF)) %>%
  select(CHROM, POS, ALT, LOF) %>%
  # Separate into additional rows when multiple per variant
  separate_longer_delim(LOF, delim = ",") %>%
  # Split LOF into named columns
  separate_wider_delim(LOF, delim = "|", names = lof_names)
nmd_df <- sub_vcf %>%
  filter(!is.na(NMD)) %>%
  select(CHROM, POS, ALT, NMD) %>%
  # Separate into additional rows when multiple per variant
  separate_longer_delim(NMD, delim = ",") %>%
  # Split NMD into named columns
  separate_wider_delim(NMD, delim = "|", names = nmd_names)
# Join LOF and NMD back to sub_vcf using Gene_ID and Protein_ID
split_sub_vcf <- split_sub_vcf %>%
  left_join(
    lof_df, by = c("CHROM", "POS", "ALT", "Gene_ID", "Protein_ID")
  ) %>%
  left_join(
    nmd_df, by = c("CHROM", "POS", "ALT", "Gene_ID", "Protein_ID")
  )
# Filter LOF & NMD cols
split_sub_vcf <- split_sub_vcf %>%
  mutate(
    isLOF = case_when(
      !is.na(Total_Transcripts_LOF) &
        !is.na(HGVS.c) &
        # Filter out 3' UTR LOF hits!grepl("c\\.\\*", HGVS.c) &
        # Keep LOF baesd on Effect criteria
        grepl(paste(lof_nmd_effs, collapse = "|"), Effect) ~ T,
      .default = F
    ),
    isNMD = case_when(
      !is.na(Total_Transcripts_NMD) &
        !is.na(HGVS.c) &
        # Filter out 3' UTR NMD hits!grepl("c\\.\\*", HGVS.c) &
        # Keep NMD baesd on Effect criteria
        grepl(paste(lof_nmd_effs, collapse = "|"), Effect) ~ T,
      .default = F
    )
  ) %>%
  mutate(
    LOF = ifelse(isLOF, LOF, NA),
    Total_Transcripts_LOF = ifelse(isLOF, Total_Transcripts_LOF, NA),
    Percent_Transcripts_LOF = ifelse(isLOF, Percent_Transcripts_LOF, NA),
    NMD = ifelse(isNMD, NMD, NA),
    Total_Transcripts_NMD = ifelse(isNMD, Total_Transcripts_NMD, NA),
    Percent_Transcripts_NMD = ifelse(isNMD, Percent_Transcripts_NMD, NA)
  ) %>%
  # Guess types for newly split character columns with readr::parse_guess
  mutate(across(.cols = where(is.character), .fns = parse_guess))
# Annotate variant types (e.g., SNP, INDEL)
split_sub_vcf <- split_sub_vcf %>%
  mutate(
    Variant_Type = case_when(
      ALT == "*" ~ "INDEL",
      nchar(REF) != nchar(ALT) ~ "INDEL",
      .default = "SNP"
    ),
    .after = ALT
  )
# Parse Error/Warning/Info into separate columns
split_sub_vcf <- split_sub_vcf %>%
  mutate(Debug =
           case_when(!is.na(Debug) ~ str_split(Debug, "&") %>%
                       map(~ {
                         key_vals <- str_split_fixed(.x, "_", 2)
                         set_names(key_vals[, 2], key_vals[, 1])
                       }), .default = NA)) %>%
  unnest_longer(
    Debug,
    values_to = "Debug_Message",
    indices_to = "Debug_Code",
    keep_empty = T
  )
debug_df <- split_sub_vcf %>%
  select(CHROM,
         POS,
         REF,
         ALT,
         Effect,
         Protein_ID,
         Debug_Message,
         Debug_Code)
# # Save debugging messages
# vroom_write(debug_df, debug_df_file, delim = "\t")
split_sub_vcf <- split_sub_vcf %>%
  select(!starts_with("Debug")) %>%
  distinct()
# # Save split table
# vroom_write(split_sub_vcf, split_sub_vcf_file, delim = "\t")


split_sub_vcf %>%
  # Decompose genotypes for each individual
  separate_wider_delim(
    cols = all_of(matches("\\d+")),
    delim = ":",
    names = format_names,
    names_sep = "_"
  )
  # Check for decomp_index == GT to mark genotypes with allele
  
  
