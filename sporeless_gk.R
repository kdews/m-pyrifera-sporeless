# Clear environment
rm(list = ls())
# Required packages
library(tidyverse)
library(vroom)
library(ggpubr)
library(ggrepel)
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
  ann_vcf_file <-
    "raw_haploid_559_indv_on_CI_03_genome_final.info_only.ann.vcf.gz"
  # Gene feature file
  gff_file <- "genes.gff"
  # CSV of genes of interest
  gene_list_file <- "gene_list.csv"
  # Output directory
  outdir <- "m-pyrifera-sporeless/"
} else {
  line_args <- commandArgs(trailingOnly = T)
  meta_file <- line_args[1]
  vcf_id_file <- line_args[2]
  ann_vcf_file <- line_args[3]
  gff_file <- line_args[4]
  gene_list_file <- line_args[5]
  outdir <- line_args[6]
}
# Output
simple_meta_file <- "simple_metadata.tsv"
# Output table
ann_vcf_base <- gsub("\\..*", "", basename(ann_vcf_file))
ann_vcf_base <-
  tools::file_path_sans_ext(basename(ann_vcf_file), compression = T)
split_ann_vcf_file <- paste0(ann_vcf_base, ".split.tsv.gz")
debug_df_file <- paste0("debug", "_", ann_vcf_base, ".split.tsv.gz")
debug_plot <- "debug_snpeff.png"
var_qc_boxplot_plot <- "all_var_qc_boxplots.png"
eff_qc_boxplot_plot <- "all_eff_qc_boxplots.png"
all_eff_plot <- "all_eff_bar.png"
all_impact_plot <- "all_impact_bar.png"
top10_eff_plot <- "top10_eff_bar.png"
all_eff_impact_type_plot <- "all_eff_impact_type_bar.png"
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

# Parse annotated VCF
if (!file.exists(split_ann_vcf_file)) {
  # Ensure bgzip is accessible
  if (system("bgzip --version") != 0) {
    stop("bgzip is not available in this environment.")
  }
  # Use pipe("bgzip ...") to unzip large VCF file
  # ann_vcf_raw <- vroom(ann_vcf_file, delim = "\t", comment = "##")
  ann_vcf_raw <-
    vroom(pipe(paste(
      "bgzip -@ $SLURM_CPUS_PER_TASK -d -c", ann_vcf_file
    )), delim = "\t", comment = "##")
  # If interactive, keep random subset of annotated VCF + LOF/NMD lines
  # Otherwise, keep entire dataset in df
  if (interactive()) {
    lof_nmd_raw <-
      vroom(pipe(
        paste(
          "bgzip -@ $SLURM_CPUS_PER_TASK -d -c",
          "raw_haploid_559_indv_on_CI_03_genome_final.info_only.lof_nmd.ann.vcf.gz"
        )
      ), delim = "\t", comment = "##")
    set.seed(42)
    ann_vcf_raw <- ann_vcf_raw %>% slice_sample(n = 50000)
    lof_nmd_raw <- lof_nmd_raw %>% slice_sample(n = 50000)
    ann_vcf <- bind_rows(ann_vcf_raw, lof_nmd_raw)
  } else {
    ann_vcf <- ann_vcf_raw
  }
  # Remove remaining "#" in colnames
  ann_vcf <- ann_vcf %>%
    rename_with(~ gsub("#", "", .)) %>%
    # Drop irrelevant cols
    select(-c(ID, FILTER))
  print(paste("Colnames:", paste(colnames(ann_vcf), collapse = " ")))
  # Separate INFO col into columns named by regular expression (e.g., AC=##)
  ann_vcf <- ann_vcf %>%
    mutate(INFO = str_split(INFO, ";") %>%
             map(~ {
               key_vals <- str_split_fixed(.x, "=", 2)
               set_names(key_vals[, 2], key_vals[, 1])
             })) %>%
    unnest_wider(INFO, names_repair = "unique")
  print(paste("Colnames (INFO parsed):", paste(colnames(ann_vcf), collapse = " ")))
  # Strip special characters
  ann_vcf <- ann_vcf %>%
    # Remove parentheses from LOF and NMD cols
    mutate(across(c(LOF, NMD), ~ gsub("\\(|\\)", "", .))) %>%
    # Keep only protein ID numbers from "jgi.p|Macpyr2|####" for downstream
    mutate(
      ANN = gsub("jgi\\.p\\|Macpyr2\\|", "", ANN),
      LOF = gsub("jgi\\.p_Macpyr2_", "", LOF),
      NMD = gsub("jgi\\.p_Macpyr2_", "", NMD)
    )
  # Second reformat: Decompose VCF
  ann_vcf <- ann_vcf %>%
    # Separate into rows when multiple ALT per variant
    separate_longer_delim(cols = all_of(contains_commas(.)), delim = ",")
  # Third reformat: Split ANN col
  split_ann_vcf <- ann_vcf %>%
    # Separate into rows when multiple ANN per variant
    separate_longer_delim(ANN, delim = ",") %>%
    # Split ANN into named columns
    separate_wider_delim(ANN, delim = "|", names = ann_names) %>%
    # Keep only matching ANN annotations
    filter(ALT == Allele)
  # Fourth reformat: Split LOF & NMD cols
  lof_df <- ann_vcf %>%
    filter(!is.na(LOF)) %>%
    select(CHROM, POS, ALT, LOF) %>%
    # Separate into additional rows when multiple per variant
    separate_longer_delim(LOF, delim = ",") %>%
    # Split LOF into named columns
    separate_wider_delim(LOF, delim = "|", names = lof_names)
  nmd_df <- ann_vcf %>%
    filter(!is.na(NMD)) %>%
    select(CHROM, POS, ALT, NMD) %>%
    # Separate into additional rows when multiple per variant
    separate_longer_delim(NMD, delim = ",") %>%
    # Split NMD into named columns
    separate_wider_delim(NMD, delim = "|", names = nmd_names)
  # Join LOF and NMD back to ann_vcf using Gene_ID and Protein_ID
  split_ann_vcf <- split_ann_vcf %>%
    left_join(lof_df, by = c("CHROM", "POS", "ALT", "Gene_ID", "Protein_ID")) %>%
    left_join(nmd_df, by = c("CHROM", "POS", "ALT", "Gene_ID", "Protein_ID"))
  # Filter LOF & NMD cols
  split_ann_vcf <- split_ann_vcf %>%
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
  split_ann_vcf <- split_ann_vcf %>%
    mutate(
      Variant_Type = case_when(
        ALT == "*" ~ "INDEL",
        nchar(REF) != nchar(ALT) ~ "INDEL",
        .default = "SNP"
      ),
      .after = ALT
    )
  # Parse Error/Warning/Info into separate columns
  split_ann_vcf <- split_ann_vcf %>%
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
  debug_df <- split_ann_vcf %>%
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
  split_ann_vcf <- split_ann_vcf %>%
    select(!starts_with("Debug")) %>%
    distinct()
  # Save split table
  vroom_write(split_ann_vcf, split_ann_vcf_file, delim = "\t")
} else {
  print(paste("Reading existing formatted table:", split_ann_vcf_file))
  split_ann_vcf <- vroom(
    split_ann_vcf_file,
    delim = "\t",
    col_types = cols(.default = col_character())
  )
  print("Guessing column types...")
  split_ann_vcf <- split_ann_vcf %>%
    mutate(across(.cols = where(is.character), .fns = parse_guess))
}

# Plots
# Save total number of variants x effects
n_tot <- dim(split_ann_vcf)[1]
n_tot_lab <- paste0("n = ", scales::comma(n_tot))
plots <- F
if (plots == T){
  print("Plotting Effect QC boxplots...")
  # Effect vs AF
  p1 <- split_ann_vcf %>%
    separate_longer_delim(Effect, delim = "&") %>%
    ggplot(aes(x = Effect, y = AF, fill = Effect)) +
    geom_boxplot(width = 0.5, outliers = F, show.legend = F) +
    scale_x_discrete(labels = as_labeller(pretty_labs)) +
    scale_fill_discrete(name = "Effect") +
    labs(x = "", y = "Allele Frequency") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Effect vs QUAL
  p2 <- split_ann_vcf %>%
    separate_longer_delim(Effect, delim = "&") %>%
    ggplot(aes(x = Effect, y = QUAL, fill = Effect)) +
    geom_boxplot(width = 0.5, outliers = F, show.legend = F) +
    scale_x_discrete(labels = as_labeller(pretty_labs)) +
    scale_fill_discrete(name = "Effect") +
    labs(x = "", y = "Phred-Scaled Quality Score (QUAL)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Effect vs Depth (DP)
  p3 <- split_ann_vcf %>%
    separate_longer_delim(Effect, delim = "&") %>%
    ggplot(aes(x = Effect, y = DP, fill = Effect)) +
    geom_boxplot(width = 0.5, outliers = F, show.legend = F) +
    scale_x_discrete(labels = as_labeller(pretty_labs)) +
    scale_fill_discrete(name = "Effect") +
    labs(x = "", y = "Read Depth (DP)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Effect vs QD
  p4 <- split_ann_vcf %>%
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
  p5 <- split_ann_vcf %>%
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
  p1 <- split_ann_vcf %>%
    ggplot(aes(x = Variant_Type, y = AF, fill = Variant_Type)) +
    geom_boxplot(width = 0.5, outliers = F) +
    scale_fill_discrete(name = "Variant Type") +
    labs(x = "", y = "Allele Frequency") +
    theme_minimal()
  
  # Variant_Type vs QUAL
  p2 <- split_ann_vcf %>%
    ggplot(aes(x = Variant_Type, y = QUAL, fill = Variant_Type)) +
    geom_boxplot(width = 0.5, outliers = F) +
    scale_fill_discrete(name = "Variant Type") +
    labs(x = "", y = "Phred-Scaled Quality Score (QUAL)") +
    theme_minimal()
  
  # Variant_Type vs Depth (DP)
  p3 <- split_ann_vcf %>%
    ggplot(aes(x = Variant_Type, y = DP, fill = Variant_Type)) +
    geom_boxplot(width = 0.5, outliers = F) +
    scale_fill_discrete(name = "Variant Type") +
    labs(x = "", y = "Read Depth (DP)") +
    theme_minimal()
  
  # Variant_Type vs QD
  p4 <- split_ann_vcf %>%
    filter(!is.na(QD)) %>%
    ggplot(aes(x = Variant_Type, y = QD, fill = Variant_Type)) +
    geom_boxplot(width = 0.5, outliers = F) +
    scale_fill_discrete(name = "Variant Type") +
    labs(x = "", y = "Quality by Depth (QD)") +
    theme_minimal()
  
  # Variant_Type vs AC
  p5 <- split_ann_vcf %>%
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
  p <- split_ann_vcf %>%
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
  p <- split_ann_vcf %>%
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
  p <- split_ann_vcf %>%
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
  rm(split_ann_vcf)
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

# Select genes of interest
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
print(paste(
  "GFF colnames (Attributes parsed):",
  paste(colnames(gff_split), collapse = " ")
))

gene_list <- read_csv(gene_list_file)
gene_list %>%
  mutate(Strand = gsub("\\(|\\)", "", Strand)) %>%
  select(
    Model_Name = Name,
    Protein_ID = `Protein Id`,
    Transcript_ID = `Transcript Id`,
    CHROM = Scaffold,
    Start,
    End,
    Strand
  )

head(gff)
head(gff_split)
