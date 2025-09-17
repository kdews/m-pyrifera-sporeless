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
  # Table of sample metadata
  sample_meta_file <- "070721_metadata.csv"
  # Annotated VCF (metadata only)
  ann_vcf_meta_file <-
    "raw_haploid_559_indv_on_CI_03_genome_final.info_only.ann.vcf.gz"
  # Gene feature file
  gff_file <- "genes.gff"
  # CSV of meiotic genes from JGI GUI
  meiotic_gene_list_file <- "jgi_gui_meiotic_genes.csv"
  # Prefix for gene annotation table filenames
  gene_annot_tab_base <- "Macpyr2_GeneCatalog_proteins_20220914"
  # Output directory
  outdir <- "m-pyrifera-sporeless/"
} else {
  line_args <- commandArgs(trailingOnly = T)
  sample_meta_file <- line_args[1]
  ann_vcf_meta_file <- line_args[2]
  gff_file <- line_args[3]
  meiotic_gene_list_file <- line_args[4]
  gene_annot_tab_base <- line_args[5]
  outdir <- line_args[6]
}
# VCF ID file
vcf_id_file <- "raw_vcf_ids.txt"
# Gene annotation tables
# Use prefix to get filenames
gene_annot_tab_files <- list.files(
  pattern = paste0(gene_annot_tab_base, ".+\\.tab\\.gz")
)
go_tab_file <- grep("_GO\\.", gene_annot_tab_files, value = T)
ipr_tab_file <- grep("_IPR\\.", gene_annot_tab_files, value = T)
kegg_tab_file <- grep("_KEGG\\.", gene_annot_tab_files, value = T)
kog_tab_file <- grep("_KOG\\.", gene_annot_tab_files, value = T)

# Output
# Sample ID index
sample_idx_file <- "sample_idx.tsv"
# Parsed VCF metadata files
ann_vcf_meta_base <-
  tools::file_path_sans_ext(basename(ann_vcf_meta_file), compression = T)
split_ann_vcf_meta_file <- paste0(ann_vcf_meta_base, ".split.tsv.gz")
split_ann_col_types_file <- paste0(ann_vcf_meta_base, ".split.col_types.tsv")
# Table of snpEff debugging messages
debug_df_file <- paste0("debug", "_", ann_vcf_meta_base, ".split.tsv.gz")
# Plot filenames
var_qc_boxplot_plot <- "all_var_qc_boxplots.png"
eff_qc_boxplot_plot <- "all_eff_qc_boxplots.png"
all_eff_plot <- "all_eff_bar.png"
all_impact_plot <- "all_impact_bar.png"
all_eff_impact_type_plot <- "all_eff_impact_type_bar.png"
debug_plot <- "debug_snpeff.png"
# Meiosis-associated protein annotations
all_meiotic_prot_ids_file <- "all_meiotic_protein_IDs.txt"
all_meiotic_prot_annot_file <- "all_meiotic_proteins.tab"
# 2-column tsv to subset variants of interest from VCF
subset_tab_file <- "high_eff_subset.tsv"
# If exists, prepend output directory to output filenames of plots
if (dir.exists(outdir)) {
  debug_plot <- paste0(outdir, debug_plot)
  var_qc_boxplot_plot <- paste0(outdir, var_qc_boxplot_plot)
  eff_qc_boxplot_plot <- paste0(outdir, eff_qc_boxplot_plot)
  all_eff_plot <- paste0(outdir, all_eff_plot)
  all_impact_plot <- paste0(outdir, all_impact_plot)
  all_eff_impact_type_plot <- paste0(outdir, all_eff_impact_type_plot)
}

# Load shared functions and objects
source("m-pyrifera-sporeless/vcf_parsing_functions.R")

# Data analysis
# Table correlating sample IDs and gametophyte codes
if (!file.exists(sample_idx_file)) {
  # Read sample metadata
  sample_meta <- read_csv(sample_meta_file)
  vcf_ids <- read_table(vcf_id_file, col_names = "vcfSampleID")
  vcf_ids <- vcf_ids %>%
    arrange(vcfSampleID) %>%
    pull(vcfSampleID)
  sample_idx <- sample_meta %>%
    select(SampleID, GametophyteCode) %>%
    unique() %>%
    arrange(SampleID)
  # Check ID index versus VCF IDs
  all(sample_idx$SampleID == vcf_ids)
  sample_idx <- sample_idx %>%
    mutate(Sex = str_extract(GametophyteCode, "(?<=\\.)[MF](?=\\.)"))
  write_tsv(sample_idx, sample_idx_file)
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
  # Read meiotic protein IDs from JGI GUI
  meiotic_gene_list <- read_csv(meiotic_gene_list_file)
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
    summarize(IPR=paste(na.omit(iprDesc), collapse = ";"), .by = proteinId)
  kegg_tab <- read_tsv(kegg_tab_file, na = c("\\N")) %>%
    # Remove "#" from colnames
    rename_with(~ gsub("#", "", .)) %>%
    # Collapse to one line per proteinId
    summarize(KEGG=paste(na.omit(definition), collapse = ";"), .by = proteinId)
  kog_tab <- read_tsv(kog_tab_file, na = c("\\N")) %>%
    # Remove "#" from colnames
    rename_with(~ gsub("#", "", .)) %>%
    # Collapse to one line per proteinId
    summarize(KOG=paste(na.omit(kogdefline), collapse = ";"), .by = proteinId)
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
    meiotic_gene_list %>%
      select(proteinId = `Protein Id`)
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
    filter(Type != "mRNA")
  # Add any user annotations
  all_meiotic_prot_annot <- all_meiotic_prot_annot %>%
    mutate(
      Protein_Product = case_when(
        as.numeric(Protein_ID) == 9336431 ~
          paste(Protein_Product, "exonuclease I (Exo1)", sep = ";"),
        .default = Protein_Product
      )
    )
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

# Write 2-column (CHR POS) table to subset VCF
# Keep only high impact variants in meiosis genes
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
  all_meiotic_prot_annot <- all_meiotic_prot_annot %>%
    mutate(Protein_ID = as.character(Protein_ID))
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
