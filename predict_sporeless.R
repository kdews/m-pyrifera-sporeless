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
  # High impact variants (subset of annotated VCF)
  sub_vcf_file <-
    "high_eff_subset_raw_haploid_559_indv_on_CI_03_genome_final.ann.vcf.gz"
  # Output directory
  outdir <- "m-pyrifera-sporeless/"
} else {
  line_args <- commandArgs(trailingOnly = T)
  sub_vcf_file <- line_args[1]
  outdir <- line_args[2]
}
# Sample ID index
sample_idx_file <- "sample_idx.tsv"
# Meiosis-associated protein annotations
all_meiotic_prot_ids_file <- "all_meiotic_protein_IDs.txt"
all_meiotic_prot_annot_file <- "all_meiotic_proteins.tab"

# Output
# Parsed VCF table filenames
sub_vcf_base <-
  tools::file_path_sans_ext(basename(sub_vcf_file), compression = T)
split_sub_vcf_file <- paste0(sub_vcf_base, ".split.tsv.gz")
split_sub_col_types_file <- paste0(sub_vcf_base, ".split.col_types.tsv")
prot_split_sub_vcf_file <- paste0(sub_vcf_base, ".prot.split.tsv.gz")
prot_split_sub_col_types_file <- paste0(
  sub_vcf_base, ".prot.split.col_types.tsv"
)
# Ranked crosses
cross_table_file <- "cross_table.tsv"
cross_list_file <- "cross_list.tsv"

# Load shared functions and objects
source("m-pyrifera-sporeless/vcf_parsing_functions.R")

# Data analysis
# Read sample ID index
if (file.exists(sample_idx_file)) {
  sample_idx <- read_tsv(sample_idx_file)
}
n_indivs <- dim(sample_idx)[1]

# Read table of meiosis-associated protein annotations
if (file.exists(all_meiotic_prot_annot_file)) {
  all_meiotic_prot_annot <- read_tsv(all_meiotic_prot_annot_file)
}

# Parse subset of full annotated VCF (including genotypes)
if (!file.exists(split_sub_vcf_file)) {
  # Verify subset VCF exists
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
    mutate(ROW_ID = row_number(), .before = "CHROM") %>%
    # Drop irrelevant cols
    select(-c(ID, FILTER, FORMAT))
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
    group_by(ROW_ID) %>%
    mutate(ALT_IDX = row_number(), .after = ALT) %>%
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
  # Keep LOF & NMD col info only on matching rows
  split_sub_vcf <- split_sub_vcf %>%
    mutate(
      isLOF = case_when(
        !is.na(Total_Transcripts_LOF) &
          !is.na(HGVS.c) &
          # # Filter out 3' UTR LOF hits
          # !grepl("c\\.\\*", HGVS.c) &
          # Keep LOF based on Effect criteria
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
      ),
      LOF = ifelse(isLOF, LOF, NA),
      Total_Transcripts_LOF = ifelse(isLOF, Total_Transcripts_LOF, NA),
      Percent_Transcripts_LOF = ifelse(isLOF, Percent_Transcripts_LOF, NA),
      NMD = ifelse(isNMD, NMD, NA),
      Total_Transcripts_NMD = ifelse(isNMD, Total_Transcripts_NMD, NA),
      Percent_Transcripts_NMD = ifelse(isNMD, Percent_Transcripts_NMD, NA)
    )
  # Guess types for newly split character columns with readr::parse_guess
  split_sub_vcf <- split_sub_vcf %>%
    mutate(across(.cols = where(is.character), .fns = parse_guess))
  # Annotate variant types (e.g., SNP, INDEL)
  split_sub_vcf <- split_sub_vcf %>%
    mutate(
      Variant_Type = case_when(
        ALT == "*" ~ "INDEL",
        nchar(REF) != nchar(ALT) ~ "INDEL",
        .default = "SNP"
      ),
      .after = ALT_IDX
    )
  split_sub_vcf <- split_sub_vcf %>%
    select(!starts_with("Debug")) %>%
    distinct()
  # Parse genotype columns for each sample
  split_sub_vcf <- split_sub_vcf %>%
    # Separate "FORMAT" specified IDs
    separate_wider_delim(
      cols = all_of(matches("^\\d+$")),
      delim = ":",
      names = format_names,
      names_sep = "_"
    ) %>%
    # Sub missing (.) with NA
    mutate(
      across(ends_with(c("_GT", "_PL", "_DP_S", "_GQ")), ~ gsub("\\.", NA, .))
    ) %>%
    # Coerce genotype (GT), sample read depth (DP_S), and genotype quality (GQ)
    # to integers
    mutate(across(ends_with(c("_GT", "_DP_S", "_GQ")), as.integer)) %>%
    # Guess types for newly split character columns with readr::parse_guess
    mutate(across(where(is.character), parse_guess)) %>%
    rowwise() %>%
    mutate(
      # Check for ALT_IDX == GT to sum number of genotypes with given allele
      N_ALT_GTS = sum(c_across(ends_with("_GT")) == ALT_IDX, na.rm = T),
      N_CALLS = sum(!is.na(c_across(ends_with("_GT")))),
      .after = ALT_IDX
    ) %>%
    ungroup() %>%
    # Mark number of possible ALT genotypes (N_ALT_ALLELES) per site
    mutate(
      N_ALT_ALLELES = n(),
      `N_ALT_GTS:REF>ALT` = paste(
        N_ALT_GTS,
        paste(REF, ALT, sep = ">"),
        sep = ":",
        collapse = ", "
      ),
      .by = c(CHROM, POS)
    )
  # Parse columns with REF and ALT, matching GT to ALT_IDX
  split_sub_vcf <- split_sub_vcf %>%
    parseRefAlt()
  # Filter for only HIGH impact variants following VCF decomposition
  split_sub_vcf <- split_sub_vcf %>%
    filter(Impact == "HIGH")
  # Save split table column types
  split_sub_col_types <- sapply(split_sub_vcf, function(x) class(x)[1])
  split_sub_col_types <- data.frame(column = names(split_sub_col_types),
                                    type = unname(split_sub_col_types))
  write.table(
    split_sub_col_types,
    split_sub_col_types_file,
    sep = "\t",
    row.names = F,
    quote = F
  )
  # Save parsed table
  vroom_write(split_sub_vcf, split_sub_vcf_file, delim = "\t")
} else {
  # Load col type metadata
  split_sub_col_types <- read.delim(
    split_sub_col_types_file,
    stringsAsFactors = F
  )
  # Create a `col_types` specification for vroom
  vroom_col_types <- do.call(
    cols, setNames(lapply(split_sub_col_types$type, function(t)
      switch(
        t,
        character = col_character(),
        numeric   = col_double(),
        integer   = col_integer(),
        logical   = col_logical(),
        factor    = col_factor(),
        Date      = col_date(),
        col_guess()  # fallback
      )), split_sub_col_types$column)
  )
  split_sub_vcf <- vroom(
    split_sub_vcf_file,
    delim = "\t",
    col_types = vroom_col_types
  )
}

# Add protein annotation to table and rank by variant severity and relevance
if (!file.exists(prot_split_sub_vcf_file)) {
  # Annotate split subset VCF lines with protein information
  prot_split_sub_vcf <- split_sub_vcf %>%
    mutate(Protein_ID = as.numeric(Protein_ID)) %>%
    left_join(
      filter(all_meiotic_prot_annot, Type == "gene")
    ) %>%
    relocate(
      all_of(matches(setdiff(
        colnames(all_meiotic_prot_annot),
        colnames(split_sub_vcf))
      )),
      .after = Gene_ID
    )
  
  # Sort variants for analysis
  prot_split_sub_vcf <- prot_split_sub_vcf %>%
    rowwise() %>%
    # Rank variant effects by severity
    mutate(
      Effect_Rank = str_split(Effect, "&"),
      Effect_Rank = min(severity_rank[Effect_Rank], na.rm = T),
      .before = 1
    ) %>%
    # Collapse all protein annotation columns into one string
    mutate(Annotations = paste(na.omit(unique(unlist(
      str_split(c_across(c(
        IPR, KEGG, GO, KOG, Protein_Product
      )), ";")
    ))), collapse = ";"),
    .after = Protein_Product) %>%
    # Rank protein annotations by significance to meiosis
    mutate(
      # Find which meiotic search terms match annotations in each row
      Meiotic_Terms = list(
        sapply(
          names(meiotic_rank), function(term) {
            str_detect(Annotations, fixed(term, ignore_case = T))
          }
        )
      ),
      Meiotic_Terms = list(names(meiotic_rank)[unlist(Meiotic_Terms)]),
      Annot_Rank = list(meiotic_rank[Meiotic_Terms]),
      # Assign min rank among matching meiotic terms (no match = lowest rank)
      Annot_Rank = if (length(Annot_Rank) > 0)
        min(Annot_Rank)
      else
        max(meiotic_rank) + 1,
      Meiotic_Terms = paste(Meiotic_Terms, collapse = ", "),
      Meiotic_Terms = case_when(Meiotic_Terms == "" ~ NA,
                                .default = Meiotic_Terms),
      .after = Effect_Rank
    ) %>%
    ungroup() %>%
    # Filter out singletons
    filter(N_ALT_GTS > 1)
  # Rank by effect severity, then annotation relevance, then allele rarity
  prot_split_sub_vcf <- prot_split_sub_vcf %>%
    mutate(
      Comb_Rank = Effect_Rank + Annot_Rank,
      .before = Effect_Rank
    ) %>%
    arrange(
      Comb_Rank,
      AF,
      N_ALT_ALLELES
    )
  # Save split table column types
  prot_split_sub_col_types <- sapply(prot_split_sub_vcf, function(x) class(x)[1])
  prot_split_sub_col_types <- data.frame(column = names(prot_split_sub_col_types),
                                    type = unname(prot_split_sub_col_types))
  write.table(
    prot_split_sub_col_types,
    prot_split_sub_col_types_file,
    sep = "\t",
    row.names = F,
    quote = F
  )
  # Save parsed table
  vroom_write(prot_split_sub_vcf, prot_split_sub_vcf_file, delim = "\t")
} else {
  # Load col type metadata
  prot_split_sub_col_types <- read.delim(
    prot_split_sub_col_types_file,
    stringsAsFactors = F
  )
  # Create a `col_types` specification for vroom
  vroom_col_types <- do.call(
    cols, setNames(lapply(prot_split_sub_col_types$type, function(t)
      switch(
        t,
        character = col_character(),
        numeric   = col_double(),
        integer   = col_integer(),
        logical   = col_logical(),
        factor    = col_factor(),
        Date      = col_date(),
        col_guess()  # fallback
      )), prot_split_sub_col_types$column)
  )
  prot_split_sub_vcf <- vroom(
    prot_split_sub_vcf_file,
    delim = "\t",
    col_types = vroom_col_types
  )
}

# Reframe data to rank crosses by variant rank
if (!file.exists(cross_table_file) & !file.exists(cross_list_file)){
  # Pivot variant table to per-sample orientation
  per_sample_variants <- prot_split_sub_vcf %>%
    # Drop old stat columns not split on REF/ALT
    select(-ends_with(c("_PL", "_AD"))) %>%
    pivot_longer(
      # Select columns that start with a number, underscore, then a field name
      cols = matches("^\\d+_"),
      # ".value" means keep field name as column
      names_to = c("SampleID", ".value"),
      names_pattern = "^(\\d+)_(.*)$",
      names_transform = list(SampleID = as.integer)
    ) %>%
    relocate(SampleID:last_col(), .after = Annot_Rank) %>%
    # Merge with sample ID index
    left_join(sample_idx, relationship = "many-to-one") %>%
    relocate(matches(colnames(sample_idx)), .after = SampleID)
  
  # Quality-filter samples
  filt_samp <- per_sample_variants %>%
    filter(
      # Keep only ALT variants
      GT == ALT_IDX,
      # Filter for rare minor alleles
      AF <= 0.15,
      # Filter out variants with high missingness
      N_CALLS/n_indivs >= 0.8,
      # Minimum coverage
      DP_S >= 10,
      # Minimum genotype quality
      GQ >= 20
    ) %>%
    # Keep only variants where at least one female and male pair exist
    filter(length(unique(Sex)) == 2, .by = c(CHROM, POS, REF, ALT)) %>%
    # Keep top 2 samples of each sex for each variant
    group_by(CHROM, POS, REF, ALT, Sex) %>%
    arrange(
      Comb_Rank,
      N_ALT_GTS,
      AF,
      desc(GQ),
      desc(DP_S),
      desc(AD_ALT),
      .by_group = T
    ) %>%
    slice_head(n = 2) %>%
    ungroup() %>%
    arrange(
      Comb_Rank,
      N_ALT_GTS,
      AF,
      desc(GQ),
      desc(DP_S),
      desc(AD_ALT)
    )
  
  # Generate ranked crosses from gametophyte codes
  cross_df <- filt_samp %>%
    reframe(
      expand.grid(
        F_ID = GametophyteCode[Sex == "F"],
        M_ID = GametophyteCode[Sex == "M"]
      ),
      .by = -c(SampleID:PL_ALT)
    ) %>%
    mutate(Cross = paste(F_ID, M_ID, sep = " x "),
           across(c(F_ID, M_ID, Cross), as.character)) %>%
    relocate(F_ID:last_col(), .before = 1)
  # Create simple list of ranked crosses
  cross_list <- cross_df %>%
    distinct(Cross, F_ID, M_ID) %>%
    rowid_to_column("Rank")
  
  # Save detailed dataframe of crosses
  write_tsv(cross_df, cross_table_file)
  # Save basic list of crosses
  write_tsv(cross_list, cross_list_file)
}


