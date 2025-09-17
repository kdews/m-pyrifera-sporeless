# R Functions and Objects to Parse VCFs

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
    "DP_S",
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
  "Spo11",
  "Mre11",
  "Rad50",
  "Rad51",
  "RecA",
  "DMC1",
  "Red1",
  "Hop1",
  "Smc5",
  "MutS",
  "meio",
  "reproduct",
  "double-strand break",
  " homologous recomb"
)
# Meiosis-associated search terms ranked by importance
meiotic_rank <- c(
  "Spo11" = 1,
  "Mre11" = 1,
  "Rad50" = 1,
  "RecA" = 1,
  "DMC1" = 1,
  "Red1" = 1,
  "Hop1" = 1,
  "Smc5" = 1,
  "Exo1" = 1,
  "meio" = 2,
  "reproduct" = 2,
  "MutS" = 2,
  "double-strand break" = 3,
  " homologous recomb" = 3,
  "SCC" = 4,
  "recombination" = 4,
  "crossover" = 4,
  "telomere" = 4
)
# snpEff ANN severity scale
severity_rank <- c(
  "stop_gained" = 1,
  "frameshift_variant" = 2,
  "stop_lost" = 3,
  "start_lost" = 4,
  "splice_acceptor_variant" = 5,
  "splice_donor_variant" = 6,
  "exon_loss_variant" = 7,
  "splice_region_variant" = 8,
  "conservative_inframe_deletion" = 9,
  "conservative_inframe_insertion" = 10,
  "5_prime_UTR_variant" = 11,
  "3_prime_UTR_variant" = 12,
  "intron_variant" = 13
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
# Parse REF and ALT for AD and PL columns
parseRefAlt <- function(.data) {
  vars <- names(.data) %>%
    str_subset("_GT$") %>%
    sub("GT$", "", .)
  for (var in vars) {
    ad_col <- paste0(var, "AD")
    pl_col <- paste0(var, "PL")
    gt_col <- paste0(var, "GT")
    ad_ref_col <- paste0(ad_col, "_REF")
    ad_alt_col <- paste0(ad_col, "_ALT")
    pl_ref_col <- paste0(pl_col, "_REF")
    pl_alt_col <- paste0(pl_col, "_ALT")
    .data <- .data %>%
      rowwise() %>%
      mutate(
        !!ad_ref_col := if (is.na(.data[[ad_col]]) ||
                            is.na(.data[[gt_col]]))
          NA
        else
          as.integer(str_split_i(.data[[ad_col]], ",", 1)),
        !!ad_alt_col := if (is.na(.data[[gt_col]]) ||
                            is.na(.data[[ad_col]]) ||
                            .data[[gt_col]] == 0)
          NA
        else
          as.integer(
            str_split_i(
              .data[[ad_col]], ",", .data[[gt_col]] + 1
            )
          ),
        !!pl_ref_col := if (is.na(.data[[pl_col]]) ||
                            is.na(.data[[gt_col]]))
          NA
        else
          as.integer(str_split_i(.data[[pl_col]], ",", 1)),
        !!pl_alt_col := if (is.na(.data[[gt_col]]) ||
                            is.na(.data[[pl_col]]) ||
                            .data[[gt_col]] == 0)
          NA
        else
          as.integer(
            str_split_i(
              .data[[pl_col]], ",", .data[[gt_col]] + 1
            )
          )
      ) %>%
      ungroup()
  }
  .data
}
