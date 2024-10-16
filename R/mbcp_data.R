#' MBCproject clinical data
#'
#' Data frame with clinical data for each RNA or Wes sample in the MBCproject
"mbcp_clinical_data"

#' MBCproject log2 TPMs
#'
#' Numeric matrix with log2 TPMs for each gene in the MBCproject.
#' Rows are genes identified by hugo symbol and columns are sample aliases
#'
"mbcp_log2_tpms"

#' MBCproject mutations
#'
#' Data frame in MAF format with mutations found in MBCproject.
#' It includes mandatory fields in a MAF file, plus columns added by ABSOLUTE
#' and ONCOKB annotator
"mbcp_mutations"

#' MBCproject copy number variations
#' Data frame in gene level seg file format, plus columns added by ABSOLUTE and
#' ONCOKB annotator
"mbcp_cnvs"

#' MBCproject enrichment scores
#'
#' Numeric matrix with normalized enrichment scores for hallmark and selected sets in the MBCproject
#' Rows are gene sets and columns are sample aliases. Scores were obtained using the
#' GSEA implementation in the fgsea package
#'
"mbcp_enrichment_scores"

#' MBCproject enrichment classification
#'
#' Classification of normalized enrichment scores for hallmark signatures and other selected sets
#' obtained using a method for cohort balancing according to receptor status
#' Categories: Lower Decile, Lower Quartile, Upper Quartile, Upper Decile.
#' `NA` values for entries that do not fall into these categories
"mbcp_enrichment_class"

#' MBCproject log2 tpms classification
#'
#' Classification of log2 tpms values of gene expression
#' obtained using a method for cohort balancing according to receptor status
#' Categories: Lower Decile, Lower Quartile, Upper Quartile, Upper Decile.
#' `NA` values for entries that do not fall into these categories
"mbcp_log2_tpms_class"

#' MBCproject germline hits
#'
#' Must remove this for public release
"mbcp_germline_hits"
