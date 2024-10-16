utils::globalVariables(c("mbcp_mutations", "mbcp_germline_hits", "mbcp_cnvs"))

#' Get mutations in MBC project
#'
#' @param samples Character vector of samples to filter variations by sample
#' @param genes Character vector of hugo symbols to filter variations by gene
#' @param min_ccf Minimum cancer cell fraction
#' @return A data frame with mutations in a MAF format
#' @export
#'
#' @examples
#' get_mbcp_mutations()
get_mbcp_mutations <- function(samples = NULL, genes = NULL, min_ccf = NULL) {

  stopifnot("samples should be a character vector" =  is.null(samples) | is.character(samples),
            "genes should be a character vector" = is.null(genes) | is.character(genes),
            "min_ccf should be a numeric value" =  is.null(min_ccf) | is.numeric(min_ccf) & length(min_ccf) == 1)

  all_samples <- unique(mbcp_mutations$Sample_ID)

  if(any(!samples %in% all_samples)) {
    warning("some samples were not found in the mbcp_mutations dataset")
  }

  mutations <- mbcp_mutations |>
    dplyr::filter(is.null(samples) | .data$Sample_ID %in% samples,
                  is.null(genes) | .data$Hugo_Symbol %in% genes)

  if(!is.null(min_ccf)) {
    mutations <- mutations |>
      dplyr::filter(.data$ccf_hat > min_ccf)
  }

  return(mutations)
}

#' Get CNV in MBC project
#' @inheritParams get_mbcp_mutations
#' @export
#'
#' @examples
#' get_mbcp_cnvs()
get_mbcp_cnvs <- function(samples = NULL, genes = NULL, min_ccf = NULL) {

  stopifnot("samples should be a character vector" =  is.null(samples) | is.character(samples),
            "genes should be a character vector" = is.null(genes) | is.character(genes),
            "min_ccf should be a numeric value" =  is.null(min_ccf) | is.numeric(min_ccf) & length(min_ccf) == 1)

  all_samples <- unique(mbcp_cnvs$sample_id)

  if(any(!samples %in% all_samples)) {
    warning("some samples were not found in the mbcp_cnvs dataset")
  }

  cnvs <- mbcp_cnvs |>
    dplyr::filter(is.null(samples) | .data$sample_id %in% samples,
                  is.null(genes) | .data$hugo_symbol %in% genes)

  if(!is.null(min_ccf)) {
    cnvs <- cnvs |>
      dplyr::filter(.data$cancer_cell_frac_a2 > min_ccf | .data$cancer_cell_frac_a1 > min_ccf)
  }

  return(cnvs)

}

#' Get dna alterations (CNVs and mutations) in the MBC project.
#' @inheritParams get_mbcp_mutations
#' @export
#'
#' @examples
#' get_mbcp_dna_alts_summary(min_ccf = 0.75)
get_mbcp_dna_alts_summary <- function(samples = NULL, genes = NULL, min_ccf = NULL) {

  stopifnot("samples should be a character vector" =  is.null(samples) | is.character(samples),
            "genes should be a character vector" = is.null(genes) | is.character(genes),
            "min_ccf should be a numeric value" =  is.null(min_ccf) | is.numeric(min_ccf) & length(min_ccf) == 1)

  amps_categories <- get_mbcp_labels("dna_alts")[3:6]

  mutations <- get_mbcp_mutations(samples, genes, min_ccf) |>
     dplyr::select(hugo_symbol = Hugo_Symbol, sample_id = Sample_ID,
                   protein_change = Protein_Change,  ccf_hat,
                   variant_type = Variant_Type, variant_classification = Variant_Classification) |>
     dplyr::group_by(hugo_symbol, sample_id) |>
     dplyr::summarise(variant_classification = paste0(variant_classification, collapse = ","),
                      n_alts = dplyr::n(),
                      protein_change = paste0(protein_change, collapse = ","),
                      ccfs = paste0(ccf_hat, collapse = ",")) |>
     dplyr::mutate(final_variant_classification = dplyr::if_else(n_alts == 1, variant_classification, "Multi_hit"),
                   variant_type = "MUT") |>
     dplyr::ungroup()

  ### Remove GAINs
  cnvs <- get_mbcp_cnvs(samples, genes, min_ccf) |>
    dplyr::select(.data$hugo_symbol, .data$sample_id, .data$cnap,  .data$variant_type,
                   .data$variant_classification, .data$cancer_cell_frac_a2, .data$cancer_cell_frac_a1) |>
    dplyr::mutate(variant_classification = gsub("GAIN,", "", variant_classification),
                   ccf_hat = pmax(.data$cancer_cell_frac_a2, .data$cancer_cell_frac_a1, na.rm = T)) |>
    dplyr::filter(variant_classification != "GAIN") |>
    tidyr::separate_wider_delim(.data$variant_classification, ",", names_sep = "", too_few = "align_start") |>
    dplyr::select(-.data$cancer_cell_frac_a1, -.data$cancer_cell_frac_a2)

  ### Pick the highest amplification
   amps <- cnvs |>
    dplyr::filter(.data$variant_classification1 %in% amps_categories)

   amps_count <-  amps |>
    dplyr::count(.data$hugo_symbol, .data$sample_id)

   all_amps <- amps |>
      dplyr::anti_join(amps_count |>
                        dplyr::filter(n == 1),
                      by = c("sample_id", "hugo_symbol")) |>
      dplyr::mutate(factor_vc = factor(variant_classification1, levels = amps_categories, ordered = TRUE)) |>
      dplyr::group_by(.data$hugo_symbol, .data$sample_id) |>
      dplyr::summarise(n_alts = dplyr::n(),
                      variant_classification1 = paste0(variant_classification1, collapse = ","),
                      variant_classification2 = unique(variant_classification2),
                      ccfs = paste0(ccf_hat, collapse = ","),
                      cnaps = paste0(cnap, collapse = ","),
                      final_variant_classification = max(factor_vc)) |>
     dplyr::ungroup() |>
     dplyr::bind_rows(
       amps |>
         dplyr::semi_join(amps_count |>
                            dplyr::filter(n == 1),
                          by = c("sample_id", "hugo_symbol")) |>
         dplyr::mutate(n_alts = 1,
                       final_variant_classification = variant_classification1,
                       cnap = as.character(cnap),
                       ccf_hat = as.character(ccf_hat)) |>
         dplyr::rename(cnaps = cnap, ccfs = ccf_hat)
     ) |>
   dplyr::mutate(variant_classification = dplyr::if_else(is.na(variant_classification2), variant_classification1,
                                                  paste(variant_classification1, variant_classification2, sep = ","))) |>
   dplyr::select(-variant_classification1, -variant_classification2)

   ## Separate lohs and deepdels
   deepdels_loh <- cnvs |>
     dplyr::filter(!variant_classification1 %in% amps_categories) |>
     dplyr::group_by(.data$hugo_symbol, .data$sample_id, .data$variant_classification1) |>
     dplyr::summarise(cnaps = paste0(cnap, collapse = ","),
                      ccfs = paste0(ccf_hat, collapse = ","),
                      variant_classification = paste0(variant_classification1, collapse = ","),
                      n_alts = dplyr::n()) |>
     dplyr::rename(final_variant_classification = variant_classification1) |>
     dplyr::ungroup()

   cnvs <- dplyr::bind_rows(all_amps, deepdels_loh) |>
     dplyr::mutate(variant_type = "CNV")

   dplyr::bind_rows(mutations, cnvs)
}

#' Adds Any Amp rows, which are rows with GAIN, AMP, HighAMP, and FocalHighAMP
#' in the alt column
#' @param .data A data frame with an 'alt' column
#'
#' @return A data frame with 'Any Amp' rows
#' @export
#'
#' @examples
#' add_any_amp(get_mbcp_dna_alts_summary())
add_any_amp <- function(.data) {

  stopifnot("'variant_classification' column missing in data" = "variant_classification" %in% colnames(.data))

  amps_categories <- get_mbcp_labels("dna_alts")[3:6]

  return(.data |>
    dplyr::bind_rows(.data |>
                       dplyr::filter(.data$variant_classification %in% amps_categories) |>
                       dplyr::mutate(variant_classification = "Any Amp")))

}

#' Get germline hits dataset
#'
#' @return A data frame with entries from the germline hits found in the MBCProject
#' @export
#'
#' @examples
#' get_mbcp_germline_hits()
get_mbcp_germline_hits <- function() {
  return(mbcp_germline_hits)
}
