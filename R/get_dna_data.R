utils::globalVariables(c("mbcp_dna_alts"))
#' Get DNA alts found in MBC project
#'
#' @param samples Character vector of samples to filter variations by sample
#' @param genes Character vector of hugo symbols to filter variations by gene
#' @param min_ccf Minimum ccf
#' @return A data frame with DNA alts. Variant types: MUT (snv) y CNV
#' @export
#'
#' @examples
#' get_mbcp_dna_alts()
get_mbcp_dna_alts <- function(samples = NULL, genes = NULL, min_ccf = 0) {

  stopifnot("samples should be a character vector" = is.character(samples) | is.null(samples),
            "genes should be a character vector" = is.character(genes) | is.null(genes),
            "min_ccf should be a numeric value" = is.numeric(min_ccf) & length(min_ccf) == 1)

  all_samples <- unique(mbcp_dna_alts$sample_id)

  if(any(!samples %in% all_samples)) {
    warning("some samples were not found in the mbcp_dna_alts dataset")
  }

  return(mbcp_dna_alts |>
           dplyr::filter(.data$sample_id %in% samples | is.null(samples),
                         .data$hugo_symbol %in% genes | is.null(genes),
                         .data$ccf_hat >= min_ccf))
}

#' Get SNVs found in MBC project
#'
#' @inheritParams get_mbcp_dna_alts
#' @return A data frame with SNVs
#' @export
#'
#' @examples
#' get_mbcp_snvs()
get_mbcp_snvs <- function(samples = NULL, genes = NULL, min_ccf = 0) {
  return(get_mbcp_dna_alts(samples = samples, genes = genes, min_ccf =  min_ccf ) |>
           dplyr::filter(.data$variant_type == "MUT"))
}

#' Get CNVs found in MBC project
#'
#' @inheritParams get_mbcp_dna_alts
#' @return A data frame with CNVs
#' @export
#'
#' @examples
#' get_mbcp_cnvs()
get_mbcp_cnvs <- function(samples = NULL, genes = NULL, min_ccf = 0) {
  return(get_mbcp_dna_alts(samples = samples, genes = genes, min_ccf =  min_ccf ) |>
           dplyr::filter(.data$variant_type == "CNV"))
}

#' Adds Any Amp rows, which are rows with AMP, HighAMP, and FocalHighAMP
#' in the alt column
#' @param .data A data frame with an 'alt' column
#'
#' @return A data frame with 'Any Amp' rows
#' @export
#'
#' @examples
#' add_any_amp(get_mbcp_dna_alts())
add_any_amp <- function(.data) {

  stopifnot("'variant_classification' column missing in data" = "variant_classification" %in% colnames(.data))

  return(.data |>
    dplyr::bind_rows(.data |>
                       dplyr::filter(.data$variant_classification %in%  c("AMP", "HighAMP", "FocalHighAMP")) |>
                       dplyr::mutate(variant_classification = "Any Amp")))

}
