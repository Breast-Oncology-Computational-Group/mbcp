utils::globalVariables(c("mbcp_enrichment_scores", "mbcp_enrichment_class"))
#' Get MBCproject enrichment scores
#'
#' @param signatures Character vector of signature ids to filter rows
#' @param samples Character vector of samples to filter columns
#' @param signatures_set Name of signatures set to select a subset of signatures. Options are:
#' hallmarks, breast, rtk
#' Only one of `signatures` or `signatures_set` must be provided
#' @return A numeric matrix with normalized enrichment scores for hallmark and selected sets in the MBCproject
#' @export
#'
#' @examples
#' get_mbcp_expression_signature_scores()
get_mbcp_expression_signature_scores <- function(signatures = NULL, samples = NULL, signatures_set = NULL) {
  scores <- mbcp_enrichment_scores

  rlang::check_exclusive(signatures, signatures_set, .require = FALSE)

  if(!is.null(samples)) {
    stopifnot("samples should be a character vector" = is.character(samples),
              "invalid sample selection" = all(samples %in% colnames(scores)))
  }

  if(!is.null(signatures)) {
    stopifnot("signatures should be a character vector" = is.character(signatures),
              "invalid signature selection" = all(signatures %in% rownames(scores)))
  } else if(!is.null(signatures_set)) {
    if(!signatures_set %in% pkg_env$signatures_set) {
      stop(paste0(signatures_set, " is not a valid signatures set"))
    }

    signatures <- switch (signatures_set,
            hallmarks = grep("HALLMARK", rownames(scores), value = T),
            rtk = c("HER2MUTvsGFP_3_UP", "RTK_ACT_UP"),
            breast = setdiff(rownames(scores), c(grep("HALLMARK", rownames(scores), value = T), "HER2MUTvsGFP_3_UP", "RTK_ACT_UP"))
    )
  }

  scores <- if(!is.null(signatures)) scores[signatures, , drop = F] else scores
  scores <- if(!is.null(samples)) scores[, samples, drop = F] else scores
  return(scores)
}

#' Get MBCproject enrichment classification
#'
#' @inheritParams get_mbcp_expression_signature_scores
#' @return A character matrix with classification of normalized enrichment scores for hallmark signatures and other selected sets.
#' Categories are: Lower Decile, Lower Quartile, Upper Quartile, and Upper Decile. `NA` values for entries that do not fall into these categories.
#' @export
#'
#' @examples
#' get_mbcp_expression_signature_class()
get_mbcp_expression_signature_class <- function(signatures = NULL, samples = NULL) {
  class <- mbcp_enrichment_class

  if(!is.null(samples)) {
    stopifnot("samples should be a character vector" = is.character(samples),
              "invalid sample selection" = all(samples %in% colnames(class)))
  }

  if(!is.null(signatures)) {
    stopifnot("signatures should be a character vector" = is.character(signatures),
              "invalid signature selection" = all(signatures %in% rownames(class)))
  }

  class <- if(!is.null(signatures)) class[signatures, , drop = F] else class
  class <- if(!is.null(samples)) class[, samples, drop = F] else class
  return(class)
}

#' Get list of signatures and their associated genes studied in the MBCProject
#'
#' @return A data frame with signature, gene, and category columns for genes functionally
#' associated with their signatures.
#' @export
#'
#' @examples
#' get_mbcp_expression_signatures()
get_mbcp_expression_signatures <- function() {
  return(mbcp_genes_signatures |>
           dplyr::select(.data$signature_id, .data$hugo_symbol, .data$category))
}
