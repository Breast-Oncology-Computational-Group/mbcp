utils::globalVariables(c("mbcp_enrichment_scores", "mbcp_enrichment_class"))
#' Get MBCproject enrichment scores
#'
#' @param signatures Character vector of signature ids to filter rows
#' @param samples Character vector of samples to filter columns
#' @return A numeric matrix with normalized enrichment scores for hallmark and selected sets in the MBCproject
#' @export
#'
#' @examples
#' get_mbcp_enrichment_scores()
get_mbcp_enrichment_scores <- function(signatures = NULL, samples = NULL) {
  scores <- mbcp_enrichment_scores

  if(!is.null(samples)) {
    stopifnot("samples should be a character vector" = is.character(samples),
              "invalid sample selection" = all(samples %in% colnames(scores)))
  }

  if(!is.null(signatures)) {
    stopifnot("signatures should be a character vector" = is.character(signatures),
              "invalid signature selection" = all(signatures %in% rownames(scores)))
  }

  scores <- if(!is.null(signatures)) scores[signatures, , drop = F] else scores
  scores <- if(!is.null(samples)) scores[, samples, drop = F] else scores
  return(scores)
}

#' Get MBCproject enrichment classification
#'
#' @inheritParams get_mbcp_enrichment_scores
#' @return A character matrix with classification of normalized enrichment scores for hallmark signatures and other selected sets.
#' Categories are: Lower Decile, Lower Quartile, Upper Quartile, and Upper Decile. `NA` values for entries that do not fall into these categories.
#' @export
#'
#' @examples
#' get_mbcp_enrichment_class()
get_mbcp_enrichment_class <- function(signatures = NULL, samples = NULL) {
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
