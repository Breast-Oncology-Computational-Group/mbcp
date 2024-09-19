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
get_mbcp_enrichment_scores <- function(samples = NULL, signatures = NULL) {

  if(!is.null(samples)) {
    stopifnot("samples should be a character vector" = is.character(samples),
              "invalid sample selection" = all(samples %in% colnames(mbcp_enrichment_scores)))
  }

  if(!is.null(signatures)) {
    stopifnot("signatures should be a character vector" = is.character(signatures),
              "invalid signature selection" = all(signatures %in% rownames(mbcp_enrichment_scores)))
  }
  return(mbcp_enrichment_scores[signatures, samples, drop = FALSE])
}

#' Get MBCproject enrichment classification
#'
#' @inheritParams get_mbcp_enrichment_scores
#' @return A character matrix with classification of normalized enrichment scores for hallmark and selected sets in the MBCproject
#' @export
#'
#' @examples
#' get_mbcp_enrichment_class()
get_mbcp_enrichment_class <- function(samples = NULL, signatures = NULL) {
  if(!is.null(samples)) {
    stopifnot("samples should be a character vector" = is.character(samples),
              "invalid sample selection" = all(samples %in% colnames(mbcp_enrichment_class)))
  }

  if(!is.null(signatures)) {
    stopifnot("signatures should be a character vector" = is.character(signatures),
              "invalid signature selection" = all(signatures %in% rownames(mbcp_enrichment_class)))
  }
  return(mbcp_enrichment_class[signatures, samples, drop = FALSE])
}
