utils::globalVariables(c("mbcp_log2p1_tpms", "mbcp_log2p1_tpms_class"))
#' Get TPMs
#'
#' @param genes Character vector of hugo symbols to filter rows
#' @param samples Character vector of samples to filter columns
#' @return A matrix of TPMs with genes as rows and samples as columns with sample_alias as column names
#' @export
#'
#' @examples
#' get_mbcp_log2p1_tpms()
get_mbcp_log2p1_tpms <- function(genes = NULL, samples = NULL) {
  tpms <- mbcp_log2p1_tpms

  if(!is.null(samples)) {
    stopifnot("samples should be a character vector" = is.character(samples),
              "invalid sample selection" = all(samples %in% colnames(tpms)))
  }

  if(!is.null(genes)) {
    stopifnot("genes should be a character vector" = is.character(genes),
              "invalid gene selection" = all(genes %in% rownames(tpms)))
  }

  tpms <- if(!is.null(genes)) tpms[genes, , drop = F] else tpms
  tpms <- if(!is.null(samples)) tpms[, samples, drop = F] else tpms
  return(tpms)
}

#' Get upper quartile normalization TPMs
#'
#' @inheritParams get_mbcp_log2p1_tpms
#' @return A matrix of UQ TPMs with genes as rows and samples as columns with sample_alias as column names
#' @export
#'
#' @examples
#' get_mbcp_log2p1_tpms_uq()
get_mbcp_log2p1_tpms_uq <- function(genes = NULL, samples = NULL) {
  tpms <- get_mbcp_log2p1_tpms(genes = NULL, samples = NULL)

  if(!is.null(samples)) {
    stopifnot("samples should be a character vector" = is.character(samples),
              "invalid sample selection" = all(samples %in% colnames(tpms)))
  }

  if(!is.null(genes)) {
    stopifnot("genes should be a character vector" = is.character(genes),
              "invalid gene selection" = all(genes %in% rownames(tpms)))
  }

  tpms_uq <- apply(tpms, 2, \(x) x/stats::quantile(x, 0.75))

  tpms_uq <- if(!is.null(genes)) tpms_uq[genes, , drop = F] else tpms_uq
  tpms_uq <- if(!is.null(samples)) tpms_uq[, samples, drop = F] else tpms_uq

  return(tpms_uq)
}

#' Get gene expression classification
#'
#' @inheritParams get_mbcp_log2p1_tpms
#' @return A character matrix with classification of log2 tpms.
#' Categories are: Lower Decile, Lower Quartile, Upper Quartile, and Upper Decile. `NA` values for entries that do not fall into these categories.
#'
#' @export
#'
#' @examples
#' get_mbcp_log2p1_tpms_class()
get_mbcp_log2p1_tpms_class <- function(genes = NULL, samples = NULL) {

  tpms_class <- mbcp_log2p1_tpms_class

  if(!is.null(samples)) {
    stopifnot("samples should be a character vector" = is.character(samples),
              "invalid sample selection" = all(samples %in% colnames(tpms_class)))
  }

  if(!is.null(genes)) {
    stopifnot("genes should be a character vector" = is.character(genes),
              "invalid gene selection" = all(genes %in% rownames(tpms_class)))
  }

  tpms_class <- if(!is.null(genes)) tpms_class[genes, , drop = F] else tpms_class
  tpms_class <- if(!is.null(samples)) tpms_class[, samples, drop = F] else tpms_class
}

#' Get matrix with mbcp data in long format
#'
#' @param mbcp_matrix Numeric matrix with TPM values or signature data, or character matrix with their classification
#' @param id_column Character value to identify rownames in matrix. Options are: "hugo_symbol" and "signature_id"
#'
#' @return A dataframe with values in long format. The number of rows equals \var{n}x\var{m}, the dimensions
#' of \code{mbcp_matrix}, with three columns:
#' \itemize{
#'    \item{Value of \code{id_column}:} {Contains rownames in the \code{mbcp_matrix}}
#'    \item{\var{sample_alias}:} {Contains colnames in the \code{mbcp_matrix}}
#'    \item{\var{value}:} {The \[\code{id_column}, \var{sample_alias}\] value in the \code{mbcp_matrix}}
#' }
#
#' @export
#'
#' @examples
#' to_longer(get_mbcp_log2p1_tpms_uq(), "hugo_symbol")
to_longer <- function(mbcp_matrix, id_column) {

  stopifnot("mbcp_matrix should be a numeric or a character matrix" = is.matrix(mbcp_matrix),
            "mbcp_matrix should be a numeric or a character matrix" = is.numeric(mbcp_matrix) | is.character(mbcp_matrix))

  if(!id_column %in% c("hugo_symbol", "signature_id")) {
    stop(paste0(id_column, " is not a valid id_column for mbcp"))
  }

  mbcp_df <- mbcp_matrix |>
    as.data.frame() |>
    tibble::rownames_to_column({{id_column}}) |>
    tidyr::pivot_longer(cols = -{{id_column}}, names_to = "sample_alias", values_to = "value")
  return(mbcp_df)
}
