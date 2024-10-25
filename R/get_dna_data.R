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
      dplyr::filter(.data$ccf_hat >= min_ccf)
  }

  return(mutations)
}

#' Get CNV in MBC project
#' @param samples Character vector of samples to filter variations by sample
#' @param genes Character vector of hugo symbols to filter variations by gene
#' @export
#'
#' @examples
#' get_mbcp_cnvs()
get_mbcp_cnvs <- function(samples = NULL, genes = NULL) {

  stopifnot("samples should be a character vector" =  is.null(samples) | is.character(samples),
            "genes should be a character vector" = is.null(genes) | is.character(genes))

  all_samples <- unique(mbcp_cnvs$sample_id)

  if(any(!samples %in% all_samples)) {
    warning("some samples were not found in the mbcp_cnvs dataset")
  }

  cnvs <- mbcp_cnvs |>
    dplyr::filter(is.null(samples) | .data$sample_id %in% samples,
                  is.null(genes) | .data$hugo_symbol %in% genes)

    return(cnvs)
}

#' Get dna alterations (CNVs and mutations) in the MBC project.
#' @inheritParams get_mbcp_mutations
#' @param variant_type Used to specify the type of variant alterations to add in the summary. Must
#'  be one of: "both", "CNV", "MUT"
#'
#'  One entry per sample and gene for mutations with Multiple_hit, if multiple entries
#'  One entry per sample and gene for amplifications,
#' @export
#'
#' @examples
#' get_mbcp_dna_alts_summary(genes = c("ESR1", "TP53"), min_ccf = 0.75)
get_mbcp_dna_alts_summary <- function(samples = NULL, genes = NULL, min_ccf = NULL,
                                      variant_type = "both") {

  stopifnot("samples should be a character vector" =  is.null(samples) | is.character(samples),
            "genes should be a character vector" = is.null(genes) | is.character(genes),
            "min_ccf should be a numeric value" =  is.null(min_ccf) | is.numeric(min_ccf) & length(min_ccf) == 1,
            "variant_type should be one of MUT, CNV, both" = variant_type %in% c("CNV", "MUT", "both"))

  amps_categories <- get_mbcp_labels("dna_alts")[3:6]

  if(variant_type == "both" | variant_type == "MUT") {
    mutations <- get_mbcp_mutations(samples = samples, genes = genes, min_ccf =  min_ccf) |>
       dplyr::select(hugo_symbol = "Hugo_Symbol", sample_id = "Sample_ID",
                     genomic_change = "Genomic_Change", protein_change = "Protein_Change",
                     hgvs_genomic_change = "HGVS_genomic_change",
                     hgvs_protein_change = "HGVS_protein_change", "ccf_hat",
                     variant_type = "Variant_Type", variant_classification = "Variant_Classification") |>
       dplyr::group_by(.data$hugo_symbol, .data$sample_id) |>
       dplyr::summarise(n_alts = dplyr::n(),
                        ccf_hat = ifelse(all(is.na(.data$ccf_hat)), NA, max(.data$ccf_hat, na.rm = T)),
                        all_variants = paste0(.data$variant_classification, collapse = ","),
                        all_genomic_change = paste0(.data$genomic_change, collapse = ","),
                        all_hgvs_genomic_change = paste0(.data$hgvs_genomic_change, collapse = ","),
                        all_protein_change = paste0(.data$protein_change, collapse = ","),
                        all_hgvs_protein_change = paste0(.data$hgvs_protein_change, collapse = ","),
                        all_ccfs = paste0(.data$ccf_hat, collapse = ",")) |>
      dplyr::ungroup() |>
      dplyr::mutate(variant_classification = ifelse(.data$n_alts == 1, .data$all_variants, "Multi_hit"),
                     variant_type = "MUT")
  }

  if(variant_type == "both" | variant_type == "CNV") {

    cnvs <- get_mbcp_cnvs(samples = samples, genes = genes) |>
      dplyr::select("hugo_symbol", "sample_id", "cnap",  "variant_type",
                     "variant_classification") |>
      tidyr::separate_wider_delim("variant_classification", ",", names = c("variant_classification1", "variant_classification2"),
                                  too_few = "align_start")

     ### Pick the most severe amplification as variant_classification
     amps <- cnvs |>
      dplyr::filter(.data$variant_classification1 %in% amps_categories)

     amps_count <- amps |>
      dplyr::count(.data$hugo_symbol, .data$sample_id)

     amps_one <-  amps |>
       dplyr::semi_join(amps_count |>
                          dplyr::filter(.data$n == 1),
                        by = c("sample_id", "hugo_symbol")) |>
       dplyr::mutate(n_alts = 1,
                     variant_classification = .data$variant_classification1) |>
       dplyr::select(dplyr::everything(), max_cnap = "cnap", all_cnaps = "cnap") |>
       dplyr::mutate(all_cnaps = as.character(round(.data$all_cnaps, 4)))

     if(any(amps_count$n > 1)) {
       all_amps <- amps |>
         dplyr::anti_join(amps_count |>
                            dplyr::filter(.data$n == 1),
                          by = c("sample_id", "hugo_symbol")) |>
         dplyr::mutate(factor_vc = factor(.data$variant_classification1, levels = amps_categories, ordered = TRUE)) |>
         dplyr::group_by(.data$hugo_symbol, .data$sample_id) |>
         dplyr::summarise(n_alts = dplyr::n(),
                          variant_classification1 = paste0(.data$variant_classification1, collapse = ","),
                          variant_classification2 = unique(.data$variant_classification2),
                          max_cnap = max(.data$cnap),
                          all_cnaps = paste0(round(.data$cnap, 4), collapse = ","),
                          variant_classification = max(.data$factor_vc)) |>
         dplyr::ungroup() |>
         dplyr::bind_rows(amps_one)
     } else {
       all_amps <- amps_one
     }

     all_amps <- all_amps |>
     dplyr::mutate(all_variants = dplyr::if_else(is.na(.data$variant_classification2), .data$variant_classification1,
                                                    paste(.data$variant_classification1, .data$variant_classification2, sep = ",")),
                   variant_classification = ifelse(is.na(.data$variant_classification2), .data$variant_classification,
                                  paste(.data$variant_classification, .data$variant_classification2, sep = ","))) |>
     dplyr::select(-"variant_classification1", -"variant_classification2")

     ## Separate lohs and deepdels
     deepdels_loh <- cnvs |>
       dplyr::filter(!.data$variant_classification1 %in% amps_categories)

     deepdels_count <- deepdels_loh |>
       dplyr::count(.data$hugo_symbol, .data$sample_id)

     deepdels_one <- deepdels_loh |>
       dplyr::semi_join(deepdels_count |>
                          dplyr::filter(.data$n == 1),
                        by = c("sample_id", "hugo_symbol")) |>
       dplyr::mutate(n_alts = 1,
                     variant_classification = .data$variant_classification1) |>
       dplyr::select(dplyr::everything(), max_cnap = "cnap", all_cnaps = "cnap") |>
       dplyr::mutate(all_cnaps = as.character(round(.data$all_cnaps, 4)))

     if(any(deepdels_count$n > 1)) {
        all_deepdels <- deepdels_loh |>
           dplyr::anti_join(deepdels_count |>
                              dplyr::filter(.data$n == 1),
                            by = c("sample_id", "hugo_symbol")) |>
         dplyr::group_by(.data$hugo_symbol, .data$sample_id, .data$variant_classification1) |>
         dplyr::summarise(all_cnaps = paste0(round(.data$cnap, 4), collapse = ","),
                          n_alts = dplyr::n(),
                          max_cnap = min(.data$cnap),
                          variant_classification = unique(.data$variant_classification1),
                          variant_classification2 = unique(.data$variant_classification2),
                          variant_classification1 = paste0(.data$variant_classification1, collapse = ",")) |>
         dplyr::ungroup() |>
         dplyr::bind_rows(deepdels_one)

     } else {
       all_deepdels <- deepdels_one
     }

     all_deepdels <- all_deepdels |>
       dplyr::mutate(all_variants = dplyr::if_else(is.na(.data$variant_classification2), .data$variant_classification1,
                                                   paste(.data$variant_classification1, .data$variant_classification2, sep = ",")),
                     variant_classification = ifelse(is.na(.data$variant_classification2), .data$variant_classification,
                                                     paste(.data$variant_classification, .data$variant_classification2, sep = ","))) |>
       dplyr::select(-"variant_classification1", -"variant_classification2")

     cnvs <- dplyr::bind_rows(all_amps, all_deepdels) |>
       dplyr::mutate(variant_type = "CNV")
  }

  switch (variant_type,
    both = dplyr::bind_rows(mutations, cnvs),
    MUT = mutations,
    CNV = cnvs
  )
}

#' Adds Any Amp rows, which are rows with GAIN, AMP, HighAMP, and FocalHighAMP
#' in the alt column
#' @param .data A data frame with an 'alt' column
#'
#' @return A data frame with 'Any Amp' rows
#' @export
#'
#' @examples
#' add_any_amp(get_mbcp_dna_alts_summary(genes = c("ESR1", "TP53"), min_ccf = 0.75))
add_any_amp <- function(.data) {

  stopifnot("'variant_classification' column missing in data" = "variant_classification" %in% colnames(.data))

  if("Any Amp" %in% unique(.data$variant_classification)) {
    stop("dataset already contains Any Amp")
  }

  amps_categories <- get_mbcp_labels("dna_alts")[3:6]
  ### add warning if any amp is there
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
