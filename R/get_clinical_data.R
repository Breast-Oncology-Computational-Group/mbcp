utils::globalVariables(c("mbcp_clinical_data"))
#' Get clinical data with all the variables available
#'
#' @return A data frame with all the variables available
#' @export
#'
#' @examples
#' get_clinical_data_all_variables()
get_clinical_data_all_variables <- function() {
  return(mbcp_clinical_data)
}

#' Get clinical data with selected variables
#'
#' @return A data frame with sample level clinical information
#' @export
#' @importFrom rlang .data
#' @examples
#' get_clinical_data()
get_clinical_data <- function() {
  return(get_clinical_data_all_variables() |>
           dplyr::select("sample_alias", "patient_id", "sample_timepoint", "wes_sample_id",
                         "receptor_status_dx_all", "hr_status_dx_all",  "dx_histology", "pam50_subtype",
                         "bx_location", "bx_time_days", "calc_time_to_mets_dx_days", "calc_met_setting",
                         "calc_primary_treat", "purity", "timepoint", "n_bypatient"))

}


#' Get clinical data for RNA samples
#'
#' @param filter A character indicating the filter to apply to the data. Options are:
#' NULL, "normals", "purity", "tissues", "breast"
#' @param wes A logical indicating if clinical data for samples with RNA and WES should be returned
#' @return A data frame with sample level clinical information
#' @export
#'
#' @examples
#' get_clinical_data_rna()
#' get_clinical_data_rna("normals")
get_clinical_data_rna <- function(filter = NULL, wes = FALSE) {

  stopifnot("wes should be a logical value" = is.logical(wes))

  clinical_data <- get_clinical_data()

  if(wes) {
    clinical_data <- clinical_data |>
      dplyr::filter(!is.na(.data$wes_sample_id))
  }

  if(is.null(filter)) {
    return(clinical_data |>
             dplyr::filter(!is.na(.data$pam50_subtype)))
  } else if(filter == "normals") {
    return(clinical_data |>
             dplyr::filter(.data$pam50_subtype != "Normal"))
  } else if(filter == "purity") {
    return(clinical_data |>
             dplyr::filter(.data$pam50_subtype != "Normal" & .data$purity >= 0.272))
  } else if(filter == "tissues") {
    return(clinical_data |>
             dplyr::filter(.data$pam50_subtype != "Normal" &
                           !.data$bx_location %in% c("ASCENDING COLON","BONE","LUNG", "PLEURAL FLUID")))
  } else if(filter == "breast") {
    return(clinical_data |>
             dplyr::filter(.data$pam50_subtype != "Normal" &
                             .data$bx_location == "BREAST"))
  } else {
    stop(paste0(filter, " is not a valid filter argument"))
  }
}


#' Get clinical data for WES samples
#'
#' @return A data frame with sample level clinical information
#' @export
#'
#' @examples
#' get_clinical_data_wes()
get_clinical_data_wes <- function() {
    get_clinical_data() |>
    dplyr::filter(!is.na(.data$wes_sample_id))
}

#' Get clinical data for first samples per patient
#'
#' @param sample_type Character value to specify the type of samples to include. Options are:
#' \itemize{
#'  \item{any:}{gets the first timepoint sample per patient}
#'  \item{wes:}{gets the minimum timepoint sample per patient, only for WES samples}
#'  \item{rna:}{gets the minimum timepoint sample per patient, only for RNA samples}
#'  \item{wes_rna:}{gets the minimum timepoint sample per patient, for samples with WES and RNA}
#' }
#'
#' @return A data frame with sample level clinical information
#' @export
#'
#' @examples
#' get_clinical_data_first_samples(sample_type = "any")
#' get_clinical_data_first_samples(sample_type = "wes_rna")
get_clinical_data_first_samples <- function(sample_type = "any") {

  stopifnot("sample_type should be one of: any, wes, rna or wes_rna" = sample_type %in% c("any", "wes", "rna", "wes_rna"))

  if(sample_type == "any") {
    return(get_clinical_data() |>
             dplyr::filter(.data$timepoint == 1))
  } else if(sample_type == "wes") {
    return(get_clinical_data_wes() |>
             dplyr::group_by(.data$patient_id) |>
             dplyr::slice_min(.data$timepoint) |>
             dplyr::ungroup())
  } else if(sample_type == "rna") {
    return(get_clinical_data_rna() |>
             dplyr::group_by(.data$patient_id) |>
             dplyr::slice_min(.data$timepoint)|>
             dplyr::ungroup())
  } else if(sample_type == "wes_rna") {
    return(get_clinical_data_rna(wes = TRUE) |>
             dplyr::group_by(.data$patient_id) |>
             dplyr::slice_min(.data$timepoint)|>
             dplyr::ungroup())
  }
}

#' Get clinical data for latest samples per patient
#'
#' @inheritParams get_clinical_data_first_samples
#' @return A data frame with sample level clinical information
#' @export
#'
#' @examples
#' get_clinical_data_latest_samples(sample_type = "wes")
#' get_clinical_data_latest_samples(sample_type = "wes_rna")
get_clinical_data_latest_samples <- function(sample_type = "wes") {

  stopifnot("sample_type should be one of: any, wes, rna or wes_rna" = sample_type %in% c("any", "wes", "rna", "wes_rna"))

  if(sample_type == "any") {
    return(get_clinical_data() |>
             dplyr::filter(.data$timepoint == .data$n_bypatient))
  } else if(sample_type == "wes") {
    return(get_clinical_data_wes() |>
             dplyr::group_by(.data$patient_id) |>
             dplyr::slice_max(.data$timepoint) |>
             dplyr::ungroup())
  } else if(sample_type == "rna") {
    return(get_clinical_data_rna() |>
             dplyr::group_by(.data$patient_id) |>
             dplyr::slice_max(.data$timepoint) |>
             dplyr::ungroup())
  } else if(sample_type == "wes_rna") {
    return(get_clinical_data_rna(wes = TRUE) |>
             dplyr::group_by(.data$patient_id) |>
             dplyr::slice_max(.data$timepoint) |>
             dplyr::ungroup())
  }
}


#' Get clinical data for patients with multiple samples in reference/sample format
#'
#' @return A data frame with one entry per reference and comparison sample pair.
#' Reference samples are first samples with WES in patients where n_bypatient > 1 and they
#' are paired with the remaining WES samples from the same patient.
#' This set has an extra filter requiring the removal of five samples from the total
#' dataset.
#' @export
#'
#' @examples
#' get_clinical_data_multisamples()
get_clinical_data_multisamples <- function() {

  first_samples <- get_clinical_data_first_samples("wes") |>
    dplyr::filter(.data$n_bypatient > 1, !.data$wes_sample_id %in% samples_to_remove$sample_id)  |>
    dplyr::select("patient_id", time_to_mets_dx = "calc_time_to_mets_dx_days",
                  patient_hr_sub = "receptor_status_dx_all",
                  patient_histology = "dx_histology",
                  reference_biopsy = "bx_location", reference_pam50 = "pam50_subtype",
                  reference_hr_sub = "hr_status_dx_all",
                  reference_days = "bx_time_days", reference_met_set = "calc_met_setting",
                  reference_primary_treat = "calc_primary_treat",
                  reference_sample_id = "wes_sample_id",
                  reference_timepoint = "timepoint")

  other_samples <- get_clinical_data_wes() |>
    dplyr::filter(.data$n_bypatient > 1, !.data$wes_sample_id %in% samples_to_remove$sample_id) |>
    dplyr::anti_join(first_samples, by = c("patient_id", "wes_sample_id" = "reference_sample_id",
                                           "timepoint" = "reference_timepoint")) |>
    dplyr::select("patient_id", biopsy = "bx_location",
           hr_sub = "hr_status_dx_all", pam50 = "pam50_subtype", days = "bx_time_days",
           met_set = "calc_met_setting", treat_naive = "calc_primary_treat",
           sample_id = "wes_sample_id",
           timepoint = "timepoint")

  ## Not so sure if they should match the timepoints  in the clinical_data after samples removal
  ## or if they should get a new order. Re-ordering for now
  first_samples <- first_samples |>
    dplyr::mutate(reference_timepoint = 1)

  other_samples <- other_samples |>
    dplyr::group_by(.data$patient_id) |>
    dplyr::arrange(.data$timepoint, by_group = TRUE) |>
    dplyr::mutate(timepoint = dplyr::row_number() + 1) |>
    dplyr::ungroup()

  multisamples <- first_samples |>
    dplyr::inner_join(other_samples, by = "patient_id") |>
    dplyr::mutate(days_between = .data$days - .data$reference_days)

  ## Four months for de novo metastatic according to the paper
  multisamples <- multisamples |>
    dplyr::mutate(months = round(.data$days_between/30),
           de_novo = ifelse(is.na(.data$time_to_mets_dx), "ABSTRACTION_PENDING",
                            ifelse(.data$time_to_mets_dx <= 120, "YES", "NO"))) |>
    dplyr::mutate(months_cat = ifelse(months < 1, "< 1",
                               ifelse(months < 3, "1-3",
                                      ifelse(months < 12, "3-12",
                                             "> 12"))))

  return(multisamples)
}


