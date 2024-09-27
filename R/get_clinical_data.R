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
           dplyr::select(.data$sample_alias, .data$participant_id, .data$sample_timepoint, .data$wes_sample_id,
                         .data$receptor_status_dx_all, .data$hr_status_dx_all,  .data$bx_location, .data$dx_histology, .data$pam50_subtype, .data$calc_time_to_mets_dx_days,
                         .data$calc_met_setting, .data$calc_primary_treat, .data$purity, .data$timepoint, .data$n_participant))

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
             dplyr::group_by(.data$participant_id) |>
             dplyr::slice_min(.data$timepoint) |>
             dplyr::ungroup())
  } else if(sample_type == "rna") {
    return(get_clinical_data_rna() |>
             dplyr::group_by(.data$participant_id) |>
             dplyr::slice_min(.data$timepoint)|>
             dplyr::ungroup())
  } else if(sample_type == "wes_rna") {
    return(get_clinical_data_rna(wes = TRUE) |>
             dplyr::group_by(.data$participant_id) |>
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
#' get_clinical_data_latest_samples(sample_type = "any")
#' get_clinical_data_latest_samples(sample_type = "wes_rna")
get_clinical_data_latest_samples <- function(sample_type = "any") {

  stopifnot("sample_type should be one of: any, wes, rna or wes_rna" = sample_type %in% c("any", "wes", "rna", "wes_rna"))

  if(sample_type == "any") {
    return(get_clinical_data() |>
             dplyr::filter(.data$timepoint == .data$n_participant))
  } else if(sample_type == "wes") {
    return(get_clinical_data_wes() |>
             dplyr::group_by(.data$participant_id) |>
             dplyr::slice_max(.data$timepoint) |>
             dplyr::ungroup())
  } else if(sample_type == "rna") {
    return(get_clinical_data_rna() |>
             dplyr::group_by(.data$participant_id) |>
             dplyr::slice_max(.data$timepoint) |>
             dplyr::ungroup())
  } else if(sample_type == "wes_rna") {
    return(get_clinical_data_rna(wes = TRUE) |>
             dplyr::group_by(.data$participant_id) |>
             dplyr::slice_max(.data$timepoint) |>
             dplyr::ungroup())
  }
}
