#' Get label vectors for MBC project data
#'
#' @param set Character string specifying the labels to return. Options are: "receptor", "pam50",
#' "dna_alts", "mutation", "histology", "treatment", "metastasis", "expression"
#'
#' @return A named character vector
#' @export
#'
#' @examples
#' get_mbcp_labels("receptor")
get_mbcp_labels <- function(set) {

  stopifnot("set must be a character" = is.character(set))

  if(!set %in% pkg_env$set_names) {
    stop(paste0(set, " is not a valid label set name for MBC project"))
  }
  receptor <- c("HR+/HER2-", "HER2+", "TNBC", "MIXED", "UNKNOWN")
  names(receptor) <- receptor

  pam50 <- c("LumA", "LumB", "Her2", "Basal", "Normal", "NC")
  names(pam50) <- pam50

  metastasis <- c("NO METASTATIC\nDISEASE PRESENT", "METASTATIC\nDISEASE PRESENT",
                  "ABSTRACTION PENDING", "NOT FOUND IN RECORD")
  names(metastasis) <- c("NO_METASTATIC_DISEASE_PRESENT", "METASTATIC_DISEASE_PRESENT",
                        "ABSTRACTION_PENDING", "NOT_FOUND_IN_RECORD")

  histology <- c("ILC", "IDC", "DCIS", "MIXED IDLC", "ADENOCARCINOMA",
                 "OTHER RARE\nSUBTYPE", "ABSTRACTION PENDING", "NOT FOUND\nIN RECORD")
  names(histology) <- c("ILC", "IDC", "DCIS", "MIXED_IDLC", "ADENOCARCINOMA",
                       "OTHER_RARE_SUBTYPE", "ABSTRACTION_PENDING", "NOT_FOUND_IN_RECORD")

  dna_alts <- c("DeepDEL", "Any Amp", "AMP", "HighAMP", "FocalHighAMP", "Gain_of_function",
                "Loss_of_function", "Putative_loss_of_function", "Biallelic_inactivation",
                "Multi_hit", "Hotspot", "Missense_Mutation",  "Other_Mutation", "WT")

  names(dna_alts) <- dna_alts
  mutation <- c("MUT", "WT")
  names(mutation) <- c("MUT", "WT")

  treatment <- c("PRIMARY\nTX NAIVE", "MET TXED")
  names(treatment) <- c("PRIMARY_TX_NAIVE", "MET_TXED")

  categories <- c("Lower Decile", "Lower Quartile", "Mid", "Upper Quartile", "Upper Decile")
  names(categories) <- categories

  germline_alts <- c("DeepDEL", "LOF", "LOH", "VUS SNV", "None Detected", "No Tumor WES")
  names(germline_alts) <- germline_alts

  dna_signatures <- c("Clock-like", "APOBEC", "HR", "SBS.12.unknown")
  names(dna_signatures) <- dna_signatures
  # return vector with same name as set argument with switch
  return(switch(set,
                receptor = receptor,
                pam50 = pam50,
                dna_alts = dna_alts,
                mutation = mutation,
                histology = histology,
                treatment = treatment,
                metastasis = metastasis,
                categories = categories,
                germline_alts = germline_alts,
                dna_signatures = dna_signatures))
}
