library(vroom)
library(dplyr)
library(janitor)

receptor_order <- c("HR+/HER2-", "HER2+", "TNBC", "MIXED", "UNKNOWN")
pam50_order <- c("LumA", "LumB", "Her2", "Basal", "Normal", "NC")
metastatic_order <- c("NO_METASTATIC_DISEASE_PRESENT", "METASTATIC_DISEASE_PRESENT",
                      "ABSTRACTION_PENDING", "NOT_FOUND_IN_RECORD")
histology_order <- c("ILC", "IDC", "DCIS", "MIXED_IDLC", "ADENOCARCINOMA",
                     "OTHER_RARE_SUBTYPE", "ABSTRACTION_PENDING", "NOT_FOUND_IN_RECORD")
primary_treat_order <- c("PRIMARY_TX_NAIVE", "MET_TXED")
dna_alts_order <- c("DeepDEL", "AMP", "HighAMP", "FocalHighAMP", "Gain_of_function",
                    "Loss_of_function", "Putative_loss_of_function", "Biallelic_inactivation",
                    "Multi_hit", "Hotspot", "Missense_Mutation", "Any Amp", "Other_Mutation", "WT")


### should we add rna sample id?
mbcp_clinical_data <- vroom("data-raw/sample_patient_metadata_HRHER2status_PAM50_MutSignatures_GenesInterest_wRNA.csv") %>%
  clean_names() %>%
  select(-x1) %>%
  mutate(patient_id = tolower(patient_id),
         wes_sample_id = tolower(wes_sample_id),
         sample_timepoint = tolower(sample_timepoint),
         receptor_status_dx_all = factor(receptor_status_dx_all, levels = receptor_order),
         ## exclude = NA by default
         pam50_subtype = factor(pam50_subtype, levels = pam50_order, exclude = NULL),
         calc_met_setting = factor(calc_met_setting, levels = metastatic_order, exclude = NULL),
         dx_histology = factor(dx_histology, levels = histology_order, exclude = NULL),
         calc_primary_treat = ifelse(calc_treatment_naive == "NO" &
                                       calc_met_setting == "METASTATIC_DISEASE_PRESENT", "MET_TXED",
                                     ifelse(calc_treatment_naive == "YES" &
                                              calc_met_setting == "NO_METASTATIC_DISEASE_PRESENT" &
                                              bx_location %in% c("BREAST","AXILLA","AXILLARY LYMPH NODE"),
                                            "PRIMARY_TX_NAIVE", NA)),
         calc_primary_treat = factor(calc_primary_treat, levels = primary_treat_order)) %>%
  select(-patient_id, -sample_id) %>%
  select(sample_alias = clean_sample_alias, patient_id = clean_participant, sample_timepoint, wes_sample_id,
         everything())


### Setting sample order by patient
mbcp_clinical_data <- mbcp_clinical_data %>%
  mutate(is_breast = ifelse(bx_location == "BREAST", 1, 2)) %>%
  group_by(patient_id) %>%
  arrange(bx_time_days, is_breast, .by_group = TRUE) %>%
  mutate(timepoint = row_number(),
         n_bypatient = n()) %>%
  ungroup() %>%
  select(-is_breast)

## MBCProject_0280
mbcp_clinical_data <- mbcp_clinical_data %>%
  mutate(timepoint = case_when(wes_sample_id == tolower("MBC-MBCProject_0280-Tumor-SM-AXGL6") ~ 2,
                           wes_sample_id == tolower("MBC-MBCProject_0280-Tumor-SM-AXGIU") ~ 1,
                           TRUE ~ timepoint))

## MBCProject_0003
mbcp_clinical_data <- mbcp_clinical_data %>%
  mutate(timepoint = case_when(wes_sample_id ==  tolower("MBC-MBCProject_0003-Tumor-SM-AZ5H9") ~ 2,
                           wes_sample_id ==  tolower("MBC-MBCProject_0003-Tumor-SM-AZ5HV") ~ 1,
                           TRUE ~ timepoint))

## MBCProject_0734
mbcp_clinical_data <- mbcp_clinical_data %>%
  mutate(timepoint = case_when(wes_sample_id == tolower("MBC-MBCProject_0734-Tumor-SM-CGM1H") ~ 2,
                           wes_sample_id == tolower("MBC-MBCProject_0734-Tumor-SM-CGMCS") ~ 1,
                           TRUE ~ timepoint))

usethis::use_data(mbcp_clinical_data, overwrite = TRUE, compress = "xz")
