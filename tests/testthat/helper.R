library(mockery)
library(dplyr)

set.seed(123)
mocked_clinical_data <- data.frame(
  sample_alias = sample(c("Alias1", "Alias2", "Alias3"), 100, replace = TRUE),
  sample_timepoint = sample(c("Timepoint1", "Timepoint2", "Timepoint3"), 100, replace = TRUE),
  patient_id = sample(1:100, 100),
  wes_sample_id = sample(1:100, 100),
  receptor_status_dx_all = sample(c("Positive", "Negative"), 100, replace = TRUE),
  hr_status_dx_all = sample(c("Positive", "Negative"), 100, replace = TRUE),
  bx_location = sample(c("Location1", "Location2"), 100, replace = TRUE),
  dx_histology = sample(c("Hist1", "Histology2"), 100, replace = TRUE),
  pam50_subtype = sample(c("Subtype1", "Subtype2"), 100, replace = TRUE),
  calc_time_to_mets_dx_days = sample(1:500, 100),
  calc_met_setting = sample(c("Setting1", "Setting2"), 100, replace = TRUE),
  calc_primary_treat = sample(c("Treatment1", "Treatment2"), 100, replace = TRUE),
  purity = runif(100), # Generate random numbers between 0 and 1
  timepoint = sample(1:5, 100, replace = TRUE),
  n_bypatient = sample(1:5, 100, replace = TRUE),
  bx_time_days =  sample(10:300, 100, replace = TRUE)
)

mocked_tpms <- matrix(runif(1000), nrow = 100, ncol = 10)
rownames(mocked_tpms) <- paste0("gene", 1:100)
colnames(mocked_tpms) <- paste0("sample", 1:10)

mocked_tpms_min <- matrix(c(10, 1, 2, 5, 30,
                15, 2, 3, 4, 15,
                20, 3, 4, 3, 20,
                25, 4, 5, 2, 25,
                30, 5, 6, 1, 10), nrow = 5, byrow = TRUE)
colnames(mocked_tpms_min) <- paste0("s", 1:5)
rownames(mocked_tpms_min) <- paste0("g", 1:5)

mocked_tpms_min_uq <- matrix(c(0.4, 0.25,  0.4, 1.25, 1.2,
                   0.6, 0.50,  0.6, 1.00,  0.6,
                   0.8, 0.75,  0.8, 0.75,  0.8,
                   1.0, 1.00,  1.0, 0.50,  1.0,
                   1.2, 1.25,  1.2, 0.25,  0.4), nrow = 5, byrow = T)
colnames(mocked_tpms_min_uq) <- paste0("s", 1:5)
rownames(mocked_tpms_min_uq) <- paste0("g", 1:5)

mocked_dna_alts <- data.frame(
  hugo_symbol = sample(paste0("g", 1:10), 80, replace = T),
  sample_id = sample(paste0("s", 20:40), 80, replace = T),
  ccf_hat = runif(80),
  variant_type = sample(c("MUT", "CNV"), 80, replace = T),
  variant_classification = sample(c("AMP", "HighAMP", "FocalHighAMP", "DeepDEL"), 80, replace = T )
)

mocked_mutations <- data.frame(
  Hugo_Symbol = sample(paste0("g", 1:10), 100, replace = T),
  Sample_ID = sample(paste0("s", 20:40), 100, replace = T),
  ccf_hat = sample(c(runif(80), rep(NA, 20))),
  Variant_Type = sample(c("DEL", "INS", "SNP"), 100, replace = T),
  Variant_Classification = sample(c("Loss_of_function", "Missense_Mutation", "Gain_of_function"), 100, replace = T)
) %>%
  dplyr::mutate(
    Genomic_Change = Sample_ID,
    HGVS_genomic_change = Sample_ID,
    HGVS_protein_change = Sample_ID,
    Protein_Change = Sample_ID
  )

loh <- sample(c(rep("", 90), rep(",LOH_LOF", 10)))

mocked_cnvs <- data.frame(
  hugo_symbol = sample(paste0("g", 1:10), 100, replace = T),
  sample_id = sample(paste0("s", 20:40), 100, replace = T),
  variant_type = "CNV",
  cnap = runif(100, min = -5, max= 5),
  variant_classification = sample(c("AMP", "HighAMP", "FocalHighAMP"), 100, replace = TRUE)
) %>%
  dplyr::mutate(variant_classification = ifelse(cnap < 0, "DeepDEL", variant_classification),
         variant_classification = paste0(variant_classification, loh))
