library(mockery)

set.seed(123)
mocked_clinical_data <- data.frame(
  sample_alias = sample(c("Alias1", "Alias2", "Alias3"), 100, replace = TRUE),
  sample_timepoint = sample(c("Timepoint1", "Timepoint2", "Timepoint3"), 100, replace = TRUE),
  participant_id = sample(1:100, 100),
  wes_sample_id = sample(1:100, 100),
  receptor_status_dx_all = sample(c("Positive", "Negative"), 100, replace = TRUE),
  hr_status_dx_all = sample(c("Positive", "Negative"), 100, replace = TRUE),
  bx_location = sample(c("Location1", "Location2"), 100, replace = TRUE),
  dx_histology = sample(c("Hist1", "Histology2"), 100, replace = TRUE),
  pam50_subtype = sample(c("Subtype1", "Subtype2"), 100, replace = TRUE),
  calc_time_to_mets_dx_days = sample(1:500, 100),
  calc_met_setting = sample(c("Setting1", "Setting2"), 100, replace = TRUE),
  calc_primary_treat = sample(c("Treatment1", "Treatment2"), 100, replace = TRUE),
  purity = runif(100) # Generate random numbers between 0 and 1
)

mocked_tpms <- matrix(runif(1000), nrow = 100, ncol = 10)
rownames(mocked_tpms) <- paste0("gene", 1:100)
colnames(mocked_tpms) <- paste0("sample", 1:10)

