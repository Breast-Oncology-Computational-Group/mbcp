library(vroom)
library(janitor)
library(dplyr)
library(tibble)

load("data/mbcp_clinical_data.rda")

cnames <- vroom("data-raw/MBCP_Hallmark_NES_UQ seed.csv", n_max = 1) %>%
  colnames()
cnames[1] <- "short_timepoint"
cnames <- gsub("NES_", "", cnames)
enrichment_scores <- vroom("data-raw/MBCP_Hallmark_NES_UQ seed.csv", col_names = cnames,
                        skip = 1)

cnames2 <- vroom("data-raw/MBCP_c2_breast_NES_UQ seed.csv", n_max = 1) %>%
  colnames()
cnames2[1] <- "short_timepoint"
cnames2 <- gsub("NES_", "", cnames2)
enrichment_scores2 <- vroom("data-raw/MBCP_c2_breast_NES_UQ seed.csv", col_names = cnames2,
                           skip = 1)


## removing HER2MUTvsGFP_3_UP and RTK_ACT_UP cause they're repeated
enrichment_scores <- enrichment_scores %>%
  select(-HER2MUTvsGFP_3_UP, -RTK_ACT_UP) %>%
  inner_join(enrichment_scores2, by = "short_timepoint")
cnames <- c(cnames, cnames2[-1])
samples_dict <- vroom("data-raw/MBCp_sample_dictionary.tsv") %>%
  clean_names()

# mbcproject_0054_t1b to mbcproject_0054_t1b_2
samples_dict <- samples_dict %>%
  mutate(short_timepoint = tolower(short_timepoint_id_2),
         sample_timepoint = tolower(timepoint_id)) %>%
  select(sample_timepoint, short_timepoint)  %>%
  mutate(sample_timepoint = ifelse(sample_timepoint == "mbcproject_0054_t1b", "mbcproject_0054_t1b_2", sample_timepoint))

enrichment_scores <- enrichment_scores %>%
  mutate(short_timepoint = tolower(short_timepoint)) %>%
  inner_join(samples_dict, by = "short_timepoint") %>%
  inner_join(mbcp_clinical_data %>%
               select(sample_timepoint, sample_alias), by = "sample_timepoint") %>%
  select(sample_alias, all_of(cnames[-1])) %>%
  column_to_rownames("sample_alias") %>%
  as.matrix() %>%
  t()

mbcp_enrichment_scores <- enrichment_scores
usethis::use_data(mbcp_enrichment_scores, overwrite = TRUE, compress = "xz")
