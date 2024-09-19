library(vroom)
library(janitor)
library(dplyr)
library(tidyr)
library(tibble)

load("data/mbcp_clinical_data.rda")

cnames <- vroom("data-raw/df_signatures_class_UQ.tsv", n_max = 1) %>%
  colnames()
cnames <- tolower(c("x1", cnames))

enrichment_class <- vroom("data-raw/df_signatures_class_UQ.tsv", col_names = cnames,
                           skip = 1) %>%
  select(-x1)

samples_dict <- vroom("data-raw/MBCp_sample_dictionary.tsv") %>%
  clean_names()

# mbcproject_0054_t1b to mbcproject_0054_t1b_2
samples_dict <- samples_dict %>%
  mutate(short_timepoint = tolower(short_timepoint_id_2),
         sample_timepoint = tolower(timepoint_id)) %>%
  select(sample_timepoint, short_timepoint)  %>%
  mutate(sample_timepoint = ifelse(sample_timepoint == "mbcproject_0054_t1b", "mbcproject_0054_t1b_2", sample_timepoint))

samples_dict <- samples_dict %>%
  inner_join(mbcp_clinical_data %>%
               select(sample_timepoint, sample_alias), by = "sample_timepoint")

enrichment_class <- enrichment_class %>%
  column_to_rownames("signature_id") %>%
  as.matrix()

samples <- samples_dict %>% pull(sample_alias, name = "short_timepoint")
## Just making sure all samples are present in the enrichment_class matrix
stopifnot(all(colnames(enrichment_class) %in% names(samples)))

## Unname it because colnames retain names and cause problems in testing
colnames(enrichment_class) <- unname(samples[colnames(enrichment_class)])

mbcp_enrichment_class <- enrichment_class
usethis::use_data(mbcp_enrichment_class, overwrite = TRUE, compress = "xz")
