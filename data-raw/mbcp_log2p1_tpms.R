library(vroom)
library(dplyr)
library(tibble)
library(janitor)
library(stringr)
library(tidyr)

load("data/mbcp_clinical_data.rda")
samples_dict <- vroom("data-raw/MBCp_sample_dictionary.tsv") %>%
  clean_names()

# mbcproject_0054_t1b to mbcproject_0054_t1b_2
samples_dict <- samples_dict %>%
  mutate(short_timepoint = tolower(short_timepoint_id_2),
         sample_timepoint = tolower(timepoint_id)) %>%
  select(sample_timepoint, short_timepoint)  %>%
  mutate(sample_timepoint = ifelse(sample_timepoint == "mbcproject_0054_t1b", "mbcproject_0054_t1b_2", sample_timepoint))

mbcp_log2p1_tpms <- vroom("data-raw/df_rsem_log2tpm_expressed2.csv") %>%
  clean_names()%>%
  select(-x1, -gene_id) %>%
  mutate(hugo_symbol = unlist(lapply(strsplit(hgnc_symbol_gene_id, split = "\\s*\\("), "[[", 1))) %>%
  select(-hgnc_symbol_gene_id) %>%
  select(hugo_symbol, everything()) %>%
  filter(hugo_symbol != "") %>%
  pivot_longer(cols = -hugo_symbol, names_to = "short_timepoint", values_to = "log2_tpm") %>%
  inner_join(samples_dict, by = "short_timepoint") %>%
  inner_join(mbcp_clinical_data %>% select(sample_timepoint, sample_alias), by = "sample_timepoint") %>%
  select(-short_timepoint, -sample_timepoint) %>%
  pivot_wider(names_from = sample_alias, values_from = log2_tpm)

mbcp_log2p1_tpms <- mbcp_log2p1_tpms %>%
  column_to_rownames("hugo_symbol") %>% as.matrix()

usethis::use_data(mbcp_log2p1_tpms, overwrite = TRUE, compress = "xz")
