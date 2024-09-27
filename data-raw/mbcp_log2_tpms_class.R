library(vroom)
library(dplyr)
library(tidyr)
library(janitor)
library(tibble)

load("data/mbcp_clinical_data.rda")

cnames <- vroom("data-raw/df_genes_class_Rank_UQ.tsv", n_max = 1) %>%
  colnames()
cnames <- tolower(c("x1", cnames))

samples_dict <- vroom("data-raw/MBCp_sample_dictionary.tsv") %>%
  clean_names()
gene_dict <- vroom("data-raw/hg19_GRCh37p13_dictionary.csv")
# mbcproject_0054_t1b to mbcproject_0054_t1b_2
samples_dict <- samples_dict %>%
  mutate(short_timepoint = tolower(short_timepoint_id_2),
         sample_timepoint = tolower(timepoint_id)) %>%
  select(sample_timepoint, short_timepoint)  %>%
  mutate(sample_timepoint = ifelse(sample_timepoint == "mbcproject_0054_t1b", "mbcproject_0054_t1b_2", sample_timepoint))

mbcp_log2_tpms_class <-vroom("data-raw/df_genes_class_Rank_UQ.tsv", col_names = cnames,
                       skip = 1) %>%
  select(-x1) %>%
  inner_join(gene_dict %>% select(gene_id, hugo_symbol), by = c("gene_id")) %>%
  select(-gene_id) %>%
  select(hugo_symbol, everything()) %>%
  filter(hugo_symbol != "") %>%
  pivot_longer(cols = -hugo_symbol, names_to = "short_timepoint", values_to = "log2_tpm") %>%
  inner_join(samples_dict, by = "short_timepoint") %>%
  inner_join(mbcp_clinical_data %>% select(sample_timepoint, sample_alias), by = "sample_timepoint") %>%
  select(-short_timepoint, -sample_timepoint) %>%
  pivot_wider(names_from = sample_alias, values_from = log2_tpm)  %>%
  column_to_rownames("hugo_symbol") %>%
  as.matrix()

usethis::use_data(mbcp_log2_tpms_class, overwrite = TRUE, compress = "xz")
