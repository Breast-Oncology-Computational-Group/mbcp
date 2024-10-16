library(vroom)
library(dplyr)
library(janitor)
### This file must be previously annotated with variant types
seg_cn <- vroom("data-raw//MBCproject_379.seg.cn.metadata.annot.mod.ForCoMut_JGTZ.tsv",
               quote = "", delim = "\t",
               na = c("", "NA", "__UNKNOWN__", "__UNKNOWN__,wex,wex", "."),
               guess_max = 1000) %>%
  clean_names()

### There are no ABSOLUTE annotations, so we take ccf_hat values for both alleles and use the higher one
seg_cn <- seg_cn %>%
  select(-new_ploidy, -new_purity, -participant, -entity_pair_id, -control_sample, -n, -key) %>%
  select(sample_id = sample, hugo_symbol = gene, everything())

### Rename variant classification.
### bi.allelic.inactivation to LOF_and_LOH
### possible.bi.allelic.inactivation and likely.bi.allelic.inactivation are renamed but later removed
seg_cn <- seg_cn %>% mutate(variant_classification =
                              gsub("likely.bi.allelic.inactivation", "Likely_biallelic_inactivation", variant_classification))
seg_cn <- seg_cn %>% mutate(variant_classification =
                              gsub("possible.bi.allelic.inactivation", "Possible_biallelic_inactivation", variant_classification))
seg_cn <- seg_cn %>% mutate(variant_classification =
                              gsub("bi.allelic.inactivation", "LOF_and_LOH", variant_classification))
seg_cn <- seg_cn %>% mutate(variant_classification =
                              gsub("deep.deletion","DeepDEL", variant_classification))

## Joining alts in MAF with alts in seg, and filtering out GAIN, Possible_biallelic_inactivation, Likely_biallelic_inactivation
seg_cn <-seg_cn %>%
  filter(!variant_classification %in% c("Possible_biallelic_inactivation","Likely_biallelic_inactivation"))

## Filtering out GAIN, Possible_biallelic_inactivation, Likely_biallelic_inactivation
## when they're listed with something else
seg_cn <- seg_cn %>% mutate(variant_classification =
                             gsub(",Possible_biallelic_inactivation|Possible_biallelic_inactivation,|,Likely_biallelic_inactivation|Likely_biallelic_inactivation,", "",
                                            variant_classification))

## "De_novo_Start_InFrame"  and   "De_novo_Start_OutOfFrame" to Other_Mutation
## LOF_and_LOH,DeepDEL to DeepDEL
seg_cn <- seg_cn %>%
  mutate(variant_classification = ifelse(variant_classification == "LOF_and_LOH,DeepDEL", "DeepDEL", variant_classification))

seg_cn <- seg_cn %>%
  mutate(sample_id = tolower(sample_id))

##############################################
##### Remove DeepDEL in BLOOD samples
load("data/mbcp_clinical_data.rda")

bbs_samples = tolower(mbcp_clinical_data %>% filter(bx_location=="BLOOD") %>% pull(wes_sample_id))

mbcp_cnvs <- seg_cn %>%
  filter(!((sample_id %in% bbs_samples) & (variant_classification =="DeepDEL")))

usethis::use_data(mbcp_cnvs, overwrite = TRUE, compress = "xz")
