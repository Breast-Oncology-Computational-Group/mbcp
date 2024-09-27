library(vroom)
library(dplyr)
library(janitor)
library(purrr)

maf <- vroom("data-raw/MBCproject_379_final_pairset.aggregated.corrected.wOncoKb.maf",
             quote = "", delim = "\t",
             na = c("", "NA", "__UNKNOWN__", "__UNKNOWN__,wex,wex", "."),
             guess_max = 1000) %>%
  clean_names()

##############################################
##### Recover ccf_hat values from ABSOLUTE's distribution output
## we only need the na ccf_hats
maf_na_ccf <- maf %>%
  filter(is.na(ccf_hat))

# we can't do anything with the ones that have no ccf distribution
# or the ones with really low counts
maf_na_ccf <- maf_na_ccf %>% filter(!is.na(x0), alt >= 2)

## the distribution is orderded in the maf, but just to make sure that
## we're not relying on starts_with to return the array in order
## and also cause adding the x1 is an extra step
dist_cols <- paste0("x", sub("\\.", "_", as.character(seq(from=0, to=1, by = 0.01))))

maf_na_ccf <- maf_na_ccf %>%
  mutate(distribution = pmap(pick(all_of(dist_cols)), c),
         imputed_ccf = (map_dbl(distribution, which.max)-1)/100) %>%
  select(hugo_symbol, start_position, end_position, chromosome, variant_classification,
         variant_type, reference_allele, tumor_seq_allele1, tumor_seq_allele2,  protein_change,
         hgvs_genomic_change, hgvs_protein_change, sample_id, -distribution, imputed_ccf)

################################################
##### Get final variant classification
# For Variant classification that has a functional annotation by oncoKB - override that and mark it on the coMut override order
# 1. Mutational effect 2. Hotspot 3. Multi hits
#For eg : any gene that has >1 non silent mutation in a sample will  be marked as multhit (eg ESR1 in 22 and ERBB2 in 32).
maf <- maf %>%
  mutate(final_variant_classification = ifelse(!mutation_effect %in% c("Unknown","Inconclusive"), mutation_effect,
                                               ifelse(!is.na(is_a_hotspot) | !is.na(is_a_3d_hotspot), paste('Hotspot'),
                                                      variant_classification))) %>%
  mutate(final_variant_classification = gsub("Likely ","",final_variant_classification)) %>%
  mutate(final_variant_classification = gsub("-of-function","_of_function",final_variant_classification))

#The "known" neutrals are classified as their variant type. This is so we are consistent with how we
# classify non-knonw missense muts. Otherwise these would be classified as "Other mutation"
maf <-  maf %>%
  mutate(final_variant_classification = if_else(final_variant_classification=="Neutral", variant_classification, final_variant_classification))

#To distinguish this from the known LOF, we use a Putative LOF classification
putative_lof_list = c("Splice_Site", "Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation",
                    "Nonstop_Mutation","In_Frame_Del","Start_Codon_SNP","In_Frame_Ins","Start_Codon_Ins",
                    "Stop_Codon_Del","Stop_Codon_Ins")
maf <- maf %>%
  mutate(final_variant_classification = if_else(final_variant_classification %in%
                                                  putative_lof_list, "Putative_loss_of_function", final_variant_classification))

## Remove silent alts
silent_remove = c("Silent","Intron","3'UTR","5'Flank","5'UTR","IGR","lincRNA","RNA")
maf_sub <- maf %>%
  filter(!(final_variant_classification %in% silent_remove))

##############################################
##### Merge imputed ccf_hat values
maf_sub <- maf_sub %>%
  left_join(maf_na_ccf %>%
              select(hugo_symbol, sample_id, hgvs_genomic_change, imputed_ccf),
            by = c("sample_id", "hugo_symbol", "hgvs_genomic_change")) %>%
  mutate(ccf_hat = if_else(is.na(ccf_hat), imputed_ccf, ccf_hat)) %>%
  select(-imputed_ccf)

##############################################
##### Add CN alts
### This file must be previously annotated with variant types
seg_cn = vroom("data-raw//MBCproject_379.seg.cn.metadata.annot.mod.ForCoMut_JGTZ.tsv",
               quote = "", delim = "\t",
               na = c("", "NA", "__UNKNOWN__", "__UNKNOWN__,wex,wex", "."),
               guess_max = 1000) %>%
  clean_names()

### There are no ABSOLUTE annotations, so we take ccf_hat values for both alleles and use the higher one
seg_sub <- seg_cn %>%
  filter(!is.na(variant_classification)) %>%
  select(hugo_symbol = gene, variant_type, sample_id = sample, variant_classification, ccf_hat_a1 = ccf_ci95_high_a1, ccf_hat_a2 = ccf_ci95_high_a2) %>%
  mutate(ccf_hat = pmax(ccf_hat_a1, ccf_hat_a2, na.rm = T)) %>%
  select(-ccf_hat_a1,-ccf_hat_a2)

### Rename variant classification
seg_sub<-seg_sub %>% mutate(variant_classification=gsub("likely.bi.allelic.inactivation","Likely_biallelic_inactivation",variant_classification))
seg_sub<-seg_sub %>% mutate(variant_classification=gsub("possible.bi.allelic.inactivation","Possible_biallelic_inactivation",variant_classification))
seg_sub<-seg_sub %>% mutate(variant_classification=gsub("bi.allelic.inactivation","LOF_and_LOH",variant_classification))
seg_sub<-seg_sub %>% mutate(variant_classification=gsub("deep.deletion","DeepDEL",variant_classification))

## Joining alts in MAF with alts in seg, and filtering out GAIN, Possible_biallelic_inactivation, Likely_biallelic_inactivation
maf_seg_sub = bind_rows(maf_sub %>%
                          mutate(variant_type = "MUT") %>%
                          select(hugo_symbol, variant_type, sample_id = tumor_sample_uuid,
                                 variant_classification = final_variant_classification, ccf_hat, protein_change, genome_change),
                        seg_sub %>% filter(!variant_classification == "GAIN")) %>%
  filter(!variant_classification %in% c("Possible_biallelic_inactivation","Likely_biallelic_inactivation"))

## Filtering out GAIN, Possible_biallelic_inactivation, Likely_biallelic_inactivation
## when they're listed with something else
maf_seg_sub = maf_seg_sub %>% mutate(variant_classification =
                                       gsub("GAIN,|,Possible_biallelic_inactivation|Possible_biallelic_inactivation,|,Likely_biallelic_inactivation|Likely_biallelic_inactivation,", "",
                                            variant_classification))

maf_seg_sub = maf_seg_sub %>% mutate(variant_classification=if_else(variant_classification=="AMP,LOF_and_LOH","LOF_and_LOH",
                                                                    variant_classification))
maf_seg_sub = maf_seg_sub  %>% filter(!variant_classification == "GAIN") %>%
  filter(!variant_classification %in% c("Possible_biallelic_inactivation","Likely_biallelic_inactivation"))

### define muts type----
match_list <- c("AMP", "DeepDEL", "FocalHighAMP", "Gain_of_function", "HighAMP", "Hotspot", "LOF_and_LOH",
                "Loss_of_function", "Missense_Mutation","Putative_loss_of_function",
                "Switch_of_function", "Multi_hit",
                unique(maf_seg_sub$variant_classification[grepl(",", maf_seg_sub$variant_classification)]))

## "De_novo_Start_InFrame"  and   "De_novo_Start_OutOfFrame" to Other_Mutation
## LOF_and_LOH,DeepDEL to DeepDEL
maf_seg_sub <- maf_seg_sub %>%
  mutate(variant_classification = ifelse(variant_classification %in% match_list, variant_classification, "Other_Mutation"),
         variant_classification = ifelse(variant_classification == "LOF_and_LOH,DeepDEL", "DeepDEL", variant_classification))

maf_seg_sub <- maf_seg_sub %>%
  mutate(sample_id = tolower(sample_id))

##############################################
##### Remove DeepDEL in BLOOD samples
load("data/mbcp_clinical_data.rda")

bbs_samples = tolower(mbcp_clinical_data %>% filter(bx_location=="BLOOD") %>% pull(wes_sample_id))

maf_seg_sub = maf_seg_sub %>%
  filter(!((sample_id %in% bbs_samples) & (variant_type == "CNV") & (variant_classification =="DeepDEL")))

# Should we remove alts with ccf NA values?
mbcp_dna_alts <-  maf_seg_sub %>%
  filter(!is.na(ccf_hat))

usethis::use_data(mbcp_dna_alts, overwrite = TRUE, compress = "xz")
