library(vroom)
library(dplyr)
library(janitor)
library(tidyr)
library(ggplot2)
library(ggpubr)

germ_hits <- vroom("data-raw/MBCproject_germline_charger_hits.all.merged.tsv",
                       na = c("", "NA", "__UNKNOWN__", "__UNKNOWN__,wex,wex", "."),
                       guess_max = 1000) %>%
  clean_names()
germ_pass <- vroom("data-raw/germline_variants_check_AF0005_FINAL.csv",
                   na = c("", "NA", "__UNKNOWN__", "__UNKNOWN__,wex,wex", ".")) %>%
  clean_names() %>%
  select(-x1)
tcga_genes <- vroom("data-raw/tcga_cancer_predisposition_genes.tsv") %>%
  clean_names()
cn_hits <- vroom("data-raw/agg_germline.ID.chrs8.CN_JGTZ.tsv",
                 na = c("", "NA", "__UNKNOWN__", "__UNKNOWN__,wex,wex", ".")) %>%
  clean_names()
maf <- vroom("data-raw/MBCproject_379_final_pairset.aggregated.corrected.wOncoKb.maf",
             quote = "", delim = "\t",
             na = c("", "NA", "__UNKNOWN__", "__UNKNOWN__,wex,wex", "."),
             guess_max = 1000) %>%
  clean_names()

sample_patient <- vroom("data-raw/sample_patient_metadata_HRHER2status_PAM50_MutSignatures_GenesInterest_wRNA.csv") %>%
  clean_names()  %>%
  filter(!is.na(wes_sample_id)) %>%
  select(patient_id, wes_sample_id, clean_participant) %>%
  mutate(patient_id = tolower(patient_id),
         wes_sample_id = tolower(wes_sample_id))
##create sub df of maf and annotate final states to plot----
#For Variant classification that has a functional annotation by oncoKB - override that and mark it on the coMut override order 1. Mutational effect 2. Hotspot 3. Multi hits
#For eg : any gene that has >1 non silent mutation in a sample will  be marked as multhit (eg ESR1 in 22 and ERBB2 in 32).
maf <- maf %>%
  mutate(final_variant_classification = ifelse(!mutation_effect %in% c("Unknown","Inconclusive"), mutation_effect,
                                               ifelse(!is.na(is_a_hotspot) | !is.na(is_a_3d_hotspot), paste('Hotspot'),
                                                      variant_classification))) %>%
  mutate(final_variant_classification = gsub("Likely ","",final_variant_classification)) %>%
  mutate(final_variant_classification = gsub("-of-function","_of_function",final_variant_classification))

#The "known" neutrals are classified as their variant type. This is so we are consistent with how we classify non-knonw missense muts. Otherwise these would be classified as "Other mutation"
maf <- maf %>%
  mutate(final_variant_classification = if_else(final_variant_classification=="Neutral", variant_classification, final_variant_classification))

#To distinguish this from the known LOF, we use a Putative LOF classification
putative_lof_list = c("Splice_Site", "Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation",
                      "Nonstop_Mutation","In_Frame_Del","Start_Codon_SNP","In_Frame_Ins","Start_Codon_Ins",
                      "Stop_Codon_Del","Stop_Codon_Ins")

maf <- maf %>% mutate(final_variant_classification = if_else(final_variant_classification %in%
                                                             putative_lof_list,"Putative_loss_of_function",final_variant_classification))

silent_remove <- c("Silent","Intron","3'UTR","5'Flank","5'UTR","IGR","lincRNA","RNA")

maf <- maf %>%
  filter(!(final_variant_classification %in% silent_remove)) %>%
  mutate(sample_id = tolower(sample_id)) %>%
  select(hugo_symbol, variant_type, sample_id, variant_classification = final_variant_classification)

##  Setting VUS SNV to anything that's not Loss_of_function because we're focusing on germline LOF mutation.
maf <- maf %>%
  mutate(variant_classification = if_else(variant_classification %in% c("Putative_loss_of_function", "Loss_of_function"),
                                          variant_classification, "VUS SNV"))

maf <- maf %>%
  left_join(sample_patient, by = c("sample_id" = "wes_sample_id")) %>%
  select(-sample_id) %>%
  distinct()

## Get only non artifact germline alts
germ_pass <- germ_pass %>%
  filter(artifact %in% c("N", "N/M")) %>%
  unite("unique_name", chr:alt, remove = FALSE)

germ_hits <- germ_hits %>%
  mutate(af_popmax = as.numeric(af_popmax)) %>%
  unite("unique_name", chr:alt, remove = FALSE) %>%
  filter(charger_score > 8, af_popmax < 0.0005| is.na(af_popmax),
         gene_ref_gene %in% tcga_genes$gene_symbol,
         unique_name %in% germ_pass$unique_name) %>%
  mutate(patient_id = tolower(patient_id)) %>%
  select(unique_name, hugo_symbol = gene_ref_gene, patient_id, func_ref_gene, exonic_func_ref_gene, charger_score) %>%
  distinct()

cn_hits <- cn_hits %>%
  mutate(af_popmax = as.numeric(af_popmax)) %>%
  unite("unique_name", chr:alt, remove = FALSE) %>%
  filter(charger_score > 8, af_popmax < 0.0005 | is.na(af_popmax),
         gene_ref_gene %in% tcga_genes$gene_symbol,
         unique_name %in% germ_pass$unique_name) %>%
  mutate(patient_id = tolower(patient_id)) %>%
  select(unique_name, hugo_symbol = gene_ref_gene, patient_id, func_ref_gene, exonic_func_ref_gene, charger_score) %>%
  mutate(variant_classification = "LOH") %>%
  distinct()

## Join with maf to get the final variant classification.
## If a patient has multiple germ_hits in the same gene, we want to keep them all,
## but we don't want to join it to multiple entries in the maf (multiple alts in the same gene)
germ_hits_maf <- germ_hits %>%
  left_join(maf, by = c("hugo_symbol", "patient_id"), multiple = "first") %>%
  filter(!is.na(variant_classification)) %>%
  select(unique_name, hugo_symbol, patient_id, func_ref_gene, exonic_func_ref_gene, charger_score, variant_classification)

mbcp_germline_hits <- germ_hits %>%
  anti_join(cn_hits, by = c("unique_name", "patient_id", "exonic_func_ref_gene")) %>%
  anti_join(germ_hits_maf, by = c("unique_name", "patient_id", "exonic_func_ref_gene")) %>%
  bind_rows(cn_hits, germ_hits_maf) %>%
  left_join(sample_patient, by = "patient_id", relationship = "many-to-many") %>%
  select(-wes_sample_id)

mbcp_germline_hits <- mbcp_germline_hits %>%
  mutate(somatic_event = ifelse(!patient_id %in% sample_patient$patient_id, "No Tumor WES",
                                ifelse(is.na(variant_classification), "None Detected",
                                       variant_classification))) %>%
  select(-variant_classification)

usethis::use_data(mbcp_germline_hits, overwrite = TRUE, compress = "xz")
