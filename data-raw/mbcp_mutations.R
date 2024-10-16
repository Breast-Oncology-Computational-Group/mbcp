library(vroom)
library(dplyr)
library(janitor)
library(purrr)

## No janitor::clean_names to keep official maf column names
maf <- vroom("data-raw/MBCproject_379_final_pairset.aggregated.corrected.wOncoKb.maf",
             quote = "", delim = "\t",
             na = c("", "NA", "__UNKNOWN__", "__UNKNOWN__,wex,wex", "."),
             guess_max = 1000)

##############################################
##### Recover ccf_hat values from ABSOLUTE's distribution output
## we only need the na ccf_hats
maf_na_ccf <- maf %>%
  filter(is.na(ccf_hat))

# ccf_hat = NA when q_hat == 0, hs_q_hat_1 = 0 and hs_q_hat_2 = 0  .
# According to ABSOLUTE docs:
# q_hat: The total copy number of the segment on which the mutation occurs (modal.a1 + modal.a2)
# fit_SSNV_model.R line 160

# Related:
# HS_q_hat_1 & HS_q_hat_2: The allele specific copy numbers of the segment on which the mutation occurs.
# HSCN_subclonal_SCNA_fit_SSNVs.R line 27 and reconcile_clonal_homdels_with_obs_SSNVs function in
# get_SCNA_cancer_cell_fractions.R line 234
# So absolute thinks the mutation is in a region with 0 copy number

# we can't do anything with the ones that have no ccf distribution
# or the ones with really low counts
maf_na_ccf <- maf_na_ccf %>% filter(!is.na(X0), alt >= 2)

## the distribution is orderded in the maf, but just to make sure that
## we're not relying on starts_with to return the array in order
## and also cause adding the x1 is an extra step
dist_cols <- paste0("X", as.character(seq(from=0, to=1, by = 0.01)))

maf_na_ccf <- maf_na_ccf %>%
  mutate(distribution = pmap(pick(all_of(dist_cols)), c),
         imputed_ccf = (map_dbl(distribution, which.max)-1)/100) %>%
  select(Hugo_Symbol, Start_Position, End_Position, Chromosome, Variant_Classification,
         Variant_Type, Reference_Allele, Tumor_Seq_Allele1, Tumor_Seq_Allele2,  Protein_Change,
         Genome_Change, Sample_ID, -distribution, imputed_ccf)

################################################
##### Get final variant classification
# For Variant classification that has a functional annotation by oncoKB - override that and mark it on the coMut override order
# 1. Mutational effect 2. Hotspot 3. Multi hits
#For eg : any gene that has >1 non silent mutation in a sample will  be marked as multhit (eg ESR1 in 22 and ERBB2 in 32).
maf <- maf %>%
  mutate(final_variant_classification = ifelse(!MUTATION_EFFECT %in% c("Unknown","Inconclusive"), MUTATION_EFFECT,
                                               ifelse(!is.na(`IS-A-HOTSPOT`) | !is.na(`IS-A-HOTSPOT`), paste('Hotspot'),
                                                      Variant_Classification))) %>%
  mutate(final_variant_classification = gsub("Likely ","", final_variant_classification)) %>%
  mutate(final_variant_classification = gsub("-of-function","_of_function",final_variant_classification))

#The "known" neutrals are classified as their variant type. This is so we are consistent with how we
# classify non-knonw missense muts. Otherwise these would be classified as "Other mutation"
maf <-  maf %>%
  mutate(final_variant_classification = if_else(final_variant_classification == "Neutral",
                                                Variant_Classification, final_variant_classification))

#To distinguish this from the known LOF, we use a Putative LOF classification
putative_lof_list = c("Splice_Site", "Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation",
                    "Nonstop_Mutation","In_Frame_Del","Start_Codon_SNP","In_Frame_Ins","Start_Codon_Ins",
                    "Stop_Codon_Del","Stop_Codon_Ins")
maf <- maf %>%
  mutate(final_variant_classification = if_else(final_variant_classification %in%
                                                  putative_lof_list, "Putative_loss_of_function", final_variant_classification))

## Remove silent alts
silent_remove = c("Silent","Intron","3'UTR","5'Flank","5'UTR","IGR","lincRNA","RNA")

maf <- maf %>%
  filter(!(final_variant_classification %in% silent_remove)) %>%
  select(-Variant_Classification) %>%
  rename(Variant_Classification = final_variant_classification)

##############################################
##### Merge imputed ccf_hat values
maf <- maf %>%
  left_join(maf_na_ccf %>%
              select(Hugo_Symbol, Sample_ID, Genome_Change, imputed_ccf),
            by = c("Sample_ID", "Hugo_Symbol", "Genome_Change")) %>%
  mutate(ccf_hat = if_else(is.na(ccf_hat), imputed_ccf, ccf_hat)) %>%
  select(-imputed_ccf)

### define muts type----
match_list <- c("Gain_of_function", "Hotspot",
                "Loss_of_function", "Missense_Mutation","Putative_loss_of_function",
                "Switch_of_function", "Multi_hit")

maf <- maf %>%
  mutate(Variant_Classification = ifelse(Variant_Classification %in% match_list, Variant_Classification, "Other_Mutation"))

start_abs <- which(colnames(maf) == "A1.ix")
end_abs <- which(colnames(maf) == "X1")

start_onco <- which(colnames(maf) == "ONCOTREE_CODE")
end_onco <- which(colnames(maf) == "PX_CITATIONS")

mbcp_mutations <-  maf %>%
  mutate(Sample_ID = tolower(Sample_ID)) %>%
  select(Hugo_Symbol, Entrez_Gene_Id, NCBI_Build, Chromosome, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele1, Tumor_Seq_Allele2,
         Variant_Classification, Variant_Type, Tumor_Sample_Barcode = tumor_sample_barcode, Sample_ID, Genome_Change, Protein_Change,
         ref, alt, n_ref_count, n_alt_count, t_depth, all_of(start_abs:end_abs), all_of(start_onco:end_onco))


usethis::use_data(mbcp_mutations, overwrite = TRUE, compress = "xz")
