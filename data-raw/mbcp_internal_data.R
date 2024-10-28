library(vroom)
library(dplyr)

ss <- vroom("data-raw/signatures_categories_for_expr.csv")
gs <- vroom("data-raw/genes_categories_for_expr.csv")

mbcp_genes_signatures <- ss %>%
  inner_join(gs, by = "category", relationship = "many-to-many") %>%
  rename(signature_id = signature)

### CCNE1 is not in our expression matrix
mbcp_genes_signatures <- mbcp_genes_signatures %>%
  filter(hugo_symbol != "CCNE1" )

mbcp_genes <- vroom("data-raw/CDK46_Seth_gene_list_16177.ordered_3.plus_genes_from_Daniel.txt", delim = "\t") %>% pull(Hugo_Symbol)

samples_to_remove <- vroom("data-raw/samples_to_remove_from_dict.tsv", delim = "\t")
usethis::use_data(mbcp_genes_signatures, mbcp_genes, samples_to_remove,
                  overwrite = TRUE, internal = TRUE, compress = "xz")

