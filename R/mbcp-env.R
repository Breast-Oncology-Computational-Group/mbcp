pkg_env <- new.env(parent = emptyenv())

pkg_env$set_names <- c("receptor", "pam50", "dna_alts", "mutation", "histology",
                       "treatment", "metastasis", "categories", "germline_alts", "dna_signatures")
