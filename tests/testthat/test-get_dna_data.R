test_that("get_mbcp_mutations returns a data frame", {
  df <- get_mbcp_mutations()
  expect_s3_class(df, "data.frame")
})

test_that("get_mbcp_mutations throws error if samples is not a character vector", {
  expect_error(get_mbcp_mutations(samples = 1), "samples should be a character vector")
  expect_error(get_mbcp_mutations(samples = c(TRUE, FALSE)), "samples should be a character vector")
  expect_error(get_mbcp_mutations(samples = NA), "samples should be a character vector")
})

test_that("get_mbcp_mutations returns data filtered by samples", {
  ss <- sample(mbcp_mutations$Sample_ID, 5)
  ss_not <- setdiff(mbcp_mutations$Sample_ID, ss)
  df <- get_mbcp_mutations(samples = ss)
  expect_true(all(ss %in% df$Sample_ID))
  expect_false(all(ss_not %in% df$Sample_ID))
})

test_that("get_mbcp_mutations throws warning if sample is missing", {
  ss <- c(sample(mbcp_mutations$Sample_ID, 5), "invalid sample")
  expect_warning(get_mbcp_mutations(samples = ss), "some samples were not found in the mbcp_mutations dataset")
})

test_that("get_mbcp_mutations still returns filtered data if sample is missing", {
  ss <- c( "invalid sample", sample(mbcp_mutations$Sample_ID, 5))
  ss_not <- setdiff(mbcp_mutations$Sample_ID, ss)
  df <- suppressWarnings(get_mbcp_mutations(samples = ss))
  expect_true(all(ss[-1] %in% df$Sample_ID))
  expect_false(all(ss_not %in% df$Sample_ID))
})

test_that("get_mbcp_mutations returns entire dataset when samples = NULL", {
  df <- get_mbcp_mutations()
  expect_equal(nrow(mbcp_mutations), nrow(df))
})

test_that("get_mbcp_data_alts throws error if genes is not a character vector", {
  expect_error(get_mbcp_mutations(genes = 1), "genes should be a character vector")
  expect_error(get_mbcp_mutations(genes = c(TRUE, FALSE)), "genes should be a character vector")
  expect_error(get_mbcp_mutations(genes = NA), "genes should be a character vector")
})

test_that("get_mbcp_mutations returns data filtered by genes", {
  gs <- sample(mbcp_mutations$Hugo_Symbol, 100)
  df <- get_mbcp_mutations(genes = gs)
  gs_not <- setdiff(unique(mbcp_mutations$Hugo_Symbol), gs)
  expect_true(all(gs %in% df$Hugo_Symbol))
  expect_false(all(gs_not %in% df$Hugo_Symbol))
})

# It's common for genes to not have alts
test_that("get_mbcp_mutations does not throw warning with invalid genes", {
  gs <- c(sample(mbcp_mutations$Hugo_Symbol, 5), "invalid gene")
  expect_no_warning(get_mbcp_mutations(genes = gs))
})

test_that("get_mbcp_mutations throws error if min_ccf is not a single numeric value", {
  expect_error(get_mbcp_mutations(min_ccf = NA), "min_ccf should be a numeric value")
  expect_error(get_mbcp_mutations(min_ccf = FALSE), "min_ccf should be a numeric value")
  expect_error(get_mbcp_mutations(min_ccf = "character"), "min_ccf should be a numeric value")
  expect_error(get_mbcp_mutations(min_ccf = c(0.5, 4.8)), "min_ccf should be a numeric value")
})

test_that("get_mbcp_mutations filters alts by ccf_hat", {
  min_ccf <- runif(1)
  df <- get_mbcp_mutations(min_ccf = min_ccf)
  expect_true(all(df$ccf_hat >= min_ccf))
})

test_that("get_mbcp_cnvs returns a data frame", {
  df <- get_mbcp_cnvs()
  expect_s3_class(df, "data.frame")
})

test_that("get_mbcp_cnvs throws error if samples is not a character vector", {
  expect_error(get_mbcp_cnvs(samples = 1), "samples should be a character vector")
  expect_error(get_mbcp_cnvs(samples = c(TRUE, FALSE)), "samples should be a character vector")
  expect_error(get_mbcp_cnvs(samples = NA), "samples should be a character vector")
})

test_that("get_mbcp_cnvs returns data filtered by samples", {
  ss <- sample(mbcp_cnvs$sample_id, 5)
  ss_not <- setdiff(mbcp_cnvs$sample_id, ss)
  df <- get_mbcp_cnvs(samples = ss)
  expect_true(all(ss %in% df$sample_id))
  expect_false(all(ss_not %in% df$sample_id))
})

test_that("get_mbcp_cnvs throws warning if sample is missing", {
  ss <- c(sample(mbcp_cnvs$sample_id, 5), "invalid sample")
  expect_warning(get_mbcp_cnvs(samples = ss), "some samples were not found in the mbcp_cnvs dataset")
})

test_that("get_mbcp_cnvs still returns filtered data if sample is missing", {
  ss <- c( "invalid sample", sample(mbcp_cnvs$sample_id, 5))
  ss_not <- setdiff(mbcp_cnvs$sample_id, ss)
  df <- suppressWarnings(get_mbcp_cnvs(samples = ss))
  expect_true(all(ss[-1] %in% df$sample_id))
  expect_false(all(ss_not %in% df$sample_id))
})

test_that("get_mbcp_cnvs returns entire dataset when samples = NULL", {
  df <- get_mbcp_cnvs()
  expect_equal(nrow(mbcp_cnvs), nrow(df))
})

test_that("get_mbcp_data_alts throws error if genes is not a character vector", {
  expect_error(get_mbcp_cnvs(genes = 1), "genes should be a character vector")
  expect_error(get_mbcp_cnvs(genes = c(TRUE, FALSE)), "genes should be a character vector")
  expect_error(get_mbcp_cnvs(genes = NA), "genes should be a character vector")
})

test_that("get_mbcp_cnvs returns data filtered by genes", {
  gs <- sample(mbcp_cnvs$hugo_symbol, 100)
  df <- get_mbcp_cnvs(genes = gs)
  gs_not <- setdiff(unique(mbcp_cnvs$hugo_symbol), gs)
  expect_true(all(gs %in% df$hugo_symbol))
  expect_false(all(gs_not %in% df$hugo_symbol))
})

# It's common for genes to not have alts
test_that("get_mbcp_cnvs does not throw warning with invalid genes", {
  gs <- c(sample(mbcp_cnvs$hugo_symbol, 5), "invalid gene")
  expect_no_warning(get_mbcp_cnvs(genes = gs))
})

test_that("get_mbcp_cnvs throws error if min_ccf is not a single numeric value", {
  expect_error(get_mbcp_cnvs(min_ccf = NA), "min_ccf should be a numeric value")
  expect_error(get_mbcp_cnvs(min_ccf = FALSE), "min_ccf should be a numeric value")
  expect_error(get_mbcp_cnvs(min_ccf = "character"), "min_ccf should be a numeric value")
  expect_error(get_mbcp_cnvs(min_ccf = c(0.5, 4.8)), "min_ccf should be a numeric value")
})

test_that("get_mbcp_cnvs filters alts by ccf_hat", {
  min_ccf <- runif(1)
  df <- get_mbcp_cnvs(min_ccf = min_ccf)
  expect_true(all(df$cancer_cell_frac_a2 > min_ccf | df$cancer_cell_frac_a1 > min_ccf))
})


###summary
## variant_type = "MUT", cnaps = NA
## no NA in other columns
## unique hugo_symbol, sample_id and final_variant_classification
## multihits -> "," in variant_classification
## n_alts > 1 -> "," in variant_classification
## proteinchange = NA for variant_type = "CNV"

test_that("add_any_amp throws error if .data does not have a variant_classification column", {

  expect_error(mocked_dna_alts %>%
    dplyr::select(-variant_classification) %>%
    add_any_amp(), "'variant_classification' column missing in data")
})

test_that("add_any_amp returns the correct number of rows", {

  any_rows <- mocked_dna_alts %>%
    dplyr::filter(variant_classification %in% c("AMP", "HighAMP", "FocalHighAMP")) %>%
    nrow()

  expect_equal(nrow(add_any_amp(mocked_dna_alts)), any_rows + nrow(mocked_dna_alts))

})
