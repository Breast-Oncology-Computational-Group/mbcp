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

#### get_mbcp_dna_alts_summary tests
test_that("get_mbcp_mutations throws error if samples is not a character vector", {
  expect_error(get_mbcp_dna_alts_summary(samples = 1), "samples should be a character vector")
  expect_error(get_mbcp_dna_alts_summary(samples = c(TRUE, FALSE)), "samples should be a character vector")
  expect_error(get_mbcp_dna_alts_summary(samples = NA), "samples should be a character vector")
})

test_that("get_mbcp_dna_alts_summary throws error if genes is not a character vector", {
  expect_error(get_mbcp_dna_alts_summary(genes = 1), "genes should be a character vector")
  expect_error(get_mbcp_dna_alts_summary(genes = c(TRUE, FALSE)), "genes should be a character vector")
  expect_error(get_mbcp_dna_alts_summary(genes = NA), "genes should be a character vector")
})

test_that("get_mbcp_dna_alts_summary throws error if variant_type is not valid", {
  expect_error(get_mbcp_dna_alts_summary(variant_type = NA), "variant_type should be one of MUT, CNV, both")
  expect_error(get_mbcp_dna_alts_summary(variant_type = FALSE), "variant_type should be one of MUT, CNV, both")
  expect_error(get_mbcp_dna_alts_summary(variant_type = "character"), "variant_type should be one of MUT, CNV, both")
  expect_error(get_mbcp_dna_alts_summary(variant_type = c(0.5, 4.8)), "variant_type should be one of MUT, CNV, both")
})

test_that("get_mbcp_dna_alts_summary calls get_mbcp_mutations and get_mbcp_cnvs with correct arguments", {
  mmut <- mock(mocked_mutations)
  mcnv <- mock(mocked_cnvs)
  ss <- unique(sample(mocked_mutations$Sample_ID, 25))
  gg <- unique(sample(mocked_mutations$Hugo_Symbol, 25))
  ccf <-runif(1)
  with_mocked_bindings(code = {
    get_mbcp_dna_alts_summary(samples = ss, genes = gg, min_ccf = ccf)
    margs <- mock_args(mmut)
    expect_called(mmut, 1)
    expect_equal(margs[[1]], list(samples = ss, genes = gg, min_ccf = ccf))
    margs <- mock_args(mcnv)
    expect_called(mcnv, 1)
    expect_equal(margs[[1]], list(samples = ss, genes = gg))
  }, get_mbcp_cnvs = mcnv, get_mbcp_mutations = mmut)
})

test_that("get_mbcp_dna_alts_summary returns correct Multi-hit entries", {
  mmut <- mock(mocked_mutations)
  mcnv <- mock(mocked_cnvs)
  with_mocked_bindings(code = {
    df <- get_mbcp_dna_alts_summary()
    df <- df %>%
      mutate(ok = variant_classification == "Multi_hit" & n_alts > 1 |
               variant_classification != "Multi_hit" & n_alts == 1 |
               variant_type == "CNV")
    expect_true(all(df$ok))
  }, get_mbcp_cnvs = mcnv, get_mbcp_mutations = mmut)
})

test_that("get_mbcp_dna_alts_summary returns NA cnaps for MUT", {
  all_alts <- get_mbcp_dna_alts_summary()
  all_alts <- all_alts %>%
    filter(variant_type == "MUT")
  expect_true(all(is.na(all_alts$all_cnaps)))
  expect_true(all(is.na(all_alts$max_cnap)))
})

test_that("get_mbcp_dna_alts_summary returns NA ccf for CNVS", {
  all_alts <- get_mbcp_dna_alts_summary()
  all_alts <- all_alts %>%
    filter(variant_type == "CNV")
  expect_true(all(is.na(all_alts$ccf_hat)))
  expect_true(all(is.na(all_alts$all_ccfs)))
  expect_true(all(is.na(all_alts$all_hgvs_protein_change)))
  expect_true(all(is.na(all_alts$all_hgvs_genomic_change)))
  expect_true(all(is.na(all_alts$all_protein_change)))
  expect_true(all(is.na(all_alts$all_genomic_change)))
  expect_false(all(is.na(all_alts$max_cnap)))
  expect_false(all(is.na(all_alts$all_cnaps)))
})

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
