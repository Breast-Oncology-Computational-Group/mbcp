mocked_dna_alts <- data.frame(
  hugo_symbol = sample(paste0("g", 1:10), 80, replace = T),
  sample_id = sample(paste0("s", 20:40), 80, replace = T),
  ccf_hat = runif(80),
  variant_type = sample(c("MUT", "CNV"), 80, replace = T),
  variant_classification = sample(c("AMP", "HighAMP", "FocalHighAMP", "DeepDEL"), 80, replace = T )
)


test_that("get_mbcp_dna_alts returns a data frame", {
  df <- get_mbcp_dna_alts()
  expect_s3_class(df, "data.frame")
})

test_that("get_mbcp_dna_alts throws error if samples is not a character vector", {
  expect_error(get_mbcp_dna_alts(samples = 1), "samples should be a character vector")
  expect_error(get_mbcp_dna_alts(samples = c(TRUE, FALSE)), "samples should be a character vector")
  expect_error(get_mbcp_dna_alts(samples = NA), "samples should be a character vector")
})

test_that("get_mbcp_dna_alts returns data filtered by samples", {
  ss <- sample(mbcp_dna_alts$sample_id, 5)
  ss_not <- setdiff(mbcp_dna_alts$sample_id, ss)
  df <- get_mbcp_dna_alts(samples = ss)
  expect_true(all(ss %in% df$sample_id))
  expect_false(all(ss_not %in% df$sample_id))
})

test_that("get_mbcp_dna_alts throws warning if sample is missing", {
  ss <- c(sample(mbcp_dna_alts$sample_id, 5), "invalid sample")
  expect_warning(get_mbcp_dna_alts(samples = ss), "some samples were not found in the mbcp_dna_alts dataset")
})

test_that("get_mbcp_dna_alts still returns filtered data if sample is missing", {
  ss <- c( "invalid sample", sample(mbcp_dna_alts$sample_id, 5))
  ss_not <- setdiff(mbcp_dna_alts$sample_id, ss)
  df <- suppressWarnings(get_mbcp_dna_alts(samples = ss))
  expect_true(all(ss[-1] %in% df$sample_id))
  expect_false(all(ss_not %in% df$sample_id))
})

test_that("get_mbcp_dna_alts returns entire dataset when samples = NULL", {
  df <- get_mbcp_dna_alts()
  expect_equal(nrow(mbcp_dna_alts), nrow(df))
})

test_that("get_mbcp_data_alts throws error if genes is not a character vector", {
  expect_error(get_mbcp_dna_alts(genes = 1), "genes should be a character vector")
  expect_error(get_mbcp_dna_alts(genes = c(TRUE, FALSE)), "genes should be a character vector")
  expect_error(get_mbcp_dna_alts(genes = NA), "genes should be a character vector")
})

test_that("get_mbcp_dna_alts returns data filtered by genes", {
  gs <- sample(mbcp_dna_alts$hugo_symbol, 100)
  df <- get_mbcp_dna_alts(genes = gs)
  gs_not <- setdiff(unique(mbcp_dna_alts$hugo_symbol), gs)
  expect_true(all(gs %in% df$hugo_symbol))
  expect_false(all(gs_not %in% df$hugo_symbol))
})

# It's common for genes to not have alts
test_that("get_mbcp_dna_alts does not throw warning with invalid genes", {
  gs <- c(sample(mbcp_dna_alts$hugo_symbol, 5), "invalid gene")
  expect_no_warning(get_mbcp_dna_alts(genes = gs))
})

test_that("get_mbcp_dna_alts throws error if min_ccf is not a single numeric value", {
  expect_error(get_mbcp_dna_alts(min_ccf = NA), "min_ccf should be a numeric value")
  expect_error(get_mbcp_dna_alts(min_ccf = FALSE), "min_ccf should be a numeric value")
  expect_error(get_mbcp_dna_alts(min_ccf = "character"), "min_ccf should be a numeric value")
  expect_error(get_mbcp_dna_alts(min_ccf = c(0.5, 4.8)), "min_ccf should be a numeric value")
})

test_that("get_mbcp_dna_alts filters alts by ccf_hat", {
  min_ccf <- runif(1)
  df <- get_mbcp_dna_alts(min_ccf = min_ccf)
  expect_true(all(df$ccf_hat >= min_ccf))
})

test_that("get_mbcp_snvs calls get_mbcp_dna_alts", {
  dna_alts <- mock(mocked_dna_alts)

  with_mocked_bindings(code = {
    get_mbcp_snvs()
    expect_called(dna_alts, 1)
  }, get_mbcp_dna_alts = dna_alts)
})

test_that("get_mbcp_snvs calls get_mbcp_dna_alts with correct arguments", {
  dna_alts <- mock(mocked_dna_alts)
  ss <- sample(mocked_dna_alts$sample_id, 5)
  gs <- sample(mocked_dna_alts$hugo_symbol, 5)
  min_ccf <- runif(1)

  with_mocked_bindings(code = {
    get_mbcp_snvs(samples = ss, genes = gs, min_ccf = min_ccf)
    margs <- mock_args(dna_alts)
    expect_equal(margs[[1]], list(samples = ss, genes = gs, min_ccf = min_ccf))
    expect_called(dna_alts, 1)
  }, get_mbcp_dna_alts = dna_alts)
})

test_that("get_mbcp_snvs returns entries with variant_type = MUT", {
  dna_alts <- mock(mocked_dna_alts)
  with_mocked_bindings(code = {
    df <- get_mbcp_snvs()
    expect_true(all(df$variant_type == "MUT"))
  }, get_mbcp_dna_alts = dna_alts)
})

test_that("get_mbcp_cnvs calls get_mbcp_dna_alts", {
  dna_alts <- mock(mocked_dna_alts)

  with_mocked_bindings(code = {
    get_mbcp_cnvs()
    expect_called(dna_alts, 1)
  }, get_mbcp_dna_alts = dna_alts)
})

test_that("get_mbcp_cnvs calls get_mbcp_dna_alts with correct arguments", {
  dna_alts <- mock(mocked_dna_alts)
  ss <- sample(mocked_dna_alts$sample_id, 5)
  gs <- sample(mocked_dna_alts$hugo_symbol, 5)
  min_ccf <- runif(1)

  with_mocked_bindings(code = {
    get_mbcp_cnvs(samples = ss, genes = gs, min_ccf = min_ccf)
    margs <- mock_args(dna_alts)
    expect_equal(margs[[1]], list(samples = ss, genes = gs, min_ccf = min_ccf))
    expect_called(dna_alts, 1)
  }, get_mbcp_dna_alts = dna_alts)
})

test_that("get_mbcp_cnvs returns entries with variant_type = MUT", {
  dna_alts <- mock(mocked_dna_alts)
  with_mocked_bindings(code = {
    df <-  get_mbcp_cnvs()
    expect_true(all(df$variant_type == "CNV"))
  }, get_mbcp_dna_alts = dna_alts)
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
