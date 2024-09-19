mat <- matrix(c(10, 1, 2, 5, 30,
                15, 2, 3, 4, 15,
                20, 3, 4, 3, 20,
                25, 4, 5, 2, 25,
                30, 5, 6, 1, 10), nrow = 5, byrow = TRUE)
colnames(mat) <- paste0("s", 1:5)
rownames(mat) <- paste0("g", 1:5)

result <- matrix(c(0.4, 0.25,  0.4, 1.25, 1.2,
                   0.6, 0.50,  0.6, 1.00,  0.6,
                   0.8, 0.75,  0.8, 0.75,  0.8,
                   1.0, 1.00,  1.0, 0.50,  1.0,
                   1.2, 1.25,  1.2, 0.25,  0.4), nrow = 5, byrow = T)
colnames(result) <- paste0("s", 1:5)
rownames(result) <- paste0("g", 1:5)

test_that("get_mbcp_tpms returns a numeric matrix", {
  tpms <- get_mbcp_tpms()
  expect_true(is.numeric(tpms))
})

test_that("get_mbcp_tpms throws error when samples argument is not a character vector", {
  expect_error(get_mbcp_tpms(samples = 1), "samples should be a character vector")
  expect_error(get_mbcp_tpms(samples = FALSE), "samples should be a character vector")
  expect_error(get_mbcp_tpms(samples = sample(1:100, 5)), "samples should be a character vector")
})

test_that("get_mbcp_tpms throws error when genes argument is not a character vector", {
  expect_error(get_mbcp_tpms(genes = 1), "genes should be a character vector")
  expect_error(get_mbcp_tpms(genes = FALSE), "genes should be a character vector")
  expect_error(get_mbcp_tpms(genes = sample(1:100, 5)), "genes should be a character vector")
})

test_that("get_mbcp_tpms throws error when samples argument includes
          none existing samples", {
  expect_error(get_mbcp_tpms(samples = c("invalid sample", "invalid sample 2")),
               "invalid sample selection")
})

test_that("get_mbcp_tpms throws error when genes argument includes
          none existing genes", {
  expect_error(get_mbcp_tpms(genes = c("invalid gene", "invalid gene 2")),
               "invalid gene selection")
})

test_that("get_mbcp_tpms returns correct dimensions with valid
          samples and genes", {
  ss <- sample(colnames(mbcp_log2_tpms), 5)
  gs <- sample(rownames(mbcp_log2_tpms), 5)
  tpms <- get_mbcp_tpms(genes = gs, samples = ss)
  expect_equal(dim(tpms), c(5, 5))
  expect_equal(colnames(tpms), ss)
  expect_equal(rownames(tpms), gs)
})

test_that("get_mbcp_tpms returns correct dimensions with valid sample and
          valid gene for 1x1", {
  ss <- sample(colnames(mbcp_log2_tpms), 1)
  gs <- sample(rownames(mbcp_log2_tpms), 1)
  tpms <- get_mbcp_tpms(genes = gs, samples = ss)
  expect_equal(dim(tpms), c(1, 1))
  expect_equal(colnames(tpms), ss)
  expect_equal(rownames(tpms), gs)
})

test_that("get_mbcp_tpms_uq returns a numeric matrix", {
  tpms <- get_mbcp_tpms_uq()
  expect_true(is.numeric(tpms))
})

test_that("get_mbcp_tpms_uq throws error when samples argument is not a character vector", {
  expect_error(get_mbcp_tpms_uq(samples = 1))
  expect_error(get_mbcp_tpms_uq(samples = FALSE))
  expect_error(get_mbcp_tpms_uq(samples = sample(1:100, 5)))
})

test_that("get_mbcp_tpms_uq throws error when genes argument is not a character vector", {
  expect_error(get_mbcp_tpms_uq(genes = 1))
  expect_error(get_mbcp_tpms_uq(genes = FALSE))
  expect_error(get_mbcp_tpms_uq(genes = sample(1:100, 5)))
})

test_that("get_mbcp_tpms_uq throws error when samples argument includes
          none existing samples", {
  expect_error(get_mbcp_tpms_uq(samples = c("invalid sample", "invalid sample 2")),
               "invalid sample selection")
})

test_that("get_mbcp_tpms_uq throws error when genes argument includes
          none existing genes", {
  expect_error(get_mbcp_tpms_uq(genes = c("invalid gene", "invalid gene 2")),
               "invalid gene selection")
})

test_that("get_mbcp_tpms_uq returns correct dimensions with valid
          samples and genes", {
  ss <- sample(colnames(mbcp_log2_tpms), 5)
  gs <- sample(rownames(mbcp_log2_tpms), 5)
  tpms <- get_mbcp_tpms_uq(genes = gs, samples = ss)
  expect_equal(dim(tpms), c(5, 5))
  expect_equal(colnames(tpms), ss)
  expect_equal(rownames(tpms), gs)
})

test_that("get_mbcp_tpms_uq calls get_mbcp_tpms", {
  gtpms <- mock(mocked_tpms)

  with_mocked_bindings(code = {
    get_mbcp_tpms_uq()
    expect_called(gtpms, 1)
  }, get_mbcp_tpms = gtpms)
})


## Should retrieve the entire matrix, apply normalization then filter
test_that("get_mbcp_tpms_uq calls get_mbcp_tpms with genes = NULL and samples = NULL", {
  gtpms <- mock(mocked_tpms)
  ss <- sample(colnames(mocked_tpms), 5)
  gs <- sample(rownames(mocked_tpms), 10)

  with_mocked_bindings(code = {
    get_mbcp_tpms_uq(genes = gs, samples = ss)
    margs <- mock_args(gtpms)
    expect_called(gtpms, 1)
    expect_equal(margs[[1]], list(genes = NULL, samples = NULL))
  }, get_mbcp_tpms = gtpms)
})

test_that("get_mbcp_tpms_up returns correct output", {
  gtpms <- mock(mat)
  with_mocked_bindings(code = {
    rtpms <- get_mbcp_tpms_uq()
    expect_identical(rtpms, result)
  }, get_mbcp_tpms = gtpms)
})

test_that("get_mbcp_tpms_uq returns correct output when filtered", {
  gtpms <- mock(mat)
  with_mocked_bindings(code = {
    rtpms <- get_mbcp_tpms_uq(c("g1", "g3", "g4"), c("s2", "s3", "s5"))
    expect_identical(rtpms, result[c("g1", "g3", "g4"), c("s2", "s3", "s5")])
  }, get_mbcp_tpms = gtpms)
})

test_that("get_mbcp_tpms_uq returns correct output when filtered for 1x1", {

  gtpms <- mock(mat)
  with_mocked_bindings(code = {
    rtpms <- get_mbcp_tpms_uq(c("g1"), c("s2"))
    expect_identical(rtpms, result[c("g1"), c("s2"), drop = F])
  }, get_mbcp_tpms = gtpms)
})

test_that("to_longer returns a data frame", {
  df <- to_longer(mocked_tpms, "hugo_symbol")
  expect_s3_class(df, "data.frame")
})

test_that("to_longer throws error when first argument is not a numeric matrix", {
  expect_error(to_longer("not a matrix", "hugo_symbol"), "mbcp_matrix should be a numeric matrix")
  expect_error(to_longer(sample(1:100, 40), "hugo_symbol"), "mbcp_matrix should be a numeric matrix")
  expect_error(to_longer(NULL, "hugo_symbol"), "mbcp_matrix should be a numeric matrix")
  expect_error(to_longer(matrix(data = letters[1:6], nrow = 2), "hugo_symbol"))
})

test_that("to_longer throws error when second argument is not valid id_column", {
  expect_error(to_longer(mocked_tpms, 3))
  expect_error(to_longer(mocked_tpms, sample(1:10, 3)))
  expect_error(to_longer(mocked_tpms, "invalid id_column"))
  expect_error(to_longer(mocked_tpms, NULL))
  expect_error(to_longer(mocked_tpms, FALSE))
})

test_that("to_longer returns a dataframe with correct columns", {
  df <- to_longer(mocked_tpms, "hugo_symbol")
  expect_equal(colnames(df), c("hugo_symbol", "sample_alias", "value"))
})

test_that("to_longer returns a dataframe with correct dimensions", {
  df <- to_longer(mocked_tpms, "hugo_symbol")
  expect_equal(nrow(df), nrow(mocked_tpms)*ncol(mocked_tpms))
  expect_equal(ncol(df), 3)
})

test_that("to_longer returns a dataframe with correct data", {
  df <- to_longer(mocked_tpms, "hugo_symbol")
  expect_true(all(rownames(mocked_tpms) %in% df$hugo_symbol))
  expect_true(all(colnames(mocked_tpms) %in% df$sample_alias))
})
