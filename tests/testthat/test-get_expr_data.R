test_that("get_mbcp_log2p1_tpms returns a numeric matrix", {
  tpms <- get_mbcp_log2p1_tpms()
  expect_true(is.numeric(tpms))
})

test_that("get_mbcp_log2p1_tpms returns the entire dataset", {
  tpms <- get_mbcp_log2p1_tpms()
  expect_equal(mbcp_log2p1_tpms, tpms)
})

test_that("get_mbcp_log2p1_tpms throws error when samples argument is not a character vector", {
  expect_error(get_mbcp_log2p1_tpms(samples = 1), "samples should be a character vector")
  expect_error(get_mbcp_log2p1_tpms(samples = FALSE), "samples should be a character vector")
  expect_error(get_mbcp_log2p1_tpms(samples = sample(1:100, 5)), "samples should be a character vector")
})

test_that("get_mbcp_log2p1_tpms throws error when genes argument is not a character vector", {
  expect_error(get_mbcp_log2p1_tpms(genes = 1), "genes should be a character vector")
  expect_error(get_mbcp_log2p1_tpms(genes = FALSE), "genes should be a character vector")
  expect_error(get_mbcp_log2p1_tpms(genes = sample(1:100, 5)), "genes should be a character vector")
})

test_that("get_mbcp_log2p1_tpms throws error when genes argument includes
          none existing genes", {
  expect_error(get_mbcp_log2p1_tpms(genes = c("invalid gene", "invalid gene 2")),
               "invalid gene selection")
})

test_that("get_mbcp_log2p1_tpms throws error when samples argument includes
          none existing samples", {
  expect_error(get_mbcp_log2p1_tpms(samples = c("invalid sample", "invalid sample 2")),
               "invalid sample selection")
})

test_that("get_mbcp_log2p1_tpms returns correct dimensions with valid
          samples and genes", {
  ss <- sample(colnames(mbcp_log2p1_tpms), 5)
  gs <- sample(rownames(mbcp_log2p1_tpms), 5)
  tpms <- get_mbcp_log2p1_tpms(genes = gs, samples = ss)
  expect_equal(dim(tpms), c(5, 5))
  expect_equal(colnames(tpms), ss)
  expect_equal(rownames(tpms), gs)
})

test_that("get_mbcp_log2p1_tpms returns correct dimensions with valid sample and
          valid gene for 1x1", {
  ss <- sample(colnames(mbcp_log2p1_tpms), 1)
  gs <- sample(rownames(mbcp_log2p1_tpms), 1)
  tpms <- get_mbcp_log2p1_tpms(genes = gs, samples = ss)
  expect_equal(dim(tpms), c(1, 1))
  expect_equal(colnames(tpms), ss)
  expect_equal(rownames(tpms), gs)
})

test_that("get_mbcp_log2p1_tpms_uq returns a numeric matrix", {
  tpms <- get_mbcp_log2p1_tpms_uq()
  expect_true(is.numeric(tpms))
})

test_that("get_mbcp_log2p1_tpms_uq throws error when samples argument is not a character vector", {
  expect_error(get_mbcp_log2p1_tpms_uq(samples = 1))
  expect_error(get_mbcp_log2p1_tpms_uq(samples = FALSE))
  expect_error(get_mbcp_log2p1_tpms_uq(samples = sample(1:100, 5)))
})

test_that("get_mbcp_log2p1_tpms_uq throws error when genes argument is not a character vector", {
  expect_error(get_mbcp_log2p1_tpms_uq(genes = 1))
  expect_error(get_mbcp_log2p1_tpms_uq(genes = FALSE))
  expect_error(get_mbcp_log2p1_tpms_uq(genes = sample(1:100, 5)))
})

test_that("get_mbcp_log2p1_tpms_uq throws error when samples argument includes
          none existing samples", {
  expect_error(get_mbcp_log2p1_tpms_uq(samples = c("invalid sample", "invalid sample 2")),
               "invalid sample selection")
})

test_that("get_mbcp_log2p1_tpms_uq throws error when genes argument includes
          none existing genes", {
  expect_error(get_mbcp_log2p1_tpms_uq(genes = c("invalid gene", "invalid gene 2")),
               "invalid gene selection")
})

test_that("get_mbcp_log2p1_tpms_uq returns correct dimensions with valid
          samples and genes", {
  ss <- sample(colnames(mbcp_log2p1_tpms), 5)
  gs <- sample(rownames(mbcp_log2p1_tpms), 5)
  tpms <- get_mbcp_log2p1_tpms_uq(genes = gs, samples = ss)
  expect_equal(dim(tpms), c(5, 5))
  expect_equal(colnames(tpms), ss)
  expect_equal(rownames(tpms), gs)
})

test_that("get_mbcp_log2p1_tpms_uq calls get_mbcp_log2p1_tpms", {
  gtpms <- mock(mocked_tpms)

  with_mocked_bindings(code = {
    get_mbcp_log2p1_tpms_uq()
    expect_called(gtpms, 1)
  }, get_mbcp_log2p1_tpms = gtpms)
})


## Should retrieve the entire matrix, apply normalization then filter
test_that("get_mbcp_log2p1_tpms_uq calls get_mbcp_log2p1_tpms with genes = NULL and samples = NULL", {
  gtpms <- mock(mocked_tpms)
  ss <- sample(colnames(mocked_tpms), 5)
  gs <- sample(rownames(mocked_tpms), 10)

  with_mocked_bindings(code = {
    get_mbcp_log2p1_tpms_uq(genes = gs, samples = ss)
    margs <- mock_args(gtpms)
    expect_called(gtpms, 1)
    expect_equal(margs[[1]], list(genes = NULL, samples = NULL))
  }, get_mbcp_log2p1_tpms = gtpms)
})

test_that("get_mbcp_log2p1_tpms_up returns correct output", {
  gtpms <- mock(mocked_tpms_min)
  with_mocked_bindings(code = {
    rtpms <- get_mbcp_log2p1_tpms_uq()
    expect_identical(rtpms, mocked_tpms_min_uq)
  }, get_mbcp_log2p1_tpms = gtpms)
})

test_that("get_mbcp_log2p1_tpms_uq returns correct output when filtered", {
  gtpms <- mock(mocked_tpms_min)
  with_mocked_bindings(code = {
    rtpms <- get_mbcp_log2p1_tpms_uq(c("g1", "g3", "g4"), c("s2", "s3", "s5"))
    expect_identical(rtpms, mocked_tpms_min_uq[c("g1", "g3", "g4"), c("s2", "s3", "s5")])
  }, get_mbcp_log2p1_tpms = gtpms)
})

test_that("get_mbcp_log2p1_tpms_uq returns correct output when filtered for 1x1", {

  gtpms <- mock(mocked_tpms_min)
  with_mocked_bindings(code = {
    rtpms <- get_mbcp_log2p1_tpms_uq(c("g1"), c("s2"))
    expect_identical(rtpms, mocked_tpms_min_uq[c("g1"), c("s2"), drop = F])
  }, get_mbcp_log2p1_tpms = gtpms)
})

test_that("to_longer returns a data frame", {
  df <- to_longer(mocked_tpms, "hugo_symbol")
  expect_s3_class(df, "data.frame")
})

test_that("to_longer throws error when first argument is not a numeric matrix", {
  expect_error(to_longer("not a matrix", "hugo_symbol"), "mbcp_matrix should be a numeric or a character matrix")
  expect_error(to_longer(sample(1:100, 40), "hugo_symbol"), "mbcp_matrix should be a numeric or a character matrix")
  expect_error(to_longer(NULL, "hugo_symbol"), "mbcp_matrix should be a numeric or a character matrix")
  expect_error(to_longer(matrix(data = sample(c(TRUE, FALSE), 6, replace = T), nrow = 2), "hugo_symbol"))
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

test_that("get_mbcp_log2p1_tpms_class returns a character matrix", {
  tpms_class <- get_mbcp_log2p1_tpms_class()
  expect_true(is.character(tpms_class))
})

test_that("get_mbcp_log2p1_tpms_class returns the entire dataset", {
  tpms_class <- get_mbcp_log2p1_tpms_class()
  expect_equal(mbcp_log2p1_tpms_class, tpms_class)
})

test_that("get_mbcp_log2p1_tpms_class throws error when samples argument is not a character vector", {
  expect_error(get_mbcp_log2p1_tpms_class(samples = 1), "samples should be a character vector")
  expect_error(get_mbcp_log2p1_tpms_class(samples = FALSE), "samples should be a character vector")
  expect_error(get_mbcp_log2p1_tpms_class(samples = sample(1:100, 5)), "samples should be a character vector")
})

test_that("get_mbcp_log2p1_tpms_class throws error when genes argument is not a character vector", {
  expect_error(get_mbcp_log2p1_tpms_class(genes = 1), "genes should be a character vector")
  expect_error(get_mbcp_log2p1_tpms_class(genes = FALSE), "genes should be a character vector")
  expect_error(get_mbcp_log2p1_tpms_class(genes = sample(1:100, 5)), "genes should be a character vector")
})

test_that("get_mbcp_log2p1_tpms_class throws error when genes argument includes
          none existing genes", {
            expect_error(get_mbcp_log2p1_tpms_class(genes = c("invalid gene", "invalid gene 2")),
                         "invalid gene selection")
})

test_that("get_mbcp_log2p1_tpms_class throws error when samples argument includes
          none existing samples", {
            expect_error(get_mbcp_log2p1_tpms_class(samples = c("invalid sample", "invalid sample 2")),
                         "invalid sample selection")
})

test_that("get_mbcp_log2p1_tpms_class returns correct dimensions with valid samples and genes", {
    ss <- sample(colnames(mbcp_log2p1_tpms_class), 5)
    gs <- sample(rownames(mbcp_log2p1_tpms_class), 5)
    tpms <- get_mbcp_log2p1_tpms_class(gs, ss)
    expect_equal(dim(tpms), c(5, 5))
    expect_equal(colnames(tpms), ss)
    expect_equal(rownames(tpms), gs)
})

test_that("get_mbcp_log2p1_tpms_class returns correct output when filtered for 1x1", {
  ss <- sample(colnames(mbcp_log2p1_tpms_class), 1)
  gs <- sample(rownames(mbcp_log2p1_tpms_class), 1)
  tpms <- get_mbcp_log2p1_tpms_class(gs, ss)
  expect_equal(dim(tpms), c(1, 1))
  expect_equal(colnames(tpms), ss)
  expect_equal(rownames(tpms), gs)
})
