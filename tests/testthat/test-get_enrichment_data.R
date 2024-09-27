test_that("get_mbcp_enrichment_scores returns a numeric matrix", {
  scores <- get_mbcp_enrichment_scores()
  expect_true(is.numeric(scores))
})

test_that("get_mbcp_enrichment_scores returns the entire dataset", {
  scores <- get_mbcp_enrichment_scores()
  expect_equal(mbcp_enrichment_scores, scores)
})

test_that("get_mbcp_enrichment_scores throws error when signatures argument is not a character vector", {
  expect_error(get_mbcp_enrichment_scores(signatures = 1), "signatures should be a character vector")
  expect_error(get_mbcp_enrichment_scores(signatures = FALSE), "signatures should be a character vector")
  expect_error(get_mbcp_enrichment_scores(signatures = sample(1:100, 5)), "signatures should be a character vector")
})

test_that("get_mbcp_enrichment_scores throws error when samples argument is not a character vector", {
  expect_error(get_mbcp_enrichment_scores(samples = 1), "samples should be a character vector")
  expect_error(get_mbcp_enrichment_scores(samples = FALSE), "samples should be a character vector")
  expect_error(get_mbcp_enrichment_scores(samples = sample(1:100, 5)), "samples should be a character vector")
})

test_that("get_mbcp_enrichment_scores throws error when samples argument includes
          none existing samples", {
            expect_error(get_mbcp_enrichment_scores(samples = c("invalid sample", "invalid sample 2")),
                         "invalid sample selection")
          })

test_that("get_mbcp_enrichment_scores throws error when signatures argument includes
          none existing signatures", {
            expect_error(get_mbcp_enrichment_scores(signatures = c("invalid signature", "invalid signature 2")),
                         "invalid signature selection")
          })

test_that("get_mbcp_enrichment_scores returns correct dimensions with valid
          samples and signatures", {
  ss <- sample(colnames(mbcp_enrichment_scores), 5)
  gs <- sample(rownames(mbcp_enrichment_scores), 5)
  scores <- get_mbcp_enrichment_scores(signatures = gs, samples = ss)
  expect_equal(dim(scores), c(5, 5))
  expect_equal(colnames(scores), ss)
  expect_equal(rownames(scores), gs)
})

test_that("get_mbcp_enrichment_scores returns correct dimensions with valid sample and
          valid signature for 1x1", {
  ss <- sample(colnames(mbcp_enrichment_scores), 1)
  gs <- sample(rownames(mbcp_enrichment_scores), 1)
  scores <- get_mbcp_enrichment_scores(signatures = gs, samples = ss)
  expect_equal(dim(scores), c(1, 1))
  expect_equal(colnames(scores), ss)
  expect_equal(rownames(scores), gs)
})

test_that("get_mbcp_enrichment_class returns a character matrix", {
  class_data <- get_mbcp_enrichment_class()
  expect_true(is.character(class_data))
})

test_that("get_mbcp_enrichment_class returns the entire dataset", {
  class_data <- get_mbcp_enrichment_class()
  expect_equal(mbcp_enrichment_class, class_data)
})

test_that("get_mbcp_enrichment_class throws error when samples argument is not a character vector", {
  expect_error(get_mbcp_enrichment_class(samples = 1), "samples should be a character vector")
  expect_error(get_mbcp_enrichment_class(samples = FALSE), "samples should be a character vector")
  expect_error(get_mbcp_enrichment_class(samples = sample(1:100, 5)), "samples should be a character vector")
})

test_that("get_mbcp_enrichment_class throws error when signatures argument is not a character vector", {
  expect_error(get_mbcp_enrichment_class(signatures = 1), "signatures should be a character vector")
  expect_error(get_mbcp_enrichment_class(signatures = FALSE), "signatures should be a character vector")
  expect_error(get_mbcp_enrichment_class(signatures = sample(1:100, 5)), "signatures should be a character vector")
})

test_that("get_mbcp_enrichment_class throws error when samples argument includes
          none existing samples", {
            expect_error(get_mbcp_enrichment_class(samples = c("invalid sample", "invalid sample 2")),
                         "invalid sample selection")
})

test_that("get_mbcp_enrichment_class throws error when signatures argument includes
          none existing signatures", {
            expect_error(get_mbcp_enrichment_class(signatures = c("invalid signature", "invalid signature 2")),
                         "invalid signature selection")
})

test_that("get_mbcp_enrichment_class returns correct dimensions with valid
          samples and signatures", {
            ss <- sample(colnames(mbcp_enrichment_class), 5)
            gs <- sample(rownames(mbcp_enrichment_class), 5)
            class_data <- get_mbcp_enrichment_class(gs, ss)
            expect_equal(dim(class_data), c(5, 5))
            expect_equal(colnames(class_data), ss)
            expect_equal(rownames(class_data), gs)
})

test_that("get_mbcp_enrichment_class returns correct dimensions with valid sample and
          valid signature for 1x1", {
            ss <- sample(colnames(mbcp_enrichment_class), 1)
            gs <- sample(rownames(mbcp_enrichment_class), 1)
            class_data <- get_mbcp_enrichment_class(gs, ss)
            expect_equal(dim(class_data), c(1, 1))
            expect_equal(colnames(class_data), ss)
            expect_equal(rownames(class_data), gs)
})


