set_names <- c("receptor", "pam50", "dna_alts", "mutation", "histology",
               "treatment", "metastasis", "expression", "germline_alts", "dna_signatures")

test_that("get_mbcp_colors returns a character", {
  for (set in set_names) {
    expect_type(get_mbcp_colors(set), "character")
  }
})

test_that("get_mbcp_colors returns a character vector", {
  for (set in set_names) {
    expect_vector(get_mbcp_colors(set), ptype = character())
  }
})

test_that("get_mbcp_colors returns a character vector with names", {
  for (set in set_names) {
    expect_true(!is.null(names(get_mbcp_colors(set))))
  }
})

test_that("get_mbcp_colors throws error with non character input", {
  expect_error(get_mbcp_colors(1), "set must be a character")
  expect_error(get_mbcp_colors(1.1), "set must be a character")
  expect_error(get_mbcp_colors(TRUE), "set must be a character")
  expect_error(get_mbcp_colors(FALSE), "set must be a character")
  expect_error(get_mbcp_colors(NULL), "set must be a character")
  expect_error(get_mbcp_colors(NA), "set must be a character")
})

test_that("get_mbcp_colors throws error with non existing set", {
  expect_error(get_mbcp_colors("non_existing_set"), "non_existing_set is not a valid palette name for MBC project")
})

