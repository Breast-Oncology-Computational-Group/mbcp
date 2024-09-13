test_that("get_mbcp_labels returns a character", {
  for(set in pkg_env$set_names) {
    expect_type(get_mbcp_labels(set), "character")
  }
})

test_that("get_mbcp_labels returns a character vector", {
  for(set in pkg_env$set_names) {
    expect_vector(get_mbcp_labels(set), ptype = character())
  }
})

test_that("get_mbcp_labels returns a named vector", {
  for(set in pkg_env$set_names) {
    expect_true(!is.null(names(get_mbcp_labels(set))))
  }
})
