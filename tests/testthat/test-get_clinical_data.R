test_that("get_clinical_data_all_variables returns a dataframe", {
  expect_s3_class(get_clinical_data_all_variables(), "data.frame")
})

test_that("get_clinical_Data returns a data.frame", {
  expect_s3_class(get_clinical_data(), "data.frame")
})

test_that("get_clinical_data calls get_clinical_data_all_variables()", {
  gcd_all = mock(mocked_clinical_data)

  with_mocked_bindings(code = {
    get_clinical_data()
    expect_called(gcd_all, 1)
  }, get_clinical_data_all_variables = gcd_all)
})

test_that("get_clinical_data returns correct columns", {
  df = get_clinical_data()
  exp_colnames = c("sample_alias", "patient_id", "sample_timepoint", "wes_sample_id",
  "receptor_status_dx_all", "hr_status_dx_all", "dx_histology", "pam50_subtype",
  "bx_location", "bx_time_days", "calc_time_to_mets_dx_days",
  "calc_met_setting", "calc_primary_treat", "purity", "timepoint", "n_bypatient")
  expect_equal(colnames(df), exp_colnames)
})


test_that("get_clinical_data_rna calls get_clinical_data()", {
  gcd = mock(mocked_clinical_data)

  with_mocked_bindings(code = {
    get_clinical_data_rna()
    expect_called(gcd, 1)
  }, get_clinical_data = gcd)
})

test_that("get_clinical_data_rna returns only entries with pam50 subtype", {
  rna_data <- get_clinical_data_rna()
  expect_true(all(!is.na(rna_data$pam50_subtype)))
})

test_that("get_clinical_data_rna returns entries with no wes_sample_id", {
  rna_data <- get_clinical_data_rna()
  expect_true(any(is.na(rna_data$wes_sample_id)))
})

test_that("get_clinical_data_rna returns only entries with wes_sample_id when wes=TRUE", {
  rna_data <- get_clinical_data_rna(wes = TRUE)
  expect_true(all(!is.na(rna_data$wes_sample_id)))
})

test_that("get_clinical_data_rna returns no pam50_subtype Normal with filter='normals'", {
  rna_data <- get_clinical_data_rna(filter = "normals")
  expect_true(all(rna_data$pam50_subtype != "Normal"))
})

test_that("get_clinical_data_rna return no pam50_subtype Normal and no NA wes_sample_id
          with filter='normals' and wes = TRUE", {
  rna_data <- get_clinical_data_rna(filter = "normals", wes = TRUE)
  expect_true(all(rna_data$pam50_subtype != "Normal" & !is.na(rna_data$wes_sample_id)))
})

test_that("get_clinical_data_rna returns no pam50_subtype Normal and no entries with
          purity < 0.272 with filter='purity'", {
            rna_data <- get_clinical_data_rna(filter = "purity")
            expect_true(all(rna_data$pam50_subtype != "Normal" &
                              !rna_data$purity < 0.272))
})

test_that("get_clinical_data_rna returns no pam50_subtype Normal, no entries with
          purity < 0.272 and no NA wes_sample_id with filter='purity' and wes = TRUE", {
  rna_data <- get_clinical_data_rna(filter = "purity", wes = TRUE)
  expect_true(all(rna_data$pam50_subtype != "Normal" & !is.na(rna_data$wes_sample_id) &
                    !rna_data$purity < 0.272))
})

test_that("get_clinical_data_rna returns no pam50_subtype Normal and no entries with
          bx_location %in% ASCENDING COLON, BONE , LUNG =, PLEURAL FLUID with filter = 'tissues'", {

  removed_tissues <- c("ASCENDING COLON","BONE","LUNG", "PLEURAL FLUID")
  rna_data <- get_clinical_data_rna(filter = "tissues")
  expect_true(all(rna_data$pam50_subtype != "Normal"))
  expect_false(any(rna_data$bx_location %in% removed_tissues))
})

test_that("get_clinical_data_rna returns no pam50_subtype Normal, no entries with
          bx_location %in% ASCENDING COLON, BONE , LUNG =, PLEURAL FLUID with filter = 'tissues'
          and no NA wes_sample_id with filter='purity' and wes = TRUE", {

  removed_tissues <- c("ASCENDING COLON","BONE","LUNG", "PLEURAL FLUID")
  rna_data <- get_clinical_data_rna(filter = "tissues", wes = TRUE)
  expect_false(any(rna_data$bx_location %in% removed_tissues))
  expect_true(all(rna_data$pam50_subtype != "Normal" & !is.na(rna_data$wes_sample_id)))
})

test_that("get_clinical_data_rna return no pam50_subtype Normal and entries
          with bx_location = BREAST with filter = 'breast'", {
  rna_data <- get_clinical_data_rna(filter = "breast")
  expect_true(all(rna_data$pam50_subtype != "Normal" & rna_data$bx_location == "BREAST"))
})

test_that("get_clinical_data_rna return no pam50_subtype Normal and entries
          with bx_location = BREAST with filter = 'breast'
          and no NA wes_sample_id with filter='purity' and wes = TRUE", {
  rna_data <- get_clinical_data_rna(filter = "breast", wes = TRUE)
  expect_true(all(rna_data$pam50_subtype != "Normal" & rna_data$bx_location == "BREAST" &
                    !is.na(rna_data$wes_sample_id)))
})

test_that("get_clinical_data_rna throws error when filter is not one of 'normals', 'purity', 'tissues', 'breast'", {
  expect_error(get_clinical_data_rna(filter = "not_a_filter"), "not_a_filter is not a valid filter argument")
})

test_that("get_clinical_data_rna throws error when wes is not logical", {
  expect_error(get_clinical_data_rna(wes = "not_logical"), "wes should be a logical value")
})

test_that("get_clinical_data_wes calls get_clinical_data()", {
  gcd = mock(mocked_clinical_data)

  with_mocked_bindings(code = {
    get_clinical_data_wes()
    expect_called(gcd, 1)
  }, get_clinical_data = gcd)
})

test_that("get_clinical_data_wes returns dataset with no NA wes_sample_id", {
  wes_data <- get_clinical_data_wes()
  expect_true(all(!is.na(wes_data$wes_sample_id)))
})

test_that("get_clinical_data_first_samples calls get_clinical_data() if sample_type == any", {
  gcd = mock(mocked_clinical_data)

  with_mocked_bindings(code = {
    get_clinical_data_first_samples(sample_type = "any")
    expect_called(gcd, 1)
  }, get_clinical_data = gcd)
})

test_that("get_clinical_data_first_samples returns one entry per patient and timepoint == 1", {
  df <- get_clinical_data_first_samples("any")
  expect_equal(length(unique(df$patient_id)), nrow(df))
  expect_true(all(df$timepoint == 1))
})

test_that("get_clinical_Data_first_samples returns ungrouped data for sample_type = any", {
  expect_false(dplyr::is_grouped_df(get_clinical_data_first_samples("any")))
})

test_that("get_clinical_data_first_samples calls get_clinical_data_wes() if sample_type == wes", {
  gcd = mock(mocked_clinical_data)

  with_mocked_bindings(code = {
    get_clinical_data_first_samples(sample_type = "wes")
    expect_called(gcd, 1)
  }, get_clinical_data_wes = gcd)
})

test_that("get_clinical_data_first_samples with sample_type == wes returns one entry per patient and correct timepoint", {
  df <- get_clinical_data_first_samples("wes")
  expect_equal(length(unique(df$patient_id)), nrow(df))

  left_out_clinical_data <- get_clinical_data() |>
    dplyr::anti_join(df, by = "sample_alias") |>
    dplyr::inner_join(df, by = dplyr::join_by(patient_id == patient_id, timepoint < timepoint))

  expect_true(all(is.na(left_out_clinical_data$wes_sample_id.x)))
})

test_that("get_clinical_data_first_samples returns ungrouped data for sample_type = wes", {
  expect_false(dplyr::is_grouped_df(get_clinical_data_first_samples("wes")))
})

test_that("get_clinical_data_first_samples calls get_clinical_data_rna() if sample_type == rna", {
  gcd = mock(mocked_clinical_data)

  with_mocked_bindings(code = {
    get_clinical_data_first_samples(sample_type = "rna")
    expect_called(gcd, 1)
    expect_equal(mock_args(gcd)[[1]], list())
  }, get_clinical_data_rna = gcd)
})

test_that("get_clinical_data_first_samples with sample_type == rna returns one entry per patient and correct timepoint", {
  df <- get_clinical_data_first_samples("rna")
  expect_equal(length(unique(df$patient_id)), nrow(df))

  left_out_clinical_data <- get_clinical_data() |>
    dplyr::anti_join(df, by = "sample_alias") |>
    dplyr::inner_join(df, by = dplyr::join_by(patient_id == patient_id, timepoint < timepoint))

  expect_true(all(is.na(left_out_clinical_data$pam50_subtype.x)))
})

test_that("get_clinical_data_first_samples returns ungrouped data for sample_type = rna", {
  expect_false(dplyr::is_grouped_df(get_clinical_data_first_samples("rna")))
})

test_that("get_clinical_data_first_samples calls get_clinical_data_rna() if sample_type == wes_rna", {
  gcd = mock(mocked_clinical_data)

  with_mocked_bindings(code = {
    get_clinical_data_first_samples(sample_type = "wes_rna")
    expect_called(gcd, 1)
    expect_equal(mock_args(gcd)[[1]], list(wes = TRUE))
  }, get_clinical_data_rna = gcd)
})

test_that("get_clinical_data_first_samples with sample_type == rna returns one entry per patient and correct timepoint", {
  df <- get_clinical_data_first_samples("wes_rna")
  expect_equal(length(unique(df$patient_id)), nrow(df))

  left_out_clinical_data <- get_clinical_data() |>
    dplyr::anti_join(df, by = "sample_alias") |>
    dplyr::inner_join(df, by = dplyr::join_by(patient_id == patient_id, timepoint < timepoint))

  expect_true(all(is.na(left_out_clinical_data$wes_sample_id.x) | is.na(left_out_clinical_data$pam50_subtype.x)))
})

test_that("get_clinical_data_first_samples returns ungrouped data for sample_type = wes_rna", {
  expect_false(dplyr::is_grouped_df(get_clinical_data_first_samples("wes_rna")))
})

test_that("get_clinical_data_first_samples throws error when invalid sample_type", {
  expect_error(get_clinical_data_first_samples("invalid sample type"), "sample_type should be one of: any, wes, rna or wes_rna")
  expect_error(get_clinical_data_first_samples(34), "sample_type should be one of: any, wes, rna or wes_rna")
})

test_that("get_clinical_data_latest_samples throws error when invalid sample_type", {
  expect_error(get_clinical_data_latest_samples("invalid sample type"), "sample_type should be one of: any, wes, rna or wes_rna")
  expect_error(get_clinical_data_latest_samples(34), "sample_type should be one of: any, wes, rna or wes_rna")
})

test_that("get_clinical_data_latest_samples calls get_clinical_data() if sample_type == any", {
  gcd = mock(mocked_clinical_data)

  with_mocked_bindings(code = {
    get_clinical_data_latest_samples(sample_type = "any")
    expect_called(gcd, 1)
  }, get_clinical_data = gcd)
})

test_that("get_clinical_data_latest_samples returns ungrouped data for sample_type = any", {
  expect_false(dplyr::is_grouped_df(get_clinical_data_latest_samples("any")))
})

test_that("get_clinical_data_latest_samples returns one entry per patient and timepoint == n_bypatient", {
  df <- get_clinical_data_latest_samples("any")
  expect_equal(length(unique(df$patient_id)), nrow(df))
  expect_true(all(df$timepoint == df$n_bypatient))
})

test_that("get_clinical_data_latest_samples calls get_clinical_data_wes() if sample_type == wes", {
  gcd = mock(mocked_clinical_data)

  with_mocked_bindings(code = {
    get_clinical_data_latest_samples(sample_type = "wes")
    expect_called(gcd, 1)
  }, get_clinical_data_wes = gcd)
})

test_that("get_clinical_data_latest_samples returns ungrouped data for sample_type = wes", {
  expect_false(dplyr::is_grouped_df(get_clinical_data_latest_samples("wes")))
})

test_that("get_clinical_data_latest_samples with sample_type == wes returns one entry per patient and correct timepoint", {
  df <- get_clinical_data_latest_samples("wes")
  expect_equal(length(unique(df$patient_id)), nrow(df))

  left_out_clinical_data <- get_clinical_data() |>
    dplyr::anti_join(df, by = "sample_alias") |>
    dplyr::inner_join(df, by = dplyr::join_by(patient_id == patient_id, timepoint > timepoint))

  expect_true(all(is.na(left_out_clinical_data$wes_sample_id.x)))
})

test_that("get_clinical_data_latest_samples calls get_clinical_data_rna() if sample_type == rna", {
  gcd = mock(mocked_clinical_data)

  with_mocked_bindings(code = {
    get_clinical_data_latest_samples(sample_type = "rna")
    expect_called(gcd, 1)
    expect_equal(mock_args(gcd)[[1]], list())
  }, get_clinical_data_rna = gcd)
})

test_that("get_clinical_data_latest_samples with sample_type == rna returns one entry per patient and correct timepoint", {
  df <- get_clinical_data_latest_samples("rna")
  expect_equal(length(unique(df$patient_id)), nrow(df))

  left_out_clinical_data <- get_clinical_data() |>
    dplyr::anti_join(df, by = "sample_alias") |>
    dplyr::inner_join(df, by = dplyr::join_by(patient_id == patient_id, timepoint > timepoint))

  expect_true(all(is.na(left_out_clinical_data$pam50_subtype.x)))
})

test_that("get_clinical_data_latest_samples returns ungrouped data for sample_type = rna", {
  expect_false(dplyr::is_grouped_df(get_clinical_data_latest_samples("rna")))
})

test_that("get_clinical_data_latest_samples calls get_clinical_data_rna() if sample_type == wes_rna", {
  gcd = mock(mocked_clinical_data)

  with_mocked_bindings(code = {
    get_clinical_data_latest_samples(sample_type = "wes_rna")
    expect_called(gcd, 1)
    expect_equal(mock_args(gcd)[[1]], list(wes = TRUE))
  }, get_clinical_data_rna = gcd)
})

test_that("get_clinical_data_latest_samples with sample_type == rna returns one entry per patient and correct timepoint", {
  df <- get_clinical_data_latest_samples("wes_rna")
  expect_equal(length(unique(df$patient_id)), nrow(df))

  left_out_clinical_data <- get_clinical_data() |>
    dplyr::anti_join(df, by = "sample_alias") |>
    dplyr::inner_join(df, by = dplyr::join_by(patient_id == patient_id, timepoint > timepoint))

  expect_true(all(is.na(left_out_clinical_data$wes_sample_id.x) | is.na(left_out_clinical_data$pam50_subtype.x)))
})

test_that("get_clinical_data_latest_samples returns ungrouped data for sample_type = wes_rna", {
  expect_false(dplyr::is_grouped_df(get_clinical_data_latest_samples("wes_rna")))
})
