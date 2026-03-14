# ---- Test data setup ----
make_test_data <- function() {
  tax <- build_tax_tree(
    species = c("sp1", "sp2", "sp3", "sp4", "sp5"),
    Genus   = c("G1", "G1", "G2", "G2", "G3"),
    Family  = c("F1", "F1", "F1", "F2", "F2"),
    Order   = c("O1", "O1", "O1", "O1", "O1")
  )
  comm <- c(sp1 = 10, sp2 = 5, sp3 = 8, sp4 = 2, sp5 = 3)
  list(tax = tax, comm = comm)
}


# ---- compare_indices ----

test_that("compare_indices returns S3 class", {

  d <- make_test_data()
  result <- compare_indices(d$comm, d$tax)
  expect_s3_class(result, "compare_indices")
  expect_s3_class(result, "data.frame")
})

test_that("compare_indices print method works", {
  d <- make_test_data()
  result <- compare_indices(d$comm, d$tax)
  expect_output(print(result), "taxdiv -- Index Comparison")
})

test_that("compare_indices summary method works", {
  d <- make_test_data()
  comms <- list(A = d$comm, B = c(sp1 = 5, sp2 = 5, sp3 = 5, sp4 = 5, sp5 = 5))
  result <- compare_indices(comms, d$tax)
  expect_output(summary(result), "Summary")
})

test_that("compare_indices backward compatible with data.frame ops", {
  d <- make_test_data()
  result <- compare_indices(d$comm, d$tax)
  expect_true(nrow(result) >= 1)
  expect_true("Shannon" %in% names(result))
  expect_true(is.numeric(result$Shannon))
  expect_true(is.data.frame(result))
  # Subsetting should work
  sub <- result[, c("Community", "Shannon")]
  expect_true(is.data.frame(sub))
})

test_that("compare_indices plot method works", {
  skip_if_not_installed("ggplot2")
  d <- make_test_data()
  result <- compare_indices(d$comm, d$tax)
  expect_silent(p <- plot(result))
  expect_s3_class(p, "gg")
})


# ---- batch_analysis ----

test_that("batch_analysis returns S3 class", {
  df <- data.frame(
    Site      = c("A", "A", "A", "B", "B", "B"),
    Species   = c("sp1", "sp2", "sp3", "sp1", "sp3", "sp4"),
    Genus     = c("G1", "G1", "G2", "G1", "G2", "G2"),
    Family    = c("F1", "F1", "F1", "F1", "F1", "F2"),
    Order     = c("O1", "O1", "O1", "O1", "O1", "O1"),
    Abundance = c(10, 20, 15, 5, 25, 10),
    stringsAsFactors = FALSE
  )
  result <- batch_analysis(df)
  expect_s3_class(result, "batch_analysis")
  expect_s3_class(result, "data.frame")
})

test_that("batch_analysis print method works", {
  df <- data.frame(
    Species   = c("sp1", "sp2", "sp3"),
    Genus     = c("G1", "G1", "G2"),
    Family    = c("F1", "F1", "F1"),
    Order     = c("O1", "O1", "O1"),
    Abundance = c(10, 20, 15),
    stringsAsFactors = FALSE
  )
  result <- batch_analysis(df)
  expect_output(print(result), "taxdiv -- Batch Analysis")
})

test_that("batch_analysis summary method works", {
  df <- data.frame(
    Site      = c("A", "A", "A", "B", "B", "B"),
    Species   = c("sp1", "sp2", "sp3", "sp1", "sp3", "sp4"),
    Genus     = c("G1", "G1", "G2", "G1", "G2", "G2"),
    Family    = c("F1", "F1", "F1", "F1", "F1", "F2"),
    Order     = c("O1", "O1", "O1", "O1", "O1", "O1"),
    Abundance = c(10, 20, 15, 5, 25, 10),
    stringsAsFactors = FALSE
  )
  result <- batch_analysis(df)
  expect_output(summary(result), "Summary")
})

test_that("batch_analysis backward compatible with data.frame ops", {
  df <- data.frame(
    Species   = c("sp1", "sp2", "sp3"),
    Genus     = c("G1", "G1", "G2"),
    Family    = c("F1", "F1", "F1"),
    Order     = c("O1", "O1", "O1"),
    Abundance = c(10, 20, 15),
    stringsAsFactors = FALSE
  )
  result <- batch_analysis(df)
  expect_true(nrow(result) >= 1)
  expect_true("Shannon" %in% names(result))
  expect_true(is.data.frame(result))
})


# ---- ozkan_pto ----

test_that("ozkan_pto returns S3 class", {
  d <- make_test_data()
  result <- ozkan_pto(d$comm, d$tax)
  expect_s3_class(result, "ozkan_pto")
})

test_that("ozkan_pto print method works", {
  d <- make_test_data()
  result <- ozkan_pto(d$comm, d$tax)
  expect_output(print(result), "taxdiv -- Ozkan pTO Result")
  expect_output(print(result), "uTO")
})

test_that("ozkan_pto single species returns S3 class", {
  tax <- build_tax_tree(
    species = "sp1",
    Genus   = "G1",
    Family  = "F1"
  )
  result <- ozkan_pto(c(sp1 = 5), tax)
  expect_s3_class(result, "ozkan_pto")
  expect_equal(result$uTO, 0)
})

test_that("ozkan_pto backward compatible with list ops", {
  d <- make_test_data()
  result <- ozkan_pto(d$comm, d$tax)
  expect_true(is.list(result))
  expect_true("uTO" %in% names(result))
  expect_true(is.numeric(result$uTO))
  expect_true(!is.null(result$Ed_levels))
})


# ---- ozkan_pto_resample ----

test_that("ozkan_pto_resample returns S3 class", {
  d <- make_test_data()
  result <- ozkan_pto_resample(d$comm, d$tax, n_iter = 101, seed = 42)
  expect_s3_class(result, "ozkan_pto_resample")
})

test_that("ozkan_pto_resample print method works", {
  d <- make_test_data()
  result <- ozkan_pto_resample(d$comm, d$tax, n_iter = 101, seed = 42)
  expect_output(print(result), "Stochastic Resampling")
  expect_output(print(result), "happy")
})

test_that("ozkan_pto_resample summary method works", {
  d <- make_test_data()
  result <- ozkan_pto_resample(d$comm, d$tax, n_iter = 101, seed = 42)
  expect_output(summary(result), "Summary")
  expect_output(summary(result), "Deterministic")
})

test_that("ozkan_pto_resample backward compatible with list ops", {
  d <- make_test_data()
  result <- ozkan_pto_resample(d$comm, d$tax, n_iter = 101, seed = 42)
  expect_true(is.list(result))
  expect_true("iteration_results" %in% names(result))
  expect_true(is.data.frame(result$iteration_results))
})


# ---- ozkan_pto_sensitivity ----

test_that("ozkan_pto_sensitivity returns S3 class", {
  d <- make_test_data()
  run2 <- ozkan_pto_resample(d$comm, d$tax, n_iter = 101, seed = 42)
  result <- ozkan_pto_sensitivity(d$comm, d$tax, run2, seed = 123)
  expect_s3_class(result, "ozkan_pto_sensitivity")
})

test_that("ozkan_pto_sensitivity print method works", {
  d <- make_test_data()
  run2 <- ozkan_pto_resample(d$comm, d$tax, n_iter = 101, seed = 42)
  result <- ozkan_pto_sensitivity(d$comm, d$tax, run2, seed = 123)
  expect_output(print(result), "Sensitivity Analysis")
})

test_that("ozkan_pto_sensitivity summary method works", {
  d <- make_test_data()
  run2 <- ozkan_pto_resample(d$comm, d$tax, n_iter = 101, seed = 42)
  result <- ozkan_pto_sensitivity(d$comm, d$tax, run2, seed = 123)
  expect_output(summary(result), "Summary")
  expect_output(summary(result), "prob")
})


# ---- rarefaction_taxonomic ----

test_that("rarefaction_taxonomic returns S3 class", {
  d <- make_test_data()
  result <- rarefaction_taxonomic(d$comm, d$tax,
                                   index = "shannon",
                                   steps = 5, n_boot = 10, seed = 1)
  expect_s3_class(result, "rarefaction_taxonomic")
  expect_s3_class(result, "data.frame")
})

test_that("rarefaction_taxonomic print method works", {
  d <- make_test_data()
  result <- rarefaction_taxonomic(d$comm, d$tax,
                                   index = "shannon",
                                   steps = 5, n_boot = 10, seed = 1)
  expect_output(print(result), "taxdiv -- Rarefaction Curve")
  expect_output(print(result), "shannon")
})

test_that("rarefaction_taxonomic summary method works", {
  d <- make_test_data()
  result <- rarefaction_taxonomic(d$comm, d$tax,
                                   index = "shannon",
                                   steps = 5, n_boot = 10, seed = 1)
  expect_output(summary(result), "Rarefaction Summary")
  expect_output(summary(result), "95%")
})

test_that("rarefaction_taxonomic backward compatible", {
  d <- make_test_data()
  result <- rarefaction_taxonomic(d$comm, d$tax,
                                   index = "shannon",
                                   steps = 5, n_boot = 10, seed = 1)
  expect_true(is.data.frame(result))
  expect_true("sample_size" %in% names(result))
  expect_true("mean" %in% names(result))
  expect_equal(attr(result, "index"), "shannon")
})
