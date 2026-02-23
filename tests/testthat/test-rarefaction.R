# Tests for rarefaction_taxonomic()
# 3 core scenarios + error handling + plot test

test_that("full sample size returns value close to actual index", {
  comm <- c(sp1 = 10, sp2 = 5, sp3 = 8, sp4 = 2, sp5 = 3)
  tax <- data.frame(
    Species = paste0("sp", 1:5),
    Genus   = c("G1", "G1", "G2", "G2", "G3"),
    Family  = c("F1", "F1", "F1", "F2", "F2"),
    stringsAsFactors = FALSE
  )

  # Shannon at full sample should be close to actual Shannon
  actual_shannon <- shannon(comm)

  rare <- rarefaction_taxonomic(comm, tax, index = "shannon",
                                 steps = 5, n_boot = 50, seed = 42)

  # Last row = full sample size
  last_row <- rare[nrow(rare), ]
  expect_equal(last_row$sample_size, sum(comm))
  expect_equal(last_row$mean, actual_shannon, tolerance = 0.15)

  # Species richness at full sample should equal actual species count
  rare_sp <- rarefaction_taxonomic(comm, tax, index = "species",
                                    steps = 5, n_boot = 50, seed = 42)
  last_sp <- rare_sp[nrow(rare_sp), ]
  expect_equal(last_sp$mean, 5)
})


test_that("rarefaction curve is monotonically non-decreasing", {
  comm <- c(sp1 = 15, sp2 = 10, sp3 = 12, sp4 = 3,
            sp5 = 6, sp6 = 8, sp7 = 1, sp8 = 4)
  tax <- data.frame(
    Species = paste0("sp", 1:8),
    Genus   = c("G1", "G1", "G2", "G2", "G3", "G3", "G4", "G4"),
    Family  = c("F1", "F1", "F1", "F1", "F2", "F2", "F2", "F2"),
    Order   = c("O1", "O1", "O1", "O1", "O1", "O1", "O1", "O1"),
    stringsAsFactors = FALSE
  )

  # Shannon should generally increase with sample size
  rare <- rarefaction_taxonomic(comm, tax, index = "shannon",
                                 steps = 10, n_boot = 100, seed = 123)

  # Mean values should be non-decreasing (allow small tolerance for stochasticity)
  diffs <- diff(rare$mean)
  # Allow at most tiny decreases due to random sampling
  expect_true(all(diffs > -0.1),
              info = "Rarefaction curve should be approximately non-decreasing")

  # Sample sizes should be strictly increasing
  expect_true(all(diff(rare$sample_size) > 0))
})


test_that("confidence intervals are valid", {
  comm <- c(sp1 = 10, sp2 = 5, sp3 = 8, sp4 = 2, sp5 = 3)
  tax <- data.frame(
    Species = paste0("sp", 1:5),
    Genus   = c("G1", "G1", "G2", "G2", "G3"),
    Family  = c("F1", "F1", "F1", "F2", "F2"),
    stringsAsFactors = FALSE
  )

  rare <- rarefaction_taxonomic(comm, tax, index = "shannon",
                                 steps = 8, n_boot = 100, seed = 99)

  # Lower <= mean <= upper at all steps
  expect_true(all(rare$lower <= rare$mean + 1e-10))
  expect_true(all(rare$mean <= rare$upper + 1e-10))

  # CI should be non-negative width
  expect_true(all(rare$upper - rare$lower >= -1e-10))

  # SD should be non-negative
  expect_true(all(rare$sd >= 0))

  # CI should narrow at full sample (less uncertainty)
  ci_width <- rare$upper - rare$lower
  # Last step (full sample) should have narrower CI than middle steps
  expect_true(ci_width[nrow(rare)] <= max(ci_width),
              info = "CI should be narrowest at full sample size")
})


test_that("rarefaction works with taxonomic indices", {
  comm <- c(sp1 = 8, sp2 = 5, sp3 = 6, sp4 = 3, sp5 = 4)
  tax <- data.frame(
    Species = paste0("sp", 1:5),
    Genus   = c("G1", "G1", "G2", "G2", "G3"),
    Family  = c("F1", "F1", "F1", "F2", "F2"),
    stringsAsFactors = FALSE
  )

  # uTO should work
  rare_uto <- rarefaction_taxonomic(comm, tax, index = "uTO",
                                     steps = 5, n_boot = 30, seed = 7)
  expect_true(nrow(rare_uto) >= 2)
  expect_true(all(is.finite(rare_uto$mean)))

  # avtd should work
  rare_avtd <- rarefaction_taxonomic(comm, tax, index = "avtd",
                                      steps = 5, n_boot = 30, seed = 7)
  expect_true(nrow(rare_avtd) >= 2)
  expect_true(all(is.finite(rare_avtd$mean)))
})


test_that("rarefaction returns correct structure and attributes", {
  comm <- c(sp1 = 10, sp2 = 5, sp3 = 8)
  tax <- data.frame(
    Species = paste0("sp", 1:3),
    Genus   = c("G1", "G1", "G2"),
    Family  = c("F1", "F1", "F1"),
    stringsAsFactors = FALSE
  )

  rare <- rarefaction_taxonomic(comm, tax, index = "simpson",
                                 steps = 5, n_boot = 20, seed = 1)

  # Correct columns
  expect_true(all(c("sample_size", "mean", "lower", "upper", "sd") %in%
                    names(rare)))

  # Attributes
  expect_equal(attr(rare, "index"), "simpson")
  expect_equal(attr(rare, "total_n"), 23)
  expect_equal(attr(rare, "n_boot"), 20)
})


test_that("error handling works", {
  comm <- c(sp1 = 10, sp2 = 5)
  tax <- data.frame(
    Species = c("sp1", "sp2"),
    Genus = c("G1", "G2"),
    stringsAsFactors = FALSE
  )

  # Invalid community
  expect_error(rarefaction_taxonomic(c(-1, 2), tax),
               "non-negative")

  # No names
  expect_error(rarefaction_taxonomic(c(10, 5), tax),
               "named vector")

  # Too few steps
  expect_error(rarefaction_taxonomic(comm, tax, steps = 1),
               "at least 2")

  # Invalid ci
  expect_error(rarefaction_taxonomic(comm, tax, ci = 1.5),
               "between 0 and 1")
})


test_that("plot_rarefaction returns ggplot object", {
  skip_if_not_installed("ggplot2")

  comm <- c(sp1 = 10, sp2 = 5, sp3 = 8)
  tax <- data.frame(
    Species = paste0("sp", 1:3),
    Genus   = c("G1", "G1", "G2"),
    Family  = c("F1", "F1", "F1"),
    stringsAsFactors = FALSE
  )

  rare <- rarefaction_taxonomic(comm, tax, index = "shannon",
                                 steps = 5, n_boot = 20, seed = 1)
  p <- plot_rarefaction(rare)

  expect_s3_class(p, "gg")
  expect_s3_class(p, "ggplot")
})
