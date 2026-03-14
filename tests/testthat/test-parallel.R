# =============================================================================
# test-parallel.R
# Tests for parallel computing support in simulate_td, rarefaction_taxonomic,
# and batch_analysis
# =============================================================================

# --- Shared test data ---
tax_10sp <- data.frame(
  Species = paste0("sp", 1:10),
  Genus   = rep(c("G1", "G2", "G3", "G4", "G5"), each = 2),
  Family  = rep(c("F1", "F1", "F2", "F2", "F3"), each = 2),
  Order   = rep(c("O1", "O1", "O2", "O2", "O2"), each = 2),
  stringsAsFactors = FALSE
)

comm_5sp <- c(sp1 = 10, sp2 = 5, sp3 = 8, sp4 = 2, sp5 = 3)
tax_5sp <- data.frame(
  Species = paste0("sp", 1:5),
  Genus   = c("G1", "G1", "G2", "G2", "G3"),
  Family  = c("F1", "F1", "F1", "F2", "F2"),
  stringsAsFactors = FALSE
)

df_multi <- data.frame(
  Site      = c("A", "A", "A", "A", "B", "B", "B", "B"),
  Species   = c("sp1", "sp2", "sp3", "sp4", "sp1", "sp2", "sp3", "sp4"),
  Genus     = c("G1", "G1", "G2", "G2", "G1", "G1", "G2", "G2"),
  Family    = c("F1", "F1", "F1", "F2", "F1", "F1", "F1", "F2"),
  Order     = c("O1", "O1", "O1", "O1", "O1", "O1", "O1", "O1"),
  Abundance = c(10, 20, 15, 5, 5, 5, 5, 5),
  stringsAsFactors = FALSE
)


# ---- simulate_td: parallel = FALSE (default) still works ----

test_that("simulate_td: parallel = FALSE works as before", {
  sim <- simulate_td(tax_10sp, n_sim = 50, seed = 1, parallel = FALSE)
  expect_s3_class(sim, "td_simulation")
  expect_equal(nrow(sim), 9)
})


# ---- simulate_td: parallel = TRUE produces valid results ----

test_that("simulate_td: parallel = TRUE returns valid structure", {
  sim <- simulate_td(tax_10sp, n_sim = 50, parallel = TRUE, n_cores = 2)

  expect_s3_class(sim, "td_simulation")
  expect_s3_class(sim, "data.frame")
  expect_equal(nrow(sim), 9)
  expect_true(all(c("s", "mean_avtd", "lower_avtd", "upper_avtd",
                     "mean_vartd", "lower_vartd", "upper_vartd") %in%
                    names(sim)))
  # Values should be reasonable
  expect_true(all(sim$lower_avtd <= sim$mean_avtd))
  expect_true(all(sim$mean_avtd <= sim$upper_avtd))
})


test_that("simulate_td: parallel results are structurally equivalent to sequential", {
  seq_sim <- simulate_td(tax_10sp, n_sim = 50, seed = 42, parallel = FALSE)
  par_sim <- simulate_td(tax_10sp, n_sim = 50, parallel = TRUE, n_cores = 2)

  # Same columns and same number of rows
  expect_equal(names(seq_sim), names(par_sim))
  expect_equal(nrow(seq_sim), nrow(par_sim))
  expect_equal(seq_sim$s, par_sim$s)

  # Values should be in same ballpark (not identical due to different RNG streams)
  expect_true(all(is.finite(par_sim$mean_avtd)))
  expect_true(all(is.finite(par_sim$mean_vartd)))
})


test_that("simulate_td: parallel with n_cores = 1 works", {
  sim <- simulate_td(tax_10sp, n_sim = 30, parallel = TRUE, n_cores = 1)
  expect_s3_class(sim, "td_simulation")
  expect_equal(nrow(sim), 9)
})


# ---- rarefaction_taxonomic: parallel = FALSE (default) still works ----

test_that("rarefaction_taxonomic: parallel = FALSE works as before", {
  rare <- rarefaction_taxonomic(comm_5sp, tax_5sp, index = "shannon",
                                 steps = 5, n_boot = 20, seed = 1,
                                 parallel = FALSE)
  expect_s3_class(rare, "rarefaction_taxonomic")
  expect_true(nrow(rare) >= 2)
})


# ---- rarefaction_taxonomic: parallel = TRUE produces valid results ----

test_that("rarefaction_taxonomic: parallel = TRUE returns valid structure", {
  rare <- rarefaction_taxonomic(comm_5sp, tax_5sp, index = "shannon",
                                 steps = 5, n_boot = 20,
                                 parallel = TRUE, n_cores = 2)

  expect_s3_class(rare, "rarefaction_taxonomic")
  expect_s3_class(rare, "data.frame")
  expect_true(all(c("sample_size", "mean", "lower", "upper", "sd") %in%
                    names(rare)))
  expect_true(all(rare$lower <= rare$mean + 1e-10))
  expect_true(all(rare$mean <= rare$upper + 1e-10))
})


test_that("rarefaction_taxonomic: parallel results are structurally equivalent", {
  seq_rare <- rarefaction_taxonomic(comm_5sp, tax_5sp, index = "shannon",
                                     steps = 5, n_boot = 30, seed = 42,
                                     parallel = FALSE)
  par_rare <- rarefaction_taxonomic(comm_5sp, tax_5sp, index = "shannon",
                                     steps = 5, n_boot = 30,
                                     parallel = TRUE, n_cores = 2)

  expect_equal(names(seq_rare), names(par_rare))
  expect_equal(nrow(seq_rare), nrow(par_rare))
  expect_equal(seq_rare$sample_size, par_rare$sample_size)
  expect_true(all(is.finite(par_rare$mean)))
})


test_that("rarefaction_taxonomic: parallel with n_cores = 1 works", {
  rare <- rarefaction_taxonomic(comm_5sp, tax_5sp, index = "shannon",
                                 steps = 3, n_boot = 10,
                                 parallel = TRUE, n_cores = 1)
  expect_s3_class(rare, "rarefaction_taxonomic")
})


# ---- batch_analysis: parallel = FALSE (default) still works ----

test_that("batch_analysis: parallel = FALSE works as before", {
  result <- batch_analysis(df_multi, parallel = FALSE)
  expect_s3_class(result, "batch_analysis")
  expect_equal(nrow(result), 2)
})


# ---- batch_analysis: parallel = TRUE produces valid results ----

test_that("batch_analysis: parallel = TRUE returns valid structure", {
  result <- batch_analysis(df_multi, parallel = TRUE, n_cores = 2)

  expect_s3_class(result, "batch_analysis")
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 2)
  expect_equal(result$Site, c("A", "B"))
  expect_true(all(is.numeric(result$Shannon)))
})


test_that("batch_analysis: parallel results match sequential results", {
  seq_result <- batch_analysis(df_multi, parallel = FALSE)
  par_result <- batch_analysis(df_multi, parallel = TRUE, n_cores = 2)

  # Results should be identical (no randomness in batch_analysis)
  expect_equal(seq_result$Shannon, par_result$Shannon)
  expect_equal(seq_result$Simpson, par_result$Simpson)
  expect_equal(seq_result$Delta, par_result$Delta)
  expect_equal(seq_result$AvTD, par_result$AvTD)
  expect_equal(seq_result$uTO, par_result$uTO)
})


test_that("batch_analysis: parallel with single site falls back to sequential", {
  df_single <- data.frame(
    Species   = c("sp1", "sp2", "sp3"),
    Genus     = c("G1", "G1", "G2"),
    Family    = c("F1", "F1", "F1"),
    Order     = c("O1", "O1", "O1"),
    Abundance = c(10, 20, 15),
    stringsAsFactors = FALSE
  )
  # parallel = TRUE with 1 site should still work (sequential fallback)
  result <- batch_analysis(df_single, parallel = TRUE, n_cores = 2)
  expect_s3_class(result, "batch_analysis")
  expect_equal(nrow(result), 1)
})
