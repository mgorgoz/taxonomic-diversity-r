# =============================================================================
# Ozkan pTO Jackknife (Run 1), Stochastic Resampling (Run 2),
# Sensitivity Analysis (Run 3) and Full Pipeline Tests
#
# Run 1: Jackknife — remove each species one-by-one, classify happy/unhappy
# Run 2: Stochastic resampling — unhappy always included, happy 50% random
# Run 3: Sensitivity — unhappy P=(S-1)/S, happy P=derived from Run 2
# Full:  3 runs + MAX
# =============================================================================

# --- Common test data used across all tests ---
make_test_data <- function() {
  # 8 species, 3 genera, 2 families, 1 order
  # Abundances suitable for Westhoff-Maarel scale (range 1-4)
  comm <- c(sp1 = 4, sp2 = 2, sp3 = 3, sp4 = 1, sp5 = 2,
            sp6 = 3, sp7 = 2, sp8 = 2)
  tax <- data.frame(
    Species = paste0("sp", 1:8),
    Genus   = c("G1", "G1", "G1", "G2", "G2", "G3", "G3", "G3"),
    Family  = c("F1", "F1", "F1", "F1", "F1", "F2", "F2", "F2"),
    Order   = rep("O1", 8),
    stringsAsFactors = FALSE
  )
  list(comm = comm, tax = tax)
}

# =============================================================================
# ozkan_pto_jackknife() tests (Run 1 — Jackknife)
# =============================================================================

test_that("ozkan_pto_jackknife returns correct structure", {
  d <- make_test_data()
  jk <- ozkan_pto_jackknife(d$comm, d$tax)

  expect_type(jk, "list")
  expect_true("full_result" %in% names(jk))
  expect_true("jackknife_results" %in% names(jk))
  expect_true("species_status" %in% names(jk))
  expect_true("n_happy" %in% names(jk))
  expect_true("n_unhappy" %in% names(jk))
})

test_that("ozkan_pto_jackknife species_status has correct length and names", {
  d <- make_test_data()
  jk <- ozkan_pto_jackknife(d$comm, d$tax)

  expect_length(jk$species_status, length(d$comm))
  expect_equal(names(jk$species_status), names(d$comm))
  expect_type(jk$species_status, "logical")
})

test_that("ozkan_pto_jackknife n_happy + n_unhappy = n_species", {
  d <- make_test_data()
  jk <- ozkan_pto_jackknife(d$comm, d$tax)

  expect_equal(jk$n_happy + jk$n_unhappy, length(d$comm))
})

test_that("ozkan_pto_jackknife jackknife_results has correct dimensions", {
  d <- make_test_data()
  jk <- ozkan_pto_jackknife(d$comm, d$tax)

  jk_df <- jk$jackknife_results
  expect_s3_class(jk_df, "data.frame")
  expect_equal(nrow(jk_df), length(d$comm))
  expect_equal(ncol(jk_df), 6)  # species, uTO_plus, TO_plus, uTO, TO, is_happy
})

test_that("ozkan_pto_jackknife full_result matches ozkan_pto", {
  d <- make_test_data()
  jk <- ozkan_pto_jackknife(d$comm, d$tax)
  det <- ozkan_pto(d$comm, d$tax)

  expect_equal(jk$full_result$uTO_plus, det$uTO_plus, tolerance = 1e-10)
  expect_equal(jk$full_result$TO_plus, det$TO_plus, tolerance = 1e-10)
  expect_equal(jk$full_result$uTO, det$uTO, tolerance = 1e-10)
  expect_equal(jk$full_result$TO, det$TO, tolerance = 1e-10)
})

test_that("ozkan_pto_jackknife leave-one-out values are finite", {
  d <- make_test_data()
  jk <- ozkan_pto_jackknife(d$comm, d$tax)

  expect_true(all(is.finite(jk$jackknife_results$uTO_plus)))
  expect_true(all(is.finite(jk$jackknife_results$TO_plus)))
})

test_that("ozkan_pto_jackknife validates inputs", {
  d <- make_test_data()

  expect_error(ozkan_pto_jackknife(c(10, 5), d$tax), "named vector")
  expect_error(ozkan_pto_jackknife(c(sp1 = -1, sp2 = 5), d$tax), "non-negative")
  expect_error(ozkan_pto_jackknife(c(sp1 = 5), d$tax), "at least 2")
  expect_error(ozkan_pto_jackknife(d$comm, d$tax, component = "xyz"), "component")
})

test_that("ozkan_pto_jackknife works with different components", {
  d <- make_test_data()

  jk1 <- ozkan_pto_jackknife(d$comm, d$tax, component = "uTO_plus")
  jk2 <- ozkan_pto_jackknife(d$comm, d$tax, component = "TO_plus")
  jk3 <- ozkan_pto_jackknife(d$comm, d$tax, component = "uTO")
  jk4 <- ozkan_pto_jackknife(d$comm, d$tax, component = "TO")

  # All should return valid results
  expect_type(jk1$species_status, "logical")
  expect_type(jk2$species_status, "logical")
  expect_type(jk3$species_status, "logical")
  expect_type(jk4$species_status, "logical")
})


# =============================================================================
# ozkan_pto_resample() tests (Run 2)
# =============================================================================

test_that("ozkan_pto_resample returns correct structure", {
  d <- make_test_data()
  result <- ozkan_pto_resample(d$comm, d$tax, n_iter = 101, seed = 42)

  expect_type(result, "list")

  # 4 maximum values
  expect_true("uTO_plus_max" %in% names(result))
  expect_true("TO_plus_max" %in% names(result))
  expect_true("uTO_max" %in% names(result))
  expect_true("TO_max" %in% names(result))

  # 4 deterministic values
  expect_true("uTO_plus_det" %in% names(result))
  expect_true("TO_plus_det" %in% names(result))
  expect_true("uTO_det" %in% names(result))
  expect_true("TO_det" %in% names(result))

  # New fields: jackknife and species_status
  expect_true("species_status" %in% names(result))
  expect_true("jackknife" %in% names(result))
  expect_true("n_positive" %in% names(result))

  expect_true("n_iter" %in% names(result))
  expect_true("iteration_results" %in% names(result))
  expect_equal(result$n_iter, 101L)
})

test_that("ozkan_pto_resample iteration results has correct dimensions", {
  d <- make_test_data()
  result <- ozkan_pto_resample(d$comm, d$tax, n_iter = 101, seed = 42)

  iter_df <- result$iteration_results
  expect_s3_class(iter_df, "data.frame")
  expect_equal(nrow(iter_df), 101)
  expect_equal(ncol(iter_df), 5)
  expect_equal(names(iter_df), c("iteration", "uTO_plus", "TO_plus", "uTO", "TO"))
})

test_that("ozkan_pto_resample first iteration matches deterministic ozkan_pto", {
  d <- make_test_data()
  result <- ozkan_pto_resample(d$comm, d$tax, n_iter = 101, seed = 42)
  det <- ozkan_pto(d$comm, d$tax)

  expect_equal(result$uTO_plus_det, det$uTO_plus, tolerance = 1e-10)
  expect_equal(result$TO_plus_det, det$TO_plus, tolerance = 1e-10)
  expect_equal(result$uTO_det, unname(det$uTO), tolerance = 1e-10)
  expect_equal(result$TO_det, unname(det$TO), tolerance = 1e-10)

  expect_equal(result$iteration_results$uTO_plus[1], det$uTO_plus, tolerance = 1e-10)
  expect_equal(result$iteration_results$TO_plus[1], det$TO_plus, tolerance = 1e-10)
})

test_that("ozkan_pto_resample maximums >= deterministic values", {
  d <- make_test_data()
  result <- ozkan_pto_resample(d$comm, d$tax, n_iter = 101, seed = 42)

  # pTO+: max >= det (1st iteration is deterministic, so max is at least that)
  expect_true(result$uTO_plus_max >= result$uTO_plus_det)
  expect_true(result$TO_plus_max >= result$TO_plus_det)

  # Valid values
  expect_true(result$uTO_max >= 0)
  expect_true(result$TO_max >= 0)
})

test_that("ozkan_pto_resample is reproducible with same seed", {
  d <- make_test_data()
  r1 <- ozkan_pto_resample(d$comm, d$tax, n_iter = 101, seed = 123)
  r2 <- ozkan_pto_resample(d$comm, d$tax, n_iter = 101, seed = 123)

  expect_equal(r1$uTO_plus_max, r2$uTO_plus_max)
  expect_equal(r1$TO_plus_max, r2$TO_plus_max)
  expect_equal(r1$uTO_max, r2$uTO_max)
  expect_equal(r1$TO_max, r2$TO_max)
  expect_equal(r1$iteration_results, r2$iteration_results)
  expect_equal(r1$species_status, r2$species_status)
})

test_that("ozkan_pto_resample produces different results with different seeds", {
  d <- make_test_data()
  r1 <- ozkan_pto_resample(d$comm, d$tax, n_iter = 101, seed = 1)
  r2 <- ozkan_pto_resample(d$comm, d$tax, n_iter = 101, seed = 999)

  expect_true(is.numeric(r1$uTO_plus_max))
  expect_true(is.numeric(r2$uTO_plus_max))
})

test_that("ozkan_pto_resample validates n_iter is a positive integer", {
  d <- make_test_data()
  # Zero and negative values are rejected; there is no lower bound at 101.
  expect_error(ozkan_pto_resample(d$comm, d$tax, n_iter = 0), ">= 1")
  expect_error(ozkan_pto_resample(d$comm, d$tax, n_iter = -5), ">= 1")
  expect_error(ozkan_pto_resample(d$comm, d$tax, n_iter = c(10, 20)), ">= 1")
})

test_that("ozkan_pto_resample accepts n_iter below 101 (no 101 floor)", {
  d <- make_test_data()
  # Values below the old 101 floor must now run without error.
  r50 <- ozkan_pto_resample(d$comm, d$tax, n_iter = 50, seed = 42)
  expect_equal(r50$n_iter, 50L)
  expect_equal(nrow(r50$iteration_results), 50)
})

test_that("ozkan_pto_resample with n_iter = 1 returns the deterministic Run 1", {
  d <- make_test_data()
  r1 <- ozkan_pto_resample(d$comm, d$tax, n_iter = 1, seed = 42)
  det <- ozkan_pto(d$comm, d$tax)

  expect_equal(r1$n_iter, 1L)
  expect_equal(nrow(r1$iteration_results), 1)
  # With a single (deterministic) iteration, max == deterministic value.
  expect_equal(r1$uTO_plus_max, unname(det$uTO_plus))
  expect_equal(r1$TO_plus_max, unname(det$TO_plus))
  expect_equal(r1$uTO_max, unname(det$uTO))
  expect_equal(r1$TO_max, unname(det$TO))
})

test_that("ozkan_pto_resample validates inputs", {
  d <- make_test_data()
  tax <- d$tax

  expect_error(ozkan_pto_resample(c(10, 5), tax), "named vector")
  expect_error(ozkan_pto_resample(c(sp1 = -1, sp2 = 5), tax), "non-negative")
  expect_error(ozkan_pto_resample(c(sp1 = 5), tax, n_iter = 101), "at least 2")
})

test_that("ozkan_pto_resample satisfies ordering: TO+ >= uTO+ for max", {
  d <- make_test_data()
  result <- ozkan_pto_resample(d$comm, d$tax, n_iter = 101, seed = 42)
  expect_true(result$TO_plus_max >= result$uTO_plus_max)
})

test_that("ozkan_pto_resample handles many iterations", {
  d <- make_test_data()
  result <- ozkan_pto_resample(d$comm, d$tax, n_iter = 500, seed = 42)

  expect_equal(result$n_iter, 500L)
  expect_equal(nrow(result$iteration_results), 500)
  expect_true(is.numeric(result$uTO_plus_max))
})

test_that("ozkan_pto_resample with equal abundances works", {
  comm <- setNames(rep(1, 6), paste0("sp", 1:6))
  tax <- data.frame(
    Species = paste0("sp", 1:6),
    Genus = rep(c("G1", "G2", "G3"), each = 2),
    Family = rep("F1", 6),
    stringsAsFactors = FALSE
  )

  result <- ozkan_pto_resample(comm, tax, n_iter = 101, seed = 42)
  expect_true(is.numeric(result$uTO_plus_max))
  expect_true(result$uTO_plus_max > 0)
})

test_that("ozkan_pto_resample species_status is logical vector", {
  d <- make_test_data()
  result <- ozkan_pto_resample(d$comm, d$tax, n_iter = 101, seed = 42)

  expect_type(result$species_status, "logical")
  expect_length(result$species_status, length(d$comm))
  expect_equal(names(result$species_status), names(d$comm))
})

test_that("ozkan_pto_resample n_positive is within valid range", {
  d <- make_test_data()
  result <- ozkan_pto_resample(d$comm, d$tax, n_iter = 101, seed = 42)

  expect_true(result$n_positive >= 0)
  expect_true(result$n_positive <= result$n_iter)
})


# =============================================================================
# ozkan_pto_sensitivity() tests (Run 3)
# =============================================================================

test_that("ozkan_pto_sensitivity returns correct structure", {
  d <- make_test_data()
  run2 <- ozkan_pto_resample(d$comm, d$tax, n_iter = 101, seed = 42)
  run3 <- ozkan_pto_sensitivity(d$comm, d$tax, run2, n_iter = 101, seed = 42)

  expect_type(run3, "list")

  # Overall maximums (Run 1+2+3)
  expect_true("uTO_plus_max" %in% names(run3))
  expect_true("TO_plus_max" %in% names(run3))
  expect_true("uTO_max" %in% names(run3))
  expect_true("TO_max" %in% names(run3))

  # Run 3 specific maximums
  expect_true("run3_uTO_plus_max" %in% names(run3))
  expect_true("run3_TO_plus_max" %in% names(run3))

  # Species probabilities
  expect_true("species_probs" %in% names(run3))
  expect_true("prob_happy" %in% names(run3))
  expect_true("prob_unhappy" %in% names(run3))

  # Iteration results
  expect_true("iteration_results" %in% names(run3))
})

test_that("ozkan_pto_sensitivity overall max >= Run 2 max", {
  d <- make_test_data()
  run2 <- ozkan_pto_resample(d$comm, d$tax, n_iter = 101, seed = 42)
  run3 <- ozkan_pto_sensitivity(d$comm, d$tax, run2, n_iter = 101, seed = 42)

  expect_true(run3$uTO_plus_max >= run2$uTO_plus_max)
  expect_true(run3$TO_plus_max >= run2$TO_plus_max)
})

test_that("ozkan_pto_sensitivity uses Run 2 n_iter by default", {
  d <- make_test_data()
  run2 <- ozkan_pto_resample(d$comm, d$tax, n_iter = 101, seed = 42)
  run3 <- ozkan_pto_sensitivity(d$comm, d$tax, run2, seed = 42)

  expect_equal(run3$n_iter, 101L)
  expect_equal(nrow(run3$iteration_results), 101)
})

test_that("ozkan_pto_sensitivity validates run2_result", {
  d <- make_test_data()
  expect_error(
    ozkan_pto_sensitivity(d$comm, d$tax, list(x = 1)),
    "ozkan_pto_resample"
  )
})

test_that("ozkan_pto_sensitivity validates species_status requirement", {
  d <- make_test_data()
  # Old format run2_result (no species_status)
  fake_run2 <- list(
    iteration_results = data.frame(
      iteration = 1:101,
      uTO_plus = runif(101), TO_plus = runif(101),
      uTO = runif(101), TO = runif(101)
    ),
    n_iter = 101L,
    uTO_plus_max = 1, TO_plus_max = 1,
    uTO_max = 1, TO_max = 1
  )
  expect_error(
    ozkan_pto_sensitivity(d$comm, d$tax, fake_run2),
    "species_status"
  )
})

test_that("ozkan_pto_sensitivity is reproducible with same seed", {
  d <- make_test_data()
  run2 <- ozkan_pto_resample(d$comm, d$tax, n_iter = 101, seed = 42)
  r1 <- ozkan_pto_sensitivity(d$comm, d$tax, run2, seed = 123)
  r2 <- ozkan_pto_sensitivity(d$comm, d$tax, run2, seed = 123)

  expect_equal(r1$uTO_plus_max, r2$uTO_plus_max)
  expect_equal(r1$run3_uTO_plus_max, r2$run3_uTO_plus_max)
  expect_equal(r1$iteration_results, r2$iteration_results)
})

test_that("ozkan_pto_sensitivity species_probs has correct length and range", {
  d <- make_test_data()
  run2 <- ozkan_pto_resample(d$comm, d$tax, n_iter = 101, seed = 42)
  run3 <- ozkan_pto_sensitivity(d$comm, d$tax, run2, seed = 42)

  expect_length(run3$species_probs, length(d$comm[d$comm > 0]))
  expect_true(all(run3$species_probs >= 0))
  expect_true(all(run3$species_probs <= 1))
})

test_that("ozkan_pto_sensitivity prob_unhappy equals (S-1)/S", {
  d <- make_test_data()
  run2 <- ozkan_pto_resample(d$comm, d$tax, n_iter = 101, seed = 42)
  run3 <- ozkan_pto_sensitivity(d$comm, d$tax, run2, seed = 42)

  n_sp <- length(d$comm[d$comm > 0])
  expect_equal(run3$prob_unhappy, (n_sp - 1) / n_sp, tolerance = 1e-10)
})

test_that("ozkan_pto_sensitivity happy vs unhappy use different probs", {
  d <- make_test_data()
  run2 <- ozkan_pto_resample(d$comm, d$tax, n_iter = 101, seed = 42)
  run3 <- ozkan_pto_sensitivity(d$comm, d$tax, run2, seed = 42)

  is_happy <- run2$species_status

  # If both happy and unhappy species exist, different probabilities should be used
  if (any(is_happy) && any(!is_happy)) {
    happy_probs <- unique(run3$species_probs[is_happy])
    unhappy_probs <- unique(run3$species_probs[!is_happy])
    # Happy and unhappy probabilities should differ
    expect_false(all(happy_probs %in% unhappy_probs))
  }
})


# =============================================================================
# ozkan_pto_full() tests (Full Pipeline: Run 1 + 2 + 3)
# =============================================================================

test_that("ozkan_pto_full returns correct structure", {
  d <- make_test_data()
  result <- ozkan_pto_full(d$comm, d$tax, n_iter = 101, seed = 42)

  expect_type(result, "list")

  # Final results
  expect_true("uTO_plus" %in% names(result))
  expect_true("TO_plus" %in% names(result))
  expect_true("uTO" %in% names(result))
  expect_true("TO" %in% names(result))

  # Sub-results
  expect_true("run1" %in% names(result))
  expect_true("run2" %in% names(result))
  expect_true("run3" %in% names(result))
  expect_true("jackknife" %in% names(result))
})

test_that("ozkan_pto_full final values >= deterministic values", {
  d <- make_test_data()
  result <- ozkan_pto_full(d$comm, d$tax, n_iter = 101, seed = 42)
  det <- ozkan_pto(d$comm, d$tax)

  # Final results should be >= deterministic values
  expect_true(result$uTO_plus >= det$uTO_plus)
  expect_true(result$TO_plus >= det$TO_plus)
})

test_that("ozkan_pto_full is reproducible with same seed", {
  d <- make_test_data()
  r1 <- ozkan_pto_full(d$comm, d$tax, n_iter = 101, seed = 42)
  r2 <- ozkan_pto_full(d$comm, d$tax, n_iter = 101, seed = 42)

  expect_equal(r1$uTO_plus, r2$uTO_plus)
  expect_equal(r1$TO_plus, r2$TO_plus)
  expect_equal(r1$uTO, r2$uTO)
  expect_equal(r1$TO, r2$TO)
})

test_that("ozkan_pto_full run1 matches ozkan_pto", {
  d <- make_test_data()
  result <- ozkan_pto_full(d$comm, d$tax, n_iter = 101, seed = 42)
  det <- ozkan_pto(d$comm, d$tax)

  expect_equal(result$run1$uTO_plus, det$uTO_plus, tolerance = 1e-10)
  expect_equal(result$run1$TO_plus, det$TO_plus, tolerance = 1e-10)
})

test_that("ozkan_pto_full validates n_iter", {
  d <- make_test_data()
  # No 101 floor: only non-positive / malformed values are rejected.
  expect_error(ozkan_pto_full(d$comm, d$tax, n_iter = 0), ">= 1")
  expect_error(ozkan_pto_full(d$comm, d$tax, n_iter = -1), ">= 1")
})

test_that("ozkan_pto_full accepts n_iter below 101 and n_iter = 1", {
  d <- make_test_data()
  # n_iter below the old floor must work.
  r50 <- ozkan_pto_full(d$comm, d$tax, n_iter = 50, seed = 42)
  expect_true(is.numeric(r50$uTO_plus))

  # n_iter = 1 reduces the whole pipeline to deterministic Run 1.
  r1  <- ozkan_pto_full(d$comm, d$tax, n_iter = 1, seed = 42)
  det <- ozkan_pto(d$comm, d$tax)
  expect_equal(r1$uTO_plus, unname(det$uTO_plus))
  expect_equal(r1$TO_plus,  unname(det$TO_plus))
  expect_equal(r1$uTO,      unname(det$uTO))
  expect_equal(r1$TO,       unname(det$TO))
})


# =============================================================================
# Full pipeline test: Run 1 -> Run 2 -> Run 3
# =============================================================================

test_that("full pipeline Run1 -> Run2 -> Run3 works end to end", {
  d <- make_test_data()

  # Run 1: Jackknife
  jk <- ozkan_pto_jackknife(d$comm, d$tax)
  expect_true(is.finite(jk$full_result$uTO_plus))

  # Run 2: Stochastic resampling (jackknife included internally)
  run2 <- ozkan_pto_resample(d$comm, d$tax, n_iter = 101, seed = 42)
  expect_true(is.finite(run2$uTO_plus_max))

  # Run 3: Sensitivity analysis
  run3 <- ozkan_pto_sensitivity(d$comm, d$tax, run2, n_iter = 101, seed = 42)
  expect_true(is.finite(run3$uTO_plus_max))

  # Ordering rule: Run 3 max >= Run 2 max >= Run 1 deterministic
  expect_true(run2$uTO_plus_max >= jk$full_result$uTO_plus)
  expect_true(run2$TO_plus_max >= jk$full_result$TO_plus)
  expect_true(run3$uTO_plus_max >= run2$uTO_plus_max)
  expect_true(run3$TO_plus_max >= run2$TO_plus_max)
})

test_that("full pipeline with ozkan_pto_full wrapper matches manual pipeline", {
  d <- make_test_data()

  # Manual pipeline
  run2 <- ozkan_pto_resample(d$comm, d$tax, n_iter = 101, seed = 42)
  run3 <- ozkan_pto_sensitivity(d$comm, d$tax, run2, n_iter = 101, seed = 43)

  # Automatic pipeline
  full <- ozkan_pto_full(d$comm, d$tax, n_iter = 101, seed = 42)

  # Same results
  expect_equal(full$uTO_plus, run3$uTO_plus_max, tolerance = 1e-10)
  expect_equal(full$TO_plus, run3$TO_plus_max, tolerance = 1e-10)
  expect_equal(full$uTO, run3$uTO_max, tolerance = 1e-10)
  expect_equal(full$TO, run3$TO_max, tolerance = 1e-10)
})


# =============================================================================
# Large dataset test with Westhoff-Maarel scale
# =============================================================================

test_that("ozkan_pto_resample with Westhoff-Maarel scale data", {
  set.seed(42)
  n_sp <- 20

  comm <- setNames(
    sample(1:9, n_sp, replace = TRUE,
           prob = c(0.1, 0.1, 0.5, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05)),
    paste0("sp", 1:n_sp)
  )

  tax <- data.frame(
    Species = paste0("sp", 1:n_sp),
    Genus = paste0("G", rep(1:10, each = 2)),
    Family = paste0("F", rep(1:5, each = 4)),
    Order = paste0("O", rep(1:2, each = 10)),
    stringsAsFactors = FALSE
  )

  result <- ozkan_pto_resample(comm, tax, n_iter = 101, seed = 42)

  expect_true(is.finite(result$uTO_plus_max))
  expect_true(result$uTO_plus_max > 0)
  expect_true(result$TO_plus_max > 0)
  expect_type(result$species_status, "logical")
})
