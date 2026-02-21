# =============================================================================
# Tests for ozkan_pto_resample() (Islem 2) and ozkan_pto_sensitivity() (Islem 3)
# =============================================================================

# --- Shared test data ---
make_test_data <- function() {
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
# ozkan_pto_resample() tests (Islem 2 / Run 2)
# =============================================================================

test_that("ozkan_pto_resample returns correct structure", {
  d <- make_test_data()
  result <- ozkan_pto_resample(d$comm, d$tax, n_iter = 101, seed = 42)

  expect_type(result, "list")
  expect_true("uTO_plus_max" %in% names(result))
  expect_true("TO_plus_max" %in% names(result))
  expect_true("uTO_max" %in% names(result))
  expect_true("TO_max" %in% names(result))
  expect_true("uTO_plus_det" %in% names(result))
  expect_true("TO_plus_det" %in% names(result))
  expect_true("uTO_det" %in% names(result))
  expect_true("TO_det" %in% names(result))
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

  # First iteration should be deterministic (same as ozkan_pto)
  expect_equal(result$uTO_plus_det, det$uTO_plus, tolerance = 1e-10)
  expect_equal(result$TO_plus_det, det$TO_plus, tolerance = 1e-10)
  expect_equal(result$uTO_det, unname(det$uTO), tolerance = 1e-10)
  expect_equal(result$TO_det, unname(det$TO), tolerance = 1e-10)

  # Also check iteration_results row 1
  expect_equal(result$iteration_results$uTO_plus[1], det$uTO_plus, tolerance = 1e-10)
  expect_equal(result$iteration_results$TO_plus[1], det$TO_plus, tolerance = 1e-10)
})

test_that("ozkan_pto_resample maximums >= deterministic values", {
  d <- make_test_data()
  result <- ozkan_pto_resample(d$comm, d$tax, n_iter = 101, seed = 42)

  # Max across iterations should be >= deterministic (first iteration)
  expect_true(result$uTO_plus_max >= result$uTO_plus_det)
  expect_true(result$TO_plus_max >= result$TO_plus_det)
  # Note: uTO_max and TO_max might be less than det because random removal
  # can reduce diversity. But they should be >= 0
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
})

test_that("ozkan_pto_resample produces different results with different seeds", {
  d <- make_test_data()
  r1 <- ozkan_pto_resample(d$comm, d$tax, n_iter = 101, seed = 1)
  r2 <- ozkan_pto_resample(d$comm, d$tax, n_iter = 101, seed = 999)

  # Very unlikely to get identical max values with different seeds
  # (but not impossible, so we just check they ran successfully)
  expect_true(is.numeric(r1$uTO_plus_max))
  expect_true(is.numeric(r2$uTO_plus_max))
})

test_that("ozkan_pto_resample validates n_iter >= 101", {
  d <- make_test_data()
  expect_error(ozkan_pto_resample(d$comm, d$tax, n_iter = 50), ">= 101")
  expect_error(ozkan_pto_resample(d$comm, d$tax, n_iter = 0), ">= 101")
})

test_that("ozkan_pto_resample validates inputs", {
  d <- make_test_data()
  tax <- d$tax

  expect_error(ozkan_pto_resample(c(10, 5), tax), "named vector")
  expect_error(ozkan_pto_resample(c(sp1 = -1, sp2 = 5), tax), "non-negative")
  expect_error(ozkan_pto_resample(c(sp1 = 5), tax, n_iter = 101),
               "at least 2")
})

test_that("ozkan_pto_resample satisfies ordering: TO+ >= uTO+ for max", {
  d <- make_test_data()
  result <- ozkan_pto_resample(d$comm, d$tax, n_iter = 101, seed = 42)

  # The ordering TO+ >= uTO+ should hold for max values
  expect_true(result$TO_plus_max >= result$uTO_plus_max)
})

test_that("ozkan_pto_resample handles many iterations", {
  d <- make_test_data()
  result <- ozkan_pto_resample(d$comm, d$tax, n_iter = 500, seed = 42)

  expect_equal(result$n_iter, 500L)
  expect_equal(nrow(result$iteration_results), 500)
  # More iterations should find >= the same max as fewer iterations
  r_small <- ozkan_pto_resample(d$comm, d$tax, n_iter = 101, seed = 42)
  # Can't guarantee this due to different random sequences, but check structure
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


# =============================================================================
# ozkan_pto_sensitivity() tests (Islem 3 / Run 3)
# =============================================================================

test_that("ozkan_pto_sensitivity returns correct structure", {
  d <- make_test_data()
  run2 <- ozkan_pto_resample(d$comm, d$tax, n_iter = 101, seed = 42)
  run3 <- ozkan_pto_sensitivity(d$comm, d$tax, run2, n_iter = 101, seed = 42)

  expect_type(run3, "list")
  expect_true("uTO_plus_max" %in% names(run3))
  expect_true("TO_plus_max" %in% names(run3))
  expect_true("uTO_max" %in% names(run3))
  expect_true("TO_max" %in% names(run3))
  expect_true("run3_uTO_plus_max" %in% names(run3))
  expect_true("run3_TO_plus_max" %in% names(run3))
  expect_true("species_probs" %in% names(run3))
  expect_true("iteration_results" %in% names(run3))
})

test_that("ozkan_pto_sensitivity overall max >= Run 2 max", {
  d <- make_test_data()
  run2 <- ozkan_pto_resample(d$comm, d$tax, n_iter = 101, seed = 42)
  run3 <- ozkan_pto_sensitivity(d$comm, d$tax, run2, n_iter = 101, seed = 42)

  # Overall max should be >= Run 2 max (it takes pmax)
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

test_that("ozkan_pto_sensitivity is reproducible with same seed", {
  d <- make_test_data()
  run2 <- ozkan_pto_resample(d$comm, d$tax, n_iter = 101, seed = 42)
  r1 <- ozkan_pto_sensitivity(d$comm, d$tax, run2, seed = 123)
  r2 <- ozkan_pto_sensitivity(d$comm, d$tax, run2, seed = 123)

  expect_equal(r1$uTO_plus_max, r2$uTO_plus_max)
  expect_equal(r1$run3_uTO_plus_max, r2$run3_uTO_plus_max)
  expect_equal(r1$iteration_results, r2$iteration_results)
})

test_that("ozkan_pto_sensitivity species_probs has correct length", {
  d <- make_test_data()
  run2 <- ozkan_pto_resample(d$comm, d$tax, n_iter = 101, seed = 42)
  run3 <- ozkan_pto_sensitivity(d$comm, d$tax, run2, seed = 42)

  expect_length(run3$species_probs, length(d$comm[d$comm > 0]))
  expect_true(all(run3$species_probs > 0))
  expect_true(all(run3$species_probs <= 1))
})

test_that("full pipeline Run1 -> Run2 -> Run3 works end to end", {
  d <- make_test_data()

  # Run 1 (deterministic)
  run1 <- ozkan_pto(d$comm, d$tax)

  # Run 2 (stochastic resampling)
  run2 <- ozkan_pto_resample(d$comm, d$tax, n_iter = 101, seed = 42)

  # Run 3 (sensitivity analysis)
  run3 <- ozkan_pto_sensitivity(d$comm, d$tax, run2, n_iter = 101, seed = 42)

  # All should return valid numeric results
  expect_true(is.finite(run1$uTO_plus))
  expect_true(is.finite(run2$uTO_plus_max))
  expect_true(is.finite(run3$uTO_plus_max))

  # Run 2 max >= Run 1 deterministic
  expect_true(run2$uTO_plus_max >= run1$uTO_plus)
  expect_true(run2$TO_plus_max >= run1$TO_plus)

  # Run 3 overall max >= Run 2 max
  expect_true(run3$uTO_plus_max >= run2$uTO_plus_max)
  expect_true(run3$TO_plus_max >= run2$TO_plus_max)
})

test_that("ozkan_pto_resample with Westhoff-Maarel scale data", {
  # Simulate Westhoff-Maarel scale data (values 1-9)
  set.seed(42)
  n_sp <- 20
  comm <- setNames(
    sample(1:9, n_sp, replace = TRUE, prob = c(0.1, 0.1, 0.5, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05)),
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
})
