test_that("shannon calculates correctly", {
  # Equal abundances: H' = ln(S)
  comm <- c(10, 10, 10, 10)
  expect_equal(shannon(comm), log(4), tolerance = 1e-10)

  # Single species: H' = 0
  expect_equal(shannon(c(10)), 0)

  # With base 2
  comm <- c(10, 10)
  expect_equal(shannon(comm, base = 2), 1, tolerance = 1e-10)
})

test_that("shannon handles edge cases", {
  # Zeros are removed
  comm <- c(10, 0, 5, 0, 8)
  expect_equal(shannon(comm), shannon(c(10, 5, 8)))

  # All zeros
  expect_equal(shannon(c(0, 0, 0)), 0)
})

test_that("shannon validates input", {
  expect_error(shannon(c(-1, 5)), "non-negative")
  expect_error(shannon("abc"), "numeric")
})

test_that("simpson calculates correctly", {
  # Equal abundances: D = 1/S, Gini = 1 - 1/S
  comm <- c(10, 10, 10, 10)
  expect_equal(simpson(comm, "dominance"), 0.25, tolerance = 1e-10)
  expect_equal(simpson(comm, "gini_simpson"), 0.75, tolerance = 1e-10)
  expect_equal(simpson(comm, "inverse"), 4, tolerance = 1e-10)

  # Single species: D = 1
  expect_equal(simpson(c(10), "dominance"), 1)
})

test_that("simpson validates input", {
  expect_error(simpson(c(-1, 5)), "non-negative")
})
