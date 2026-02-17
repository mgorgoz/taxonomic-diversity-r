test_that("deng_entropy_level returns Shannon entropy at species level", {
  # Equal abundances: H = ln(S)
  result <- deng_entropy_level(c(10, 10, 10, 10))
  expect_equal(result, log(4), tolerance = 1e-10)

  # Single species: H = 0
  expect_equal(deng_entropy_level(c(10)), 0)

  # Two equal species: H = ln(2)
  expect_equal(deng_entropy_level(c(5, 5)), log(2), tolerance = 1e-10)
})

test_that("deng_entropy_level with group_sizes differs from Shannon", {
  # 3 groups with different sizes
  abund <- c(9, 3, 7)
  gs <- c(3, 2, 3)

  ed_deng <- deng_entropy_level(abund, group_sizes = gs)
  ed_shannon <- deng_entropy_level(abund)  # No group_sizes = Shannon

  # Deng should be higher than Shannon when |Fi| > 1
  expect_true(ed_deng > ed_shannon)
})

test_that("deng_entropy_level handles zeros", {
  result <- deng_entropy_level(c(10, 0, 5, 0, 8))
  expected <- deng_entropy_level(c(10, 5, 8))
  expect_equal(result, expected)
})

test_that("deng_entropy_level validates input", {
  expect_error(deng_entropy_level(c(-1, 5)), "non-negative")
  expect_error(deng_entropy_level("abc"), "numeric")
})
