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

test_that("ozkan_pto returns all four components", {
  comm <- c(sp1 = 4, sp2 = 2, sp3 = 3, sp4 = 1, sp5 = 2,
            sp6 = 3, sp7 = 2, sp8 = 2)
  tax <- data.frame(
    Species = paste0("sp", 1:8),
    Genus   = c("G1", "G1", "G1", "G2", "G2", "G3", "G3", "G3"),
    Family  = c("F1", "F1", "F1", "F1", "F1", "F1", "F1", "F1"),
    stringsAsFactors = FALSE
  )

  result <- ozkan_pto(comm, tax)

  expect_type(result, "list")
  expect_true("uTO" %in% names(result))
  expect_true("TO" %in% names(result))
  expect_true("uTO_plus" %in% names(result))
  expect_true("TO_plus" %in% names(result))
  expect_true("Ed_levels" %in% names(result))
})

test_that("ozkan_pto satisfies ordering constraints", {
  # From Ozkan (2018): TO+ >= TO, uTO+ >= uTO, TO >= uTO, TO+ >= uTO+
  comm <- c(sp1 = 4, sp2 = 2, sp3 = 3, sp4 = 1, sp5 = 2,
            sp6 = 3, sp7 = 2, sp8 = 2)
  tax <- data.frame(
    Species = paste0("sp", 1:8),
    Genus   = c("G1", "G1", "G1", "G2", "G2", "G3", "G3", "G3"),
    Family  = c("F1", "F1", "F1", "F1", "F1", "F1", "F1", "F1"),
    stringsAsFactors = FALSE
  )

  r <- ozkan_pto(comm, tax)

  expect_true(r$TO_plus >= r$TO)
  expect_true(r$uTO_plus >= r$uTO)
  expect_true(r$TO >= r$uTO)
  expect_true(r$TO_plus >= r$uTO_plus)
})

test_that("ozkan_pto returns zero for single species", {
  comm <- c(sp1 = 10)
  tax <- data.frame(
    Species = "sp1",
    Genus = "G1",
    stringsAsFactors = FALSE
  )

  r <- ozkan_pto(comm, tax)
  expect_equal(r$uTO, 0)
  expect_equal(r$TO, 0)
})

test_that("ozkan_pto validates input", {
  tax <- data.frame(
    Species = c("sp1", "sp2"),
    Genus = c("G1", "G2"),
    stringsAsFactors = FALSE
  )

  expect_error(ozkan_pto(c(10, 5), tax), "named vector")
  expect_error(ozkan_pto(c(sp1 = -1, sp2 = 5), tax), "non-negative")
  expect_error(ozkan_pto(c(sp_x = 10), tax), "not found")
})

test_that("pto_components returns named numeric vector", {
  comm <- c(sp1 = 4, sp2 = 2, sp3 = 3)
  tax <- data.frame(
    Species = paste0("sp", 1:3),
    Genus   = c("G1", "G1", "G2"),
    Family  = c("F1", "F1", "F1"),
    stringsAsFactors = FALSE
  )

  result <- pto_components(comm, tax)
  expect_type(result, "double")
  expect_length(result, 4)
  expect_equal(names(result), c("uTO", "TO", "uTO_plus", "TO_plus"))
})
