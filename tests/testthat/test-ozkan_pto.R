# =============================================================================
# Ozkan pTO (Run 1) Function Tests
# Validates the correctness of ozkan_pto() and pto_components() functions.
#
# These tests cover the deterministic (Run 1) computation from Ozkan (2018).
# 4 components are produced: uTO, TO, uTO+, TO+
# =============================================================================

test_that("ozkan_pto returns all four components", {
  # Sample community with 8 species, 3 genera, 1 family
  comm <- c(sp1 = 4, sp2 = 2, sp3 = 3, sp4 = 1, sp5 = 2,
            sp6 = 3, sp7 = 2, sp8 = 2)
  tax <- data.frame(
    Species = paste0("sp", 1:8),
    Genus   = c("G1", "G1", "G1", "G2", "G2", "G3", "G3", "G3"),
    Family  = c("F1", "F1", "F1", "F1", "F1", "F1", "F1", "F1"),
    stringsAsFactors = FALSE
  )

  result <- ozkan_pto(comm, tax)

  # Function should return a list
  expect_type(result, "list")

  # List should contain 10 components:
  # uTO, TO, uTO_plus, TO_plus = all levels
  # uTO_max, TO_max, uTO_plus_max, TO_plus_max = informative levels
  # Ed_levels, max_informative_level
  expect_true("uTO" %in% names(result))
  expect_true("TO" %in% names(result))
  expect_true("uTO_plus" %in% names(result))
  expect_true("TO_plus" %in% names(result))
  expect_true("uTO_max" %in% names(result))
  expect_true("TO_max" %in% names(result))
  expect_true("uTO_plus_max" %in% names(result))
  expect_true("TO_plus_max" %in% names(result))
  expect_true("Ed_levels" %in% names(result))
  expect_true("max_informative_level" %in% names(result))
})

test_that("ozkan_pto satisfies ordering constraints", {
  # According to Ozkan (2018), the following ordering should always hold:
  # TO+ >= TO    (distance >= diversity, weighted)
  # uTO+ >= uTO  (distance >= diversity, unweighted)
  # TO >= uTO    (weighted >= unweighted, diversity)
  # TO+ >= uTO+  (weighted >= unweighted, distance)
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
  # Diversity should be zero for a single-species community
  # Because there is no other species to compare with
  comm <- c(sp1 = 10)
  tax <- data.frame(
    Species = "sp1",
    Genus = "G1",
    stringsAsFactors = FALSE
  )

  r <- ozkan_pto(comm, tax)
  expect_equal(r$uTO, 0)
  expect_equal(r$TO, 0)
  expect_equal(r$uTO_max, 0)
  expect_equal(r$TO_max, 0)
  expect_equal(r$uTO_plus_max, 0)
  expect_equal(r$TO_plus_max, 0)
  expect_equal(r$max_informative_level, 0L)
})

test_that("ozkan_pto validates input", {
  # Function should provide meaningful error messages for invalid inputs
  tax <- data.frame(
    Species = c("sp1", "sp2"),
    Genus = c("G1", "G2"),
    stringsAsFactors = FALSE
  )

  # Unnamed vector: taxonomy cannot be matched without species names
  expect_error(ozkan_pto(c(10, 5), tax), "named vector")

  # Negative abundance: abundance cannot be negative in nature
  expect_error(ozkan_pto(c(sp1 = -1, sp2 = 5), tax), "non-negative")

  # Species name not in taxonomy: no match can be found
  expect_error(ozkan_pto(c(sp_x = 10), tax), "not found")
})

test_that("pto_components returns named numeric vector with 8 elements", {
  # pto_components() is a shortcut version of ozkan_pto()
  # Returns a named numeric vector with 8 elements (Run 1+2+3 from the Excel macro)
  comm <- c(sp1 = 4, sp2 = 2, sp3 = 3)
  tax <- data.frame(
    Species = paste0("sp", 1:3),
    Genus   = c("G1", "G1", "G2"),
    Family  = c("F1", "F1", "F1"),
    stringsAsFactors = FALSE
  )

  result <- pto_components(comm, tax)

  # Should be of numeric (double) type
  expect_type(result, "double")

  # Should have 8 elements
  expect_length(result, 8)

  # Names should be in the correct order
  expect_equal(names(result), c("uTO", "TO", "uTO_plus", "TO_plus",
                                "uTO_max", "TO_max", "uTO_plus_max", "TO_plus_max"))
})

test_that("ozkan_pto uses presence-based (equal weight) entropy at species level", {
  # In the pTO formula, each species is given EQUAL WEIGHT at the species level (1/S)
  # Abundance information is not used directly — it only determines
  # which species survive during the "slicing" stage
  #
  # Therefore, species-level Deng entropy should equal ln(S)
  # Regardless of abundances (100, 1, 1, 1, 1) -> Ed_S = ln(5)
  comm <- c(sp1 = 100, sp2 = 1, sp3 = 1, sp4 = 1, sp5 = 1)
  tax <- data.frame(
    Species = paste0("sp", 1:5),
    Genus   = c("G1", "G1", "G2", "G2", "G3"),
    Family  = rep("F1", 5),
    stringsAsFactors = FALSE
  )
  r <- ozkan_pto(comm, tax)

  # In the nk=0 slice there are 5 species -> Ed_S = ln(5) = 1.6094
  expect_equal(r$Ed_levels[["Species"]], log(5), tolerance = 1e-10)
})

test_that("ozkan_pto genus-level Deng entropy with equal proportions", {
  # 6 species, 3 genera (2 species per genus): genus-level Deng entropy
  # m(Fi) = 2/6 = 1/3 (proportion of each genus)
  # |Fi| = 2 (2 species per genus)
  #
  # Ed_genus = -3 x (1/3) x ln((1/3) / (2^2-1))
  #          = -3 x (1/3) x ln(1/9)
  #          = -ln(1/9) = ln(9) = 2.197
  comm <- setNames(rep(1, 6), paste0("sp", 1:6))
  tax <- data.frame(
    Species = paste0("sp", 1:6),
    Genus = rep(c("G1", "G2", "G3"), each = 2),
    Family = rep("F1", 6),
    stringsAsFactors = FALSE
  )
  r <- ozkan_pto(comm, tax)

  # Genus-level Deng entropy should equal ln(9)
  expect_equal(r$Ed_levels[["Genus"]], log(9), tolerance = 1e-10)
})

test_that("ozkan_pto matches Ozkan (2018) hypothetical examples", {
  # Hypothetical example communities from Section 4 of Ozkan (2018)
  # These values have been verified against the tables in the paper

  # Community A: 12 species, 3 genera (4 species per genus), 1 family
  # Expected uTO+ = ln(54.6) ~ 4.00003
  comm_A <- setNames(rep(1, 12), paste0("sp", 1:12))
  tax_A <- data.frame(
    Species = paste0("sp", 1:12),
    Genus = rep(c("G1", "G2", "G3"), each = 4),
    Family = rep("F1", 12),
    stringsAsFactors = FALSE
  )
  r_A <- ozkan_pto(comm_A, tax_A)
  expect_equal(r_A$uTO_plus, log(54.6), tolerance = 1e-4)

  # Community B: 6 species, 3 genera (2 species per genus), 1 family
  # Expected uTO+ = ln(35) ~ 3.5554
  comm_B <- setNames(rep(1, 6), paste0("sp", 1:6))
  tax_B <- data.frame(
    Species = paste0("sp", 1:6),
    Genus = rep(c("G1", "G2", "G3"), each = 2),
    Family = rep("F1", 6),
    stringsAsFactors = FALSE
  )
  r_B <- ozkan_pto(comm_B, tax_B)
  expect_equal(r_B$uTO_plus, log(35), tolerance = 1e-4)

  # Community A should be more diverse than B (12 species > 6 species)
  expect_true(r_A$uTO_plus > r_B$uTO_plus)
})


# --- max_level parameter tests ---

test_that("ozkan_pto max_level='auto' detects informative levels", {
  # 5 species, 2 genera, 1 family -- Ed=0 at the Family level (single group)
  # Therefore max_informative_level = 1 (only Genus is informative)
  comm <- c(sp1 = 4, sp2 = 2, sp3 = 3, sp4 = 1, sp5 = 2)
  tax <- data.frame(
    Species = paste0("sp", 1:5),
    Genus   = c("G1", "G1", "G1", "G2", "G2"),
    Family  = c("F1", "F1", "F1", "F1", "F1"),
    stringsAsFactors = FALSE
  )

  r <- ozkan_pto(comm, tax, max_level = "auto")

  # There is only one group at the Family level (F1) -> Ed = 0
  # max_informative_level should be 1 (Genus)
  expect_equal(r$max_informative_level, 1L)

  # max versions should be <= full versions
  # (fewer levels in the product = smaller or equal value)
  expect_true(r$uTO_plus_max <= r$uTO_plus + 1e-10)
  expect_true(r$TO_plus_max <= r$TO_plus + 1e-10)
})


test_that("ozkan_pto max_level integer restricts levels", {
  # 5 species, 3 levels (Genus, Family, Order)
  comm <- c(sp1 = 4, sp2 = 2, sp3 = 3, sp4 = 1, sp5 = 2)
  tax <- data.frame(
    Species = paste0("sp", 1:5),
    Genus   = c("G1", "G1", "G2", "G2", "G3"),
    Family  = c("F1", "F1", "F1", "F2", "F2"),
    Order   = c("O1", "O1", "O1", "O1", "O1"),
    stringsAsFactors = FALSE
  )

  r_full <- ozkan_pto(comm, tax)
  r_lim1 <- ozkan_pto(comm, tax, max_level = 1)
  r_lim2 <- ozkan_pto(comm, tax, max_level = 2)

  # max_level=1 -> only Genus level
  # max versions should be smaller or equal
  expect_true(r_lim1$uTO_plus_max <= r_full$uTO_plus + 1e-10)

  # max_level=2 >= max_level=1 (more levels = larger value)
  expect_true(r_lim2$uTO_plus_max >= r_lim1$uTO_plus_max - 1e-10)

  # Full versions should be independent of max_level -- they should remain the same
  expect_equal(r_lim1$uTO, r_full$uTO, tolerance = 1e-10)
  expect_equal(r_lim2$uTO_plus, r_full$uTO_plus, tolerance = 1e-10)
})


test_that("ozkan_pto max_level validates input", {
  comm <- c(sp1 = 4, sp2 = 2, sp3 = 3)
  tax <- data.frame(
    Species = paste0("sp", 1:3),
    Genus   = c("G1", "G1", "G2"),
    Family  = c("F1", "F1", "F1"),
    stringsAsFactors = FALSE
  )

  # max_level = 0 -> error (must be at least 1)
  expect_error(ozkan_pto(comm, tax, max_level = 0), "must be between")

  # max_level > number of available levels -> error
  expect_error(ozkan_pto(comm, tax, max_level = 10), "must be between")

  # Invalid type
  expect_error(ozkan_pto(comm, tax, max_level = TRUE), "must be NULL")
})


test_that("ozkan_pto max versions equal full when all levels informative", {
  # If all levels are informative (Ed > 0), max should equal full
  comm <- c(sp1 = 4, sp2 = 2, sp3 = 3, sp4 = 1, sp5 = 2, sp6 = 1)
  tax <- data.frame(
    Species = paste0("sp", 1:6),
    Genus   = c("G1", "G1", "G2", "G2", "G3", "G3"),
    Family  = c("F1", "F1", "F1", "F2", "F2", "F2"),
    stringsAsFactors = FALSE
  )

  r <- ozkan_pto(comm, tax, max_level = "auto")

  # 2 families (F1, F2) -> Ed > 0 at the Family level as well
  # max should equal full
  expect_equal(r$uTO_max, r$uTO, tolerance = 1e-10)
  expect_equal(r$TO_max, r$TO, tolerance = 1e-10)
  expect_equal(r$uTO_plus_max, r$uTO_plus, tolerance = 1e-10)
  expect_equal(r$TO_plus_max, r$TO_plus, tolerance = 1e-10)
})


test_that("pto_components returns max values matching ozkan_pto auto", {
  comm <- c(sp1 = 4, sp2 = 2, sp3 = 3, sp4 = 1, sp5 = 2)
  tax <- data.frame(
    Species = paste0("sp", 1:5),
    Genus   = c("G1", "G1", "G1", "G2", "G2"),
    Family  = c("F1", "F1", "F1", "F1", "F1"),
    stringsAsFactors = FALSE
  )

  pto <- pto_components(comm, tax)
  full <- ozkan_pto(comm, tax, max_level = "auto")

  expect_equal(unname(pto["uTO_max"]), full$uTO_max, tolerance = 1e-10)
  expect_equal(unname(pto["TO_max"]), full$TO_max, tolerance = 1e-10)
  expect_equal(unname(pto["uTO_plus_max"]), full$uTO_plus_max, tolerance = 1e-10)
  expect_equal(unname(pto["TO_plus_max"]), full$TO_plus_max, tolerance = 1e-10)
})
