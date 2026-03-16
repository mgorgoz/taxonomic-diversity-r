# =============================================================================
# Deng Entropy Function Tests
# Tests that the deng_entropy_level() function works correctly.
#
# Deng entropy is a generalization of Shannon entropy.
# At the species level (|Fi|=1), it equals Shannon entropy.
# At higher taxonomic levels (|Fi|>1), it is greater than Shannon entropy.
# =============================================================================

test_that("deng_entropy_level returns Shannon entropy at species level", {
  # At species level, Deng should equal Shannon (because |Fi|=1)

  # 4 species with equal abundances: H = ln(4) = 1.386
  # Each species proportion p = 10/40 = 0.25
  # H = -4 × 0.25 × ln(0.25) = ln(4)
  result <- deng_entropy_level(c(10, 10, 10, 10))
  expect_equal(result, log(4), tolerance = 1e-10)

  # Single species: no diversity, H = 0
  expect_equal(deng_entropy_level(c(10)), 0)

  # 2 equal species: H = ln(2) = 0.693
  expect_equal(deng_entropy_level(c(5, 5)), log(2), tolerance = 1e-10)
})

test_that("deng_entropy_level with group_sizes differs from Shannon", {
  # At higher taxonomic levels (genus, family, etc.) Deng ≠ Shannon
  # Because each group contains more than one species (|Fi| > 1)
  # In this case, the 2^|Fi|-1 denominator in the formula comes into play

  # 3 groups: abundances [9, 3, 7], group sizes [3, 2, 3]
  # That is, group 1 has 3 species, group 2 has 2 species, group 3 has 3 species
  abund <- c(9, 3, 7)
  gs <- c(3, 2, 3)

  ed_deng <- deng_entropy_level(abund, group_sizes = gs)
  ed_shannon <- deng_entropy_level(abund)  # no group_sizes = Shannon

  # Deng should be greater than Shannon (when |Fi| > 1)
  # Because the 2^|Fi|-1 term increases the entropy value
  expect_true(ed_deng > ed_shannon)
})

test_that("deng_entropy_level handles zeros", {
  # Species with zero abundances should be automatically excluded from the calculation
  # [10, 0, 5, 0, 8] and [10, 5, 8] should give the same result
  result <- deng_entropy_level(c(10, 0, 5, 0, 8))
  expected <- deng_entropy_level(c(10, 5, 8))
  expect_equal(result, expected)
})

test_that("deng_entropy_level validates input", {
  # Should throw an error when given negative values
  expect_error(deng_entropy_level(c(-1, 5)), "non-negative")

  # Should throw an error when given non-numeric values
  expect_error(deng_entropy_level("abc"), "numeric")
})


# =============================================================================
# Deng Entropy — Bias Correction Tests
# =============================================================================

test_that("deng_entropy_level passes correction to species level", {
  comm <- c(10, 5, 8, 3, 12)
  # At species level, Deng(correction=X) should equal Shannon(correction=X)
  for (corr in c("none", "miller_madow", "grassberger", "chao_shen")) {
    expect_equal(
      deng_entropy_level(comm, correction = corr),
      shannon(comm, correction = corr),
      tolerance = 1e-10
    )
  }
})

test_that("deng_entropy_level warns when correction used with group_sizes", {
  expect_warning(
    deng_entropy_level(c(9, 3, 7), group_sizes = c(3, 2, 3),
                       correction = "chao_shen"),
    "species level"
  )
})

test_that("deng_entropy_level with group_sizes and correction falls back to none", {
  # Even if it gives a warning, the result should be the same as correction = "none"
  result_none <- deng_entropy_level(c(9, 3, 7), group_sizes = c(3, 2, 3),
                                     correction = "none")
  result_cs <- suppressWarnings(
    deng_entropy_level(c(9, 3, 7), group_sizes = c(3, 2, 3),
                       correction = "chao_shen")
  )
  expect_equal(result_cs, result_none)
})
