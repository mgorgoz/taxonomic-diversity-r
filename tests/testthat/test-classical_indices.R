# =============================================================================
# Tests for Classical Diversity Indices and Clarke & Warwick Indices
#
# This file tests the following functions:
#   shannon()     — Shannon diversity index (H')
#   simpson()     — Simpson diversity index (D, 1-D, 1/D)
#   delta()       — Taxonomic diversity (Δ) — abundance-weighted
#   delta_star()  — Taxonomic distinctness (Δ*) — excluding same-species pairs
#   avtd()        — Average taxonomic distinctness (Δ+/AvTD) — presence/absence
#   vartd()       — Variance in taxonomic distinctness (Λ+/VarTD)
# =============================================================================


# =============================================================================
# Shannon Diversity Index (H') Tests
# Formula: H' = -Σ pᵢ × ln(pᵢ)
# =============================================================================

test_that("shannon calculates correctly", {
  # 4 species with equal abundance: H' = ln(4) = 1.386
  # Each species' proportion p = 10/40 = 0.25
  # H' = -4 × 0.25 × ln(0.25) = ln(4)
  comm <- c(10, 10, 10, 10)
  expect_equal(shannon(comm), log(4), tolerance = 1e-10)

  # Single species: no diversity, H' = 0
  expect_equal(shannon(c(10)), 0)

  # Base-2 logarithm (in bits): 2 equal species → H' = 1 bit
  comm <- c(10, 10)
  expect_equal(shannon(comm, base = 2), 1, tolerance = 1e-10)
})

test_that("shannon handles edge cases", {
  # Species with zero abundance should be automatically excluded
  # [10, 0, 5, 0, 8] and [10, 5, 8] should give the same result
  comm <- c(10, 0, 5, 0, 8)
  expect_equal(shannon(comm), shannon(c(10, 5, 8)))

  # All zeros: no species, diversity = 0
  expect_equal(shannon(c(0, 0, 0)), 0)
})

test_that("shannon validates input", {
  # Negative abundance: abundance cannot be negative in nature
  expect_error(shannon(c(-1, 5)), "non-negative")

  # Non-numeric value: abundance must be numeric
  expect_error(shannon("abc"), "numeric")
})


# =============================================================================
# Simpson Diversity Index Tests
# Dominance: D = Σ pᵢ²
# Gini-Simpson: 1 - D
# Inverse Simpson: 1/D
# =============================================================================

test_that("simpson calculates correctly", {
  # 4 species with equal abundance:
  # D = 4 × (0.25)² = 4 × 0.0625 = 0.25
  # Gini-Simpson = 1 - 0.25 = 0.75
  # Inverse Simpson = 1/0.25 = 4
  comm <- c(10, 10, 10, 10)
  expect_equal(simpson(comm, "dominance"), 0.25, tolerance = 1e-10)
  expect_equal(simpson(comm, "gini_simpson"), 0.75, tolerance = 1e-10)
  expect_equal(simpson(comm, "inverse"), 4, tolerance = 1e-10)

  # Single species: D = 1 (complete dominance — only 1 species present)
  expect_equal(simpson(c(10), "dominance"), 1)
})

test_that("simpson validates input", {
  # Negative abundance: should throw an error
  expect_error(simpson(c(-1, 5)), "non-negative")
})


# =============================================================================
# Shannon Bias Correction Tests
# =============================================================================

test_that("shannon miller_madow adds positive correction", {
  comm <- c(10, 5, 8, 3, 12)
  H_naive <- shannon(comm)
  H_mm <- shannon(comm, correction = "miller_madow")
  expect_true(H_mm > H_naive)
  N <- sum(comm)
  S_obs <- length(comm)
  expected <- H_naive + (S_obs - 1) / (2 * N)
  expect_equal(H_mm, expected, tolerance = 1e-10)
})

test_that("shannon miller_madow respects base parameter", {
  comm <- c(10, 5, 8, 3, 12)
  H_mm_nat <- shannon(comm, correction = "miller_madow")
  H_mm_2 <- shannon(comm, base = 2, correction = "miller_madow")
  expect_equal(H_mm_2, H_mm_nat / log(2), tolerance = 1e-10)
})

test_that("shannon grassberger uses digamma", {
  comm <- c(10, 5, 8, 3, 12)
  H_g <- shannon(comm, correction = "grassberger")
  N <- sum(comm)
  expected <- log(N) - (1 / N) * sum(comm * digamma(comm))
  expect_equal(H_g, expected, tolerance = 1e-10)
})

test_that("shannon grassberger respects base parameter", {
  comm <- c(10, 5, 8, 3, 12)
  H_g_nat <- shannon(comm, correction = "grassberger")
  H_g_2 <- shannon(comm, base = 2, correction = "grassberger")
  expect_equal(H_g_2, H_g_nat / log(2), tolerance = 1e-10)
})

test_that("shannon chao_shen applies coverage correction", {
  comm <- c(10, 5, 8, 3, 12)
  H_cs <- shannon(comm, correction = "chao_shen")
  expect_type(H_cs, "double")
  expect_true(is.finite(H_cs))
  expect_true(H_cs > 0)
})

test_that("shannon chao_shen handles many singletons", {
  comm <- c(1, 1, 1, 1, 5, 10)
  H_cs <- shannon(comm, correction = "chao_shen")
  expect_true(is.finite(H_cs))
  expect_true(H_cs > 0)
})

test_that("shannon chao_shen handles all singletons", {
  comm <- c(1, 1, 1, 1, 1)
  H_cs <- shannon(comm, correction = "chao_shen")
  expect_true(is.finite(H_cs))
  expect_true(H_cs > 0)
})

test_that("all corrections return 0 for single species", {
  for (corr in c("none", "miller_madow", "grassberger", "chao_shen")) {
    expect_equal(shannon(c(10), correction = corr), 0)
  }
})

test_that("all corrections converge for large equal samples", {
  comm <- rep(1000, 5)
  H_naive <- shannon(comm)
  H_mm <- shannon(comm, correction = "miller_madow")
  H_g <- shannon(comm, correction = "grassberger")
  H_cs <- shannon(comm, correction = "chao_shen")
  expect_equal(H_mm, H_naive, tolerance = 0.01)
  expect_equal(H_g, H_naive, tolerance = 0.01)
  expect_equal(H_cs, H_naive, tolerance = 0.01)
})

test_that("shannon warns for non-integer input with correction", {
  comm <- c(0.5, 0.3, 0.2)
  expect_warning(shannon(comm, correction = "miller_madow"),
                 "non-integer|Non-integer")
})

test_that("shannon default correction is none (backward compatible)", {
  comm <- c(10, 5, 8, 3, 12)
  expect_equal(shannon(comm), shannon(comm, correction = "none"))
})


# =============================================================================
# Clarke & Warwick Indices
# =============================================================================

# --- Δ (Delta): Taxonomic Diversity ---
# Average taxonomic distance between all pairs of individuals
# Pairs of individuals from the same species are also included (distance = 0)

test_that("delta calculates taxonomic diversity correctly", {
  # 3 species, 2 genera, 1 family, all with equal abundance
  comm <- c(sp1 = 3, sp2 = 3, sp3 = 3)
  tax <- data.frame(
    Species = c("sp1", "sp2", "sp3"),
    Genus = c("G1", "G1", "G2"),       # sp1-sp2 same genus, sp3 different
    Family = c("F1", "F1", "F1"),       # all in the same family
    stringsAsFactors = FALSE
  )

  d <- delta(comm, tax)

  # Should return a numeric value and be positive
  expect_type(d, "double")
  expect_true(d > 0)

  # Compare with only same-genus species: distance should be lower
  # Because the taxonomic spread is smaller
  comm2 <- c(sp1 = 3, sp2 = 3)
  tax2 <- data.frame(
    Species = c("sp1", "sp2"),
    Genus = c("G1", "G1"),              # both in the same genus
    Family = c("F1", "F1"),
    stringsAsFactors = FALSE
  )
  d2 <- delta(comm2, tax2)

  # 2 species in the same genus should yield a lower Δ than 3 species spanning different genera
  expect_true(d2 < d)
})

test_that("delta returns 0 for single species", {
  # Single species: no other species to compare with → Δ = 0
  comm <- c(sp1 = 10)
  tax <- data.frame(
    Species = "sp1", Genus = "G1", Family = "F1",
    stringsAsFactors = FALSE
  )
  expect_equal(delta(comm, tax), 0)
})


# --- Δ* (Delta star): Taxonomic Distinctness ---
# Same logic as Δ but pairs of individuals from the same species are EXCLUDED
# As a result, Δ* is always greater than or equal to Δ

test_that("delta_star calculates taxonomic distinctness correctly", {
  # 3 species, different abundances
  comm <- c(sp1 = 5, sp2 = 3, sp3 = 2)
  tax <- data.frame(
    Species = c("sp1", "sp2", "sp3"),
    Genus = c("G1", "G1", "G2"),
    Family = c("F1", "F1", "F1"),
    stringsAsFactors = FALSE
  )

  ds <- delta_star(comm, tax)
  d <- delta(comm, tax)

  # Should be numeric and positive
  expect_type(ds, "double")
  expect_true(ds > 0)

  # MATHEMATICAL RULE: Δ* is always greater than or equal to Δ
  # Because the denominator of Δ also includes same-species pairs (contributing distance = 0)
  # which pulls the average down
  expect_true(ds >= d)
})

test_that("delta_star returns 0 for single species", {
  # Single species: no comparison possible → Δ* = 0
  comm <- c(sp1 = 10)
  tax <- data.frame(
    Species = "sp1", Genus = "G1",
    stringsAsFactors = FALSE
  )
  expect_equal(delta_star(comm, tax), 0)
})


# --- Δ+ (AvTD): Average Taxonomic Distinctness ---
# Uses only presence/absence data (abundance information is not needed)
# Formula: Δ+ = Σ Σ_{i<j} ω_ij / [S(S-1)/2]

test_that("avtd calculates average taxonomic distinctness", {
  # 4 species, 2 genera, 2 families
  tax <- data.frame(
    Species = c("sp1", "sp2", "sp3", "sp4"),
    Genus = c("G1", "G1", "G2", "G2"),
    Family = c("F1", "F1", "F1", "F2"),   # sp4 in a different family
    stringsAsFactors = FALSE
  )

  spp <- c("sp1", "sp2", "sp3", "sp4")
  result <- avtd(spp, tax)

  # Should be numeric and positive
  expect_type(result, "double")
  expect_true(result > 0)
})

test_that("avtd with all same-genus species equals weights[1]", {
  # All species in the same genus → all pairwise distances = weights[1] = 1
  # Therefore the average is also = 1
  tax <- data.frame(
    Species = c("sp1", "sp2", "sp3"),
    Genus = c("G1", "G1", "G1"),     # all in the same genus
    Family = c("F1", "F1", "F1"),
    stringsAsFactors = FALSE
  )

  result <- avtd(c("sp1", "sp2", "sp3"), tax)

  # All distances are 1 → average = 1
  expect_equal(result, 1)
})

test_that("avtd requires at least 2 species", {
  # At least 2 species required (to form pairs)
  tax <- data.frame(
    Species = "sp1", Genus = "G1",
    stringsAsFactors = FALSE
  )
  expect_error(avtd("sp1", tax), "At least 2")
})


# --- Λ+ (VarTD): Variance in Taxonomic Distinctness ---
# Variance of pairwise distances around Δ+
# High Λ+ = some pairs very close, others very distant (uneven distribution)
# Low Λ+ = all pairwise distances are similar (homogeneous distribution)

test_that("vartd calculates variation in taxonomic distinctness", {
  # 4 species: some close, others distant → variance > 0
  tax <- data.frame(
    Species = c("sp1", "sp2", "sp3", "sp4"),
    Genus = c("G1", "G1", "G2", "G2"),
    Family = c("F1", "F1", "F1", "F2"),
    stringsAsFactors = FALSE
  )

  spp <- c("sp1", "sp2", "sp3", "sp4")
  result <- vartd(spp, tax)

  # Should be numeric and greater than or equal to zero (variance cannot be negative)
  expect_type(result, "double")
  expect_true(result >= 0)
})

test_that("vartd is zero when all pairwise distances are equal", {
  # When all pairwise distances are equal, variance = 0
  # 3 species, all in different genera but the same family → all distances = 2
  tax <- data.frame(
    Species = c("sp1", "sp2", "sp3"),
    Genus = c("G1", "G2", "G3"),       # all in different genera
    Family = c("F1", "F1", "F1"),       # all in the same family
    stringsAsFactors = FALSE
  )

  result <- vartd(c("sp1", "sp2", "sp3"), tax)

  # All distances are equal (= 2) → variance = 0
  expect_equal(result, 0, tolerance = 1e-10)
})

test_that("vartd requires at least 2 species", {
  # At least 2 species required
  tax <- data.frame(
    Species = "sp1", Genus = "G1",
    stringsAsFactors = FALSE
  )
  expect_error(vartd("sp1", tax), "At least 2")
})


# =============================================================================
# Error Handling Tests — Invalid Input Cases
# =============================================================================

test_that("delta errors on species not in tax_tree", {
  # Should throw an error when a species name not found in the taxonomy is provided
  comm <- c(sp1 = 5, sp_unknown = 3)
  tax <- data.frame(
    Species = c("sp1", "sp2"),
    Genus = c("G1", "G2"),
    Family = c("F1", "F1"),
    stringsAsFactors = FALSE
  )
  expect_error(delta(comm, tax), "Species not found in tax_tree: sp_unknown")
})

test_that("delta_star errors on species not in tax_tree", {
  # If there are multiple missing species, list them all
  comm <- c(sp1 = 5, sp_missing = 3, sp_also = 2)
  tax <- data.frame(
    Species = c("sp1", "sp2"),
    Genus = c("G1", "G2"),
    Family = c("F1", "F1"),
    stringsAsFactors = FALSE
  )
  expect_error(delta_star(comm, tax), "Species not found in tax_tree: sp_missing, sp_also")
})

test_that("delta errors on wrong weights length", {
  # The length of the weight vector must equal the number of taxonomic levels
  # Here there are 2 levels (Genus, Family) but 3 weights were provided
  comm <- c(sp1 = 5, sp2 = 3)
  tax <- data.frame(
    Species = c("sp1", "sp2"),
    Genus = c("G1", "G2"),
    Family = c("F1", "F1"),
    stringsAsFactors = FALSE
  )
  expect_error(delta(comm, tax, weights = c(1, 2, 3)), "must have length 2")
})

test_that("delta_star errors on wrong weights length", {
  # The same weight length check applies to delta_star as well
  # There are 2 levels but only 1 weight was provided
  comm <- c(sp1 = 5, sp2 = 3)
  tax <- data.frame(
    Species = c("sp1", "sp2"),
    Genus = c("G1", "G2"),
    Family = c("F1", "F1"),
    stringsAsFactors = FALSE
  )
  expect_error(delta_star(comm, tax, weights = c(1)), "must have length 2")
})
