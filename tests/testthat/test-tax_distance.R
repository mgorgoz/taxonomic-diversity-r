# =============================================================================
# Taxonomic Distance Matrix and Taxonomy Tree Builder Tests
#
# This file tests the following functions:
#   tax_distance_matrix() — Pairwise taxonomic distance matrix between species
#   build_tax_tree()      — Classification table (data.frame) builder
#
# Distance calculation logic:
#   Distance between two species = weight of the first matching taxonomic level
#   Same genus → ω=1, Same family but different genus → ω=2, etc.
# =============================================================================

test_that("tax_distance_matrix computes correct distances", {
  # 3 species: sp_a and sp_b are in the same genus, sp_c is in a different genus but the same family
  tax <- data.frame(
    Species = c("sp_a", "sp_b", "sp_c"),
    Genus = c("Gen1", "Gen1", "Gen2"),
    Family = c("Fam1", "Fam1", "Fam1"),
    stringsAsFactors = FALSE
  )

  # Default weights: Genus=1, Family=2
  d <- tax_distance_matrix(tax)

  # sp_a and sp_b are the same genus (Gen1) → first match at Genus level → ω = 1
  expect_equal(d["sp_a", "sp_b"], 1)

  # sp_a and sp_c are different genera (Gen1 ≠ Gen2), same family → match at Family level → ω = 2
  expect_equal(d["sp_a", "sp_c"], 2)

  # Matrix should be symmetric: ω(a,c) = ω(c,a)
  expect_equal(d["sp_a", "sp_c"], d["sp_c", "sp_a"])

  # Diagonal should be zero: a species' distance to itself = 0
  expect_equal(d["sp_a", "sp_a"], 0)
})

test_that("tax_distance_matrix works with custom weights", {
  # 2 species, no match at any level (different genus, different family)
  tax <- data.frame(
    Species = c("sp_a", "sp_b"),
    Genus = c("Gen1", "Gen2"),
    Family = c("Fam1", "Fam2"),
    stringsAsFactors = FALSE
  )

  # Default weights [1, 2]: no match at any level → maximum weight = 2
  d1 <- tax_distance_matrix(tax)
  expect_equal(d1["sp_a", "sp_b"], 2)

  # Custom weights [1, 3]: no match → maximum weight = 3
  d2 <- tax_distance_matrix(tax, weights = c(1, 3))
  expect_equal(d2["sp_a", "sp_b"], 3)
})


# =============================================================================
# build_tax_tree() Tests — Classification Table Builder
# =============================================================================

test_that("build_tax_tree creates correct structure", {
  # 2 species, 2 taxonomic levels (Genus, Family)
  tree <- build_tax_tree(
    species = c("sp_a", "sp_b"),
    Genus = c("Gen1", "Gen2"),
    Family = c("Fam1", "Fam1")
  )

  # Should have 3 columns: Species + Genus + Family
  expect_equal(ncol(tree), 3)
  expect_equal(names(tree), c("Species", "Genus", "Family"))

  # Should have 2 rows (1 for each species)
  expect_equal(nrow(tree), 2)
})

test_that("build_tax_tree validates input", {
  # Number of species must match the length of taxonomic levels
  # 2 species provided but only 1 genus name given → error
  expect_error(
    build_tax_tree(species = c("a", "b"), Genus = c("G1")),
    "does not match"
  )

  # At least 1 taxonomic level must be provided (species names alone are not sufficient)
  expect_error(build_tax_tree(species = c("a")), "At least one")
})

test_that("build_tax_tree rejects duplicate species", {
  # The same species name is given twice → error
  # Because species names must be unique (they serve as row/column names in the distance matrix)
  expect_error(
    build_tax_tree(species = c("sp1", "sp2", "sp1"), Genus = c("G1", "G2", "G1")),
    "Duplicate species names: sp1"
  )
})


# =============================================================================
# Error Handling — Incorrect Weights Length
# =============================================================================

test_that("tax_distance_matrix errors on wrong weights length", {
  # There are 2 taxonomic levels (Genus, Family) → 2 weights required
  tax <- data.frame(
    Species = c("sp_a", "sp_b"),
    Genus = c("Gen1", "Gen2"),
    Family = c("Fam1", "Fam1"),
    stringsAsFactors = FALSE
  )

  # 1 weight provided but 2 required → error
  expect_error(tax_distance_matrix(tax, weights = c(1)), "must have length 2")

  # 4 weights provided but 2 required → error
  expect_error(tax_distance_matrix(tax, weights = c(1, 2, 3, 4)), "must have length 2")
})
