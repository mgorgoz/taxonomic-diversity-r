test_that("tax_distance_matrix computes correct distances", {
  tax <- data.frame(
    Species = c("sp_a", "sp_b", "sp_c"),
    Genus = c("Gen1", "Gen1", "Gen2"),
    Family = c("Fam1", "Fam1", "Fam1"),
    stringsAsFactors = FALSE
  )

  # Default weights = c(1, 2) for 2 taxonomic levels
  d <- tax_distance_matrix(tax)

  # Same genus -> first match at Genus (level 1) -> weights[1] = 1
  expect_equal(d["sp_a", "sp_b"], 1)

  # Different genus, same family -> first match at Family (level 2) -> weights[2] = 2
  expect_equal(d["sp_a", "sp_c"], 2)

  # Symmetric
  expect_equal(d["sp_a", "sp_c"], d["sp_c", "sp_a"])

  # Diagonal is zero
  expect_equal(d["sp_a", "sp_a"], 0)
})

test_that("tax_distance_matrix works with custom weights", {
  tax <- data.frame(
    Species = c("sp_a", "sp_b"),
    Genus = c("Gen1", "Gen2"),
    Family = c("Fam1", "Fam2"),
    stringsAsFactors = FALSE
  )

  # Default weights = c(1, 2); no match at any level -> max weight = 2
  d1 <- tax_distance_matrix(tax)
  expect_equal(d1["sp_a", "sp_b"], 2)

  # Custom weights c(1, 3); no match -> max weight = 3
  d2 <- tax_distance_matrix(tax, weights = c(1, 3))
  expect_equal(d2["sp_a", "sp_b"], 3)
})

test_that("build_tax_tree creates correct structure", {
  tree <- build_tax_tree(
    species = c("sp_a", "sp_b"),
    Genus = c("Gen1", "Gen2"),
    Family = c("Fam1", "Fam1")
  )

  expect_equal(ncol(tree), 3)
  expect_equal(names(tree), c("Species", "Genus", "Family"))
  expect_equal(nrow(tree), 2)
})

test_that("build_tax_tree validates input", {
  expect_error(
    build_tax_tree(species = c("a", "b"), Genus = c("G1")),
    "does not match"
  )
  expect_error(build_tax_tree(species = c("a")), "At least one")
})

test_that("build_tax_tree rejects duplicate species", {
  expect_error(
    build_tax_tree(species = c("sp1", "sp2", "sp1"), Genus = c("G1", "G2", "G1")),
    "Duplicate species names: sp1"
  )
})

test_that("tax_distance_matrix errors on wrong weights length", {
  tax <- data.frame(
    Species = c("sp_a", "sp_b"),
    Genus = c("Gen1", "Gen2"),
    Family = c("Fam1", "Fam1"),
    stringsAsFactors = FALSE
  )
  # 2 levels but 1 weight
  expect_error(tax_distance_matrix(tax, weights = c(1)), "must have length 2")
  # 2 levels but 4 weights
  expect_error(tax_distance_matrix(tax, weights = c(1, 2, 3, 4)), "must have length 2")
})
