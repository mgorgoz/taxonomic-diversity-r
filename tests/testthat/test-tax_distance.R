test_that("tax_distance_matrix computes correct distances", {
  tax <- data.frame(
    Species = c("sp_a", "sp_b", "sp_c"),
    Genus = c("Gen1", "Gen1", "Gen2"),
    Family = c("Fam1", "Fam1", "Fam1"),
    stringsAsFactors = FALSE
  )

  d <- tax_distance_matrix(tax)

  # Same genus, same family -> differ at 0 levels
  expect_equal(d["sp_a", "sp_b"], 0)

  # Different genus, same family -> differ at genus level (weight 0.5)
  expect_equal(d["sp_a", "sp_c"], 0.5)

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

  # With equal weights (default)
  d1 <- tax_distance_matrix(tax)
  expect_equal(d1["sp_a", "sp_b"], 1)

  # With custom weights
  d2 <- tax_distance_matrix(tax, weights = c(1, 3))
  expect_equal(d2["sp_a", "sp_b"], 1)  # Normalized, so still sums to 1
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
