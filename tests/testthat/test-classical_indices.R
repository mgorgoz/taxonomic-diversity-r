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

# --- Clarke & Warwick Indices ---

test_that("delta calculates taxonomic diversity correctly", {
  # 3 species, 2 genera, 1 family
  comm <- c(sp1 = 3, sp2 = 3, sp3 = 3)
  tax <- data.frame(
    Species = c("sp1", "sp2", "sp3"),
    Genus = c("G1", "G1", "G2"),
    Family = c("F1", "F1", "F1"),
    stringsAsFactors = FALSE
  )

  d <- delta(comm, tax)
  expect_type(d, "double")
  expect_true(d > 0)

  # With only same-genus species, delta should be lower
  comm2 <- c(sp1 = 3, sp2 = 3)
  tax2 <- data.frame(
    Species = c("sp1", "sp2"),
    Genus = c("G1", "G1"),
    Family = c("F1", "F1"),
    stringsAsFactors = FALSE
  )
  d2 <- delta(comm2, tax2)
  expect_true(d2 < d)
})

test_that("delta returns 0 for single species", {
  comm <- c(sp1 = 10)
  tax <- data.frame(
    Species = "sp1", Genus = "G1", Family = "F1",
    stringsAsFactors = FALSE
  )
  expect_equal(delta(comm, tax), 0)
})

test_that("delta_star calculates taxonomic distinctness correctly", {
  comm <- c(sp1 = 5, sp2 = 3, sp3 = 2)
  tax <- data.frame(
    Species = c("sp1", "sp2", "sp3"),
    Genus = c("G1", "G1", "G2"),
    Family = c("F1", "F1", "F1"),
    stringsAsFactors = FALSE
  )

  ds <- delta_star(comm, tax)
  d <- delta(comm, tax)

  expect_type(ds, "double")
  expect_true(ds > 0)
  # Delta* >= Delta always
  expect_true(ds >= d)
})

test_that("delta_star returns 0 for single species", {
  comm <- c(sp1 = 10)
  tax <- data.frame(
    Species = "sp1", Genus = "G1",
    stringsAsFactors = FALSE
  )
  expect_equal(delta_star(comm, tax), 0)
})

test_that("avtd calculates average taxonomic distinctness", {
  tax <- data.frame(
    Species = c("sp1", "sp2", "sp3", "sp4"),
    Genus = c("G1", "G1", "G2", "G2"),
    Family = c("F1", "F1", "F1", "F2"),
    stringsAsFactors = FALSE
  )

  spp <- c("sp1", "sp2", "sp3", "sp4")
  result <- avtd(spp, tax)

  expect_type(result, "double")
  expect_true(result > 0)
})

test_that("avtd with all same-genus species equals weights[1]", {
  tax <- data.frame(
    Species = c("sp1", "sp2", "sp3"),
    Genus = c("G1", "G1", "G1"),
    Family = c("F1", "F1", "F1"),
    stringsAsFactors = FALSE
  )

  # All species in same genus -> first match at Genus (level 1)
  # -> all pairwise distances = weights[1] = 1 -> avtd = 1
  result <- avtd(c("sp1", "sp2", "sp3"), tax)
  expect_equal(result, 1)
})

test_that("avtd requires at least 2 species", {
  tax <- data.frame(
    Species = "sp1", Genus = "G1",
    stringsAsFactors = FALSE
  )
  expect_error(avtd("sp1", tax), "At least 2")
})

test_that("vartd calculates variation in taxonomic distinctness", {
  tax <- data.frame(
    Species = c("sp1", "sp2", "sp3", "sp4"),
    Genus = c("G1", "G1", "G2", "G2"),
    Family = c("F1", "F1", "F1", "F2"),
    stringsAsFactors = FALSE
  )

  spp <- c("sp1", "sp2", "sp3", "sp4")
  result <- vartd(spp, tax)

  expect_type(result, "double")
  expect_true(result >= 0)
})

test_that("vartd is zero when all pairwise distances are equal", {
  # All species different genus, same family -> all distances equal
  tax <- data.frame(
    Species = c("sp1", "sp2", "sp3"),
    Genus = c("G1", "G2", "G3"),
    Family = c("F1", "F1", "F1"),
    stringsAsFactors = FALSE
  )

  result <- vartd(c("sp1", "sp2", "sp3"), tax)
  expect_equal(result, 0, tolerance = 1e-10)
})

test_that("vartd requires at least 2 species", {
  tax <- data.frame(
    Species = "sp1", Genus = "G1",
    stringsAsFactors = FALSE
  )
  expect_error(vartd("sp1", tax), "At least 2")
})

test_that("delta errors on species not in tax_tree", {
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
  comm <- c(sp1 = 5, sp_missing = 3, sp_also = 2)
  tax <- data.frame(
    Species = c("sp1", "sp2"),
    Genus = c("G1", "G2"),
    Family = c("F1", "F1"),
    stringsAsFactors = FALSE
  )
  expect_error(delta_star(comm, tax), "Species not found in tax_tree: sp_missing, sp_also")
})
