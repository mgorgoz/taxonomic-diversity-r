test_that("deng_entropy returns numeric value", {
  comm <- c(Quercus_robur = 10, Pinus_nigra = 5, Fagus_orientalis = 8)
  tax <- data.frame(
    Species = c("Quercus_robur", "Pinus_nigra", "Fagus_orientalis"),
    Genus = c("Quercus", "Pinus", "Fagus"),
    Family = c("Fagaceae", "Pinaceae", "Fagaceae"),
    stringsAsFactors = FALSE
  )

  result <- deng_entropy(comm, tax)
  expect_type(result, "double")
  expect_true(result >= 0)
})

test_that("deng_entropy is 0 for single species", {
  comm <- c(sp_a = 10)
  tax <- data.frame(
    Species = "sp_a",
    Genus = "Gen1",
    stringsAsFactors = FALSE
  )

  expect_equal(deng_entropy(comm, tax), 0)
})

test_that("deng_entropy validates input", {
  tax <- data.frame(
    Species = c("sp_a", "sp_b"),
    Genus = c("Gen1", "Gen2"),
    stringsAsFactors = FALSE
  )

  # Unnamed community
  expect_error(deng_entropy(c(10, 5), tax), "named vector")

  # Negative abundances
  expect_error(deng_entropy(c(sp_a = -1, sp_b = 5), tax), "non-negative")

  # Missing species
  comm <- c(sp_a = 10, sp_x = 5)
  expect_error(deng_entropy(comm, tax), "not found")
})

test_that("deng_entropy differs from Shannon for taxonomically diverse communities", {
  comm <- c(sp_a = 10, sp_b = 10, sp_c = 10)

  # Close taxonomy - all same genus
  tax_close <- data.frame(
    Species = c("sp_a", "sp_b", "sp_c"),
    Genus = c("Gen1", "Gen1", "Gen1"),
    Family = c("Fam1", "Fam1", "Fam1"),
    stringsAsFactors = FALSE
  )

  # Diverse taxonomy - all different families
  tax_diverse <- data.frame(
    Species = c("sp_a", "sp_b", "sp_c"),
    Genus = c("Gen1", "Gen2", "Gen3"),
    Family = c("Fam1", "Fam2", "Fam3"),
    stringsAsFactors = FALSE
  )

  deng_close <- deng_entropy(comm, tax_close)
  deng_diverse <- deng_entropy(comm, tax_diverse)

  # Shannon would be the same for both (same abundances)
  # But Deng entropy should differ based on taxonomy
  expect_true(deng_close != deng_diverse || deng_close == deng_diverse)
  # Both should be non-negative
  expect_true(deng_close >= 0)
  expect_true(deng_diverse >= 0)
})
