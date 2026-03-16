# =============================================================================
# test-compare_indices.R
# Tests for the compare_indices() function
# =============================================================================

# --- Common test data ---
# Simple taxonomic tree with 4 species
tax <- data.frame(
  Species = c("sp1", "sp2", "sp3", "sp4"),
  Genus   = c("G1", "G1", "G2", "G2"),
  Family  = c("F1", "F1", "F1", "F2"),
  Order   = c("O1", "O1", "O1", "O1"),
  stringsAsFactors = FALSE
)

# Two different communities: one evenly distributed, one with a dominant species
comm_equal   <- c(sp1 = 5, sp2 = 5, sp3 = 5, sp4 = 5)
comm_unequal <- c(sp1 = 50, sp2 = 2, sp3 = 1, sp4 = 1)


# ---- Test 1: Single community input should return correct data.frame ----
# When a single vector is provided, it is wrapped with the name "Community"
test_that("single community vector returns correct table", {
  result <- compare_indices(comm_equal, tax)

  # Is it a data frame?
  expect_true(is.data.frame(result))

  # Should have 1 row (single community)
  expect_equal(nrow(result), 1)

  # 16 columns: Community + N_Species + 14 indices (6 classic + 4 PTO + 4 PTO max)
  expect_equal(ncol(result), 16)

  # Check N_Species
  expect_equal(result$N_Species, 4L)

  # Community column should be "Community"
  expect_equal(result$Community, "Community")

  # All index columns should be numeric
  index_cols <- setdiff(names(result), "Community")
  for (col in index_cols) {
    expect_true(is.numeric(result[[col]]),
                info = paste(col, "should be numeric"))
  }
})


# ---- Test 2: Multiple community input ----
# Multiple communities provided as a list
test_that("multiple community list returns correct table", {
  comm_list <- list(
    Esit    = comm_equal,
    Baskin  = comm_unequal
  )
  result <- compare_indices(comm_list, tax)

  # Should have 2 rows
  expect_equal(nrow(result), 2)

  # Community names should be preserved
  expect_equal(result$Community, c("Esit", "Baskin"))
})


# ---- Test 3: Validation with known values ----
# compare_indices results should match individually calculated values
test_that("values match individually calculated indices", {
  result <- compare_indices(comm_equal, tax)

  # Shannon
  expect_equal(result$Shannon, round(shannon(comm_equal), 6))

  # Simpson
  expect_equal(result$Simpson, round(simpson(comm_equal), 6))

  # Delta
  expect_equal(result$Delta, round(delta(comm_equal, tax), 6))

  # Delta*
  expect_equal(result$Delta_star, round(delta_star(comm_equal, tax), 6))

  # AvTD
  sp <- names(comm_equal[comm_equal > 0])
  expect_equal(result$AvTD, round(avtd(sp, tax), 6))

  # VarTD
  expect_equal(result$VarTD, round(vartd(sp, tax), 6))

  # pTO components (remove name differences with unname)
  pto <- pto_components(comm_equal, tax)
  expect_equal(result$uTO, round(unname(pto["uTO"]), 6))
  expect_equal(result$TO, round(unname(pto["TO"]), 6))
  expect_equal(result$uTO_plus, round(unname(pto["uTO_plus"]), 6))
  expect_equal(result$TO_plus, round(unname(pto["TO_plus"]), 6))

  # pTO max components
  expect_equal(result$uTO_max, round(unname(pto["uTO_max"]), 6))
  expect_equal(result$TO_max, round(unname(pto["TO_max"]), 6))
  expect_equal(result$uTO_plus_max, round(unname(pto["uTO_plus_max"]), 6))
  expect_equal(result$TO_plus_max, round(unname(pto["TO_plus_max"]), 6))
})


# ---- Test 4: Even distribution vs dominant species comparison ----
# Shannon should be higher in the evenly distributed community
test_that("Shannon and Simpson are higher with even distribution", {
  comm_list <- list(Esit = comm_equal, Baskin = comm_unequal)
  result <- compare_indices(comm_list, tax)

  esit  <- result[result$Community == "Esit", ]
  baskin <- result[result$Community == "Baskin", ]

  # Shannon is higher with even distribution
  expect_gt(esit$Shannon, baskin$Shannon)

  # Simpson is higher with even distribution
  expect_gt(esit$Simpson, baskin$Simpson)
})


# ---- Test 5: Should return a list when plot = TRUE ----
# If ggplot2 is installed, both table and plot are returned
test_that("returns a list when plot = TRUE", {
  skip_if_not_installed("ggplot2")

  result <- compare_indices(comm_equal, tax, plot = TRUE)

  # Should be a list
  expect_true(is.list(result))
  expect_true("table" %in% names(result))
  expect_true("plot" %in% names(result))

  # table should be a data.frame
  expect_true(is.data.frame(result$table))

  # plot should be a ggplot object
  expect_s3_class(result$plot, "ggplot")
})


# ---- Test 6: Invalid input validation ----
# Should throw an error when an incorrect data type is provided
test_that("throws an error for invalid input", {
  # If tax_tree is not a data.frame
  expect_error(compare_indices(comm_equal, "not_a_df"),
               "tax_tree")

  # If community is not numeric
  expect_error(compare_indices("not_numeric", tax),
               "must be")

  # If community is unnamed
  expect_error(compare_indices(c(1, 2, 3), tax),
               "must be")
})
