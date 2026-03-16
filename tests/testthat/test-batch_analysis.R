# =============================================================================
# test-batch_analysis.R
# Tests for the batch_analysis() function
# =============================================================================

# --- Common test data ---
# Single site data (no Site column)
df_single <- data.frame(
  Species   = c("sp1", "sp2", "sp3", "sp4"),
  Genus     = c("G1", "G1", "G2", "G2"),
  Family    = c("F1", "F1", "F1", "F2"),
  Order     = c("O1", "O1", "O1", "O1"),
  Abundance = c(10, 20, 15, 5),
  stringsAsFactors = FALSE
)

# Multi-site data (with Site column)
df_multi <- data.frame(
  Site      = c("A", "A", "A", "A", "B", "B", "B", "B"),
  Species   = c("sp1", "sp2", "sp3", "sp4", "sp1", "sp2", "sp3", "sp4"),
  Genus     = c("G1", "G1", "G2", "G2", "G1", "G1", "G2", "G2"),
  Family    = c("F1", "F1", "F1", "F2", "F1", "F1", "F1", "F2"),
  Order     = c("O1", "O1", "O1", "O1", "O1", "O1", "O1", "O1"),
  Abundance = c(10, 20, 15, 5, 5, 5, 5, 5),
  stringsAsFactors = FALSE
)

# Taxonomic tree (compatible with compare_indices tests)
tax <- data.frame(
  Species = c("sp1", "sp2", "sp3", "sp4"),
  Genus   = c("G1", "G1", "G2", "G2"),
  Family  = c("F1", "F1", "F1", "F2"),
  Order   = c("O1", "O1", "O1", "O1"),
  stringsAsFactors = FALSE
)


# ---- Test 1: Single site — without Site column ----
test_that("single site data returns correct results", {
  result <- batch_analysis(df_single)

  # Is it a data frame?
  expect_true(is.data.frame(result))

  # Should have 1 row (single site)
  expect_equal(nrow(result), 1)

  # 16 columns: Site + N_Species + 14 indices (6 classic + 4 PTO + 4 PTO max)
  expect_equal(ncol(result), 16)

  # Site column should be "All" (automatic)
  expect_equal(result$Site, "All")

  # Should have N_Species
  expect_equal(result$N_Species, 4L)

  # All index columns should be numeric
  index_cols <- setdiff(names(result), "Site")
  for (col in index_cols) {
    expect_true(is.numeric(result[[col]]),
                info = paste(col, "should be numeric"))
  }
})


# ---- Test 2: Multiple sites — with Site column ----
test_that("multi-site data returns correct results", {
  result <- batch_analysis(df_multi)

  # Should have 2 rows
  expect_equal(nrow(result), 2)

  # Site names should be preserved
  expect_equal(result$Site, c("A", "B"))

  # Each row should have 16 columns
  expect_equal(ncol(result), 16)

  # N_Species for each site
  expect_equal(result$N_Species, c(4L, 4L))
})


# ---- Test 3: Consistency with compare_indices ----
# batch_analysis results should match compare_indices results for the same data
test_that("single site results match compare_indices", {
  batch_result <- batch_analysis(df_single)

  # Run compare_indices with the same data
  comm <- c(sp1 = 10, sp2 = 20, sp3 = 15, sp4 = 5)
  ci_result <- compare_indices(comm, tax)

  # Does Shannon match?
  expect_equal(batch_result$Shannon, ci_result$Shannon)

  # Does Simpson match?
  expect_equal(batch_result$Simpson, ci_result$Simpson)

  # Does Delta match?
  expect_equal(batch_result$Delta, ci_result$Delta)

  # Does Delta_star match?
  expect_equal(batch_result$Delta_star, ci_result$Delta_star)

  # Does AvTD match?
  expect_equal(batch_result$AvTD, ci_result$AvTD)

  # Does VarTD match?
  expect_equal(batch_result$VarTD, ci_result$VarTD)

  # pTO components
  expect_equal(batch_result$uTO, ci_result$uTO)
  expect_equal(batch_result$TO, ci_result$TO)
  expect_equal(batch_result$uTO_plus, ci_result$uTO_plus)
  expect_equal(batch_result$TO_plus, ci_result$TO_plus)

  # pTO max components
  expect_equal(batch_result$uTO_max, ci_result$uTO_max)
  expect_equal(batch_result$TO_max, ci_result$TO_max)
  expect_equal(batch_result$uTO_plus_max, ci_result$uTO_plus_max)
  expect_equal(batch_result$TO_plus_max, ci_result$TO_plus_max)
})


# ---- Test 4: Multiple sites — match compare_indices for each site ----
test_that("multi-site results match compare_indices", {
  batch_result <- batch_analysis(df_multi)

  # Site A: comm_A = c(sp1=10, sp2=20, sp3=15, sp4=5)
  comm_A <- c(sp1 = 10, sp2 = 20, sp3 = 15, sp4 = 5)
  ci_A <- compare_indices(comm_A, tax)
  expect_equal(batch_result$Shannon[1], ci_A$Shannon)
  expect_equal(batch_result$Simpson[1], ci_A$Simpson)

  # Site B: comm_B = c(sp1=5, sp2=5, sp3=5, sp4=5)
  comm_B <- c(sp1 = 5, sp2 = 5, sp3 = 5, sp4 = 5)
  ci_B <- compare_indices(comm_B, tax)
  expect_equal(batch_result$Shannon[2], ci_B$Shannon)
  expect_equal(batch_result$Simpson[2], ci_B$Simpson)
})


# ---- Test 5: Automatic Site column detection ----
# "Site", "Alan", "Plot" names should be automatically detected
test_that("Site column is automatically detected", {
  # With the name "Site"
  result_site <- batch_analysis(df_multi)
  expect_equal(nrow(result_site), 2)

  # With the name "Alan"
  df_alan <- df_multi
  names(df_alan)[1] <- "Alan"
  result_alan <- batch_analysis(df_alan)
  expect_equal(nrow(result_alan), 2)

  # With the name "Plot"
  df_plot <- df_multi
  names(df_plot)[1] <- "Plot"
  result_plot <- batch_analysis(df_plot)
  expect_equal(nrow(result_plot), 2)
})


# ---- Test 6: With site_column parameter ----
test_that("specifying site_column parameter works", {
  df_custom <- df_multi
  names(df_custom)[1] <- "Lokasyon"
  result <- batch_analysis(df_custom, site_column = "Lokasyon")
  expect_equal(nrow(result), 2)
  expect_equal(result$Site, c("A", "B"))
})


# ---- Test 7: Empty Site column — should work as single site ----
test_that("empty Site column works as single site", {
  df_empty_site <- df_single
  df_empty_site$Site <- ""
  result <- batch_analysis(df_empty_site)
  expect_equal(nrow(result), 1)
  expect_equal(result$Site, "All")
})


# ---- Test 8: NA Site column — should work as single site ----
test_that("NA Site column works as single site", {
  df_na_site <- df_single
  df_na_site$Site <- NA
  result <- batch_analysis(df_na_site)
  expect_equal(nrow(result), 1)
  expect_equal(result$Site, "All")
})


# ---- Test 9: Case-insensitive abundance column ----
test_that("abundance column matches case-insensitively", {
  df_lower <- df_single
  names(df_lower)[names(df_lower) == "Abundance"] <- "abundance"
  result <- batch_analysis(df_lower)
  expect_equal(nrow(result), 1)
})


# ---- Test 10: Invalid input checks ----
test_that("throws error on invalid input", {
  # Not a data.frame
  expect_error(batch_analysis("not_a_df"), "data.*must be a data frame")

  # Empty data.frame
  expect_error(batch_analysis(data.frame()), "no rows")

  # No Abundance column
  df_no_abd <- df_single[, -5]
  expect_error(batch_analysis(df_no_abd), "Abundance.*not found")

  # Insufficient taxonomic columns
  df_no_tax <- data.frame(
    Species   = c("sp1", "sp2"),
    Abundance = c(10, 20),
    stringsAsFactors = FALSE
  )
  expect_error(batch_analysis(df_no_tax), "auto-detect")

  # Specified site_column not found
  expect_error(batch_analysis(df_single, site_column = "Nonexistent"),
               "not found")
})


# ---- Test 11: With tax_columns parameter ----
test_that("works with tax_columns parameter", {
  df_custom_names <- data.frame(
    Tur     = c("sp1", "sp2", "sp3", "sp4"),
    Cins    = c("G1", "G1", "G2", "G2"),
    Familya = c("F1", "F1", "F1", "F2"),
    Takim   = c("O1", "O1", "O1", "O1"),
    Bolluk  = c(10, 20, 15, 5),
    stringsAsFactors = FALSE
  )
  result <- batch_analysis(df_custom_names,
                           tax_columns = c("Tur", "Cins", "Familya", "Takim"),
                           abundance_column = "Bolluk")
  expect_equal(nrow(result), 1)
  expect_true(is.numeric(result$Shannon))
})


# ---- Test 12: Diversity is higher with even distribution ----
test_that("diversity is higher with even distribution", {
  result <- batch_analysis(df_multi)

  # Site B (even distribution) should have higher Shannon
  expect_gt(result$Shannon[result$Site == "B"],
            result$Shannon[result$Site == "A"])
})


# ---- Test 13: Automatic detection with Turkish site name ----
test_that("automatic detection works with site name", {
  df_alan <- df_multi
  names(df_alan)[1] <- "alan"  # lowercase
  result <- batch_analysis(df_alan)
  expect_equal(nrow(result), 2)
})


# ---- Test 14: 3 or more sites ----
test_that("works with 3 sites", {
  df_three <- rbind(
    data.frame(Site = "X", df_single, stringsAsFactors = FALSE),
    data.frame(Site = "Y", df_single, stringsAsFactors = FALSE),
    data.frame(Site = "Z", df_single, stringsAsFactors = FALSE)
  )
  # Different abundance at site Y
  df_three$Abundance[df_three$Site == "Y"] <- c(5, 5, 5, 5)
  df_three$Abundance[df_three$Site == "Z"] <- c(50, 1, 1, 1)

  result <- batch_analysis(df_three)
  expect_equal(nrow(result), 3)
  expect_equal(result$Site, c("X", "Y", "Z"))
})


# ---- Test 15: Zero abundance species should be filtered ----
test_that("zero abundance species are filtered", {
  df_zeros <- data.frame(
    Species   = c("sp1", "sp2", "sp3", "sp4", "sp5"),
    Genus     = c("G1", "G1", "G2", "G2", "G3"),
    Family    = c("F1", "F1", "F1", "F2", "F2"),
    Order     = c("O1", "O1", "O1", "O1", "O1"),
    Abundance = c(10, 20, 15, 5, 0),
    stringsAsFactors = FALSE
  )
  # Should not throw an error — sp5 (abundance = 0) should be skipped
  result <- batch_analysis(df_zeros)
  expect_equal(nrow(result), 1)
  expect_true(is.numeric(result$Shannon))
})
