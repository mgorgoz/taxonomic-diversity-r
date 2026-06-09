# =============================================================================
# test-batch_analysis.R
# Tests for the batch_analysis() function
# =============================================================================

# --- Common test data ---
# Single site data (no Site column) — Westhoff-Maarel scale (1-9)
df_single <- data.frame(
  Species   = c("sp1", "sp2", "sp3", "sp4"),
  Genus     = c("G1", "G1", "G2", "G2"),
  Family    = c("F1", "F1", "F1", "F2"),
  Order     = c("O1", "O1", "O1", "O1"),
  Abundance = c(7, 3, 5, 2),
  stringsAsFactors = FALSE
)

# Multi-site data (with Site column) — Westhoff-Maarel scale (1-9)
df_multi <- data.frame(
  Site      = c("A", "A", "A", "A", "B", "B", "B", "B"),
  Species   = c("sp1", "sp2", "sp3", "sp4", "sp1", "sp2", "sp3", "sp4"),
  Genus     = c("G1", "G1", "G2", "G2", "G1", "G1", "G2", "G2"),
  Family    = c("F1", "F1", "F1", "F2", "F1", "F1", "F1", "F2"),
  Order     = c("O1", "O1", "O1", "O1", "O1", "O1", "O1", "O1"),
  Abundance = c(7, 3, 5, 2, 4, 4, 4, 4),
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
  comm <- c(sp1 = 7, sp2 = 3, sp3 = 5, sp4 = 2)
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

  # Site A: comm_A = c(sp1=7, sp2=3, sp3=5, sp4=2)
  comm_A <- c(sp1 = 7, sp2 = 3, sp3 = 5, sp4 = 2)
  ci_A <- compare_indices(comm_A, tax)
  expect_equal(batch_result$Shannon[1], ci_A$Shannon)
  expect_equal(batch_result$Simpson[1], ci_A$Simpson)

  # Site B: comm_B = c(sp1=4, sp2=4, sp3=4, sp4=4)
  comm_B <- c(sp1 = 4, sp2 = 4, sp3 = 4, sp4 = 4)
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

  # With the name "Plot"
  df_plot <- df_multi
  names(df_plot)[1] <- "Plot"
  result_plot <- batch_analysis(df_plot)
  expect_equal(nrow(result_plot), 2)

  # Lowercase "site"
  df_lower <- df_multi
  names(df_lower)[1] <- "site"
  result_lower <- batch_analysis(df_lower)
  expect_equal(nrow(result_lower), 2)
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
    Bolluk  = c(7, 3, 5, 2),
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


# ---- Test 13: 3 or more sites ----
test_that("works with 3 sites", {
  df_three <- rbind(
    data.frame(Site = "X", df_single, stringsAsFactors = FALSE),
    data.frame(Site = "Y", df_single, stringsAsFactors = FALSE),
    data.frame(Site = "Z", df_single, stringsAsFactors = FALSE)
  )
  # Different abundance at site Y and Z (Westhoff-Maarel 1-9)
  df_three$Abundance[df_three$Site == "Y"] <- c(4, 4, 4, 4)
  df_three$Abundance[df_three$Site == "Z"] <- c(9, 1, 1, 1)

  result <- batch_analysis(df_three)
  expect_equal(nrow(result), 3)
  expect_equal(result$Site, c("X", "Y", "Z"))
})


# ---- Test 15: indices parameter — classical only ----
test_that("indices = 'classical' returns only Shannon and Simpson", {
  result <- batch_analysis(df_single, indices = "classical")

  # Site + N_Species + Shannon + Simpson = 4 columns

  expect_equal(ncol(result), 4)
  expect_equal(names(result), c("Site", "N_Species", "Shannon", "Simpson"))
  expect_true(is.numeric(result$Shannon))
  expect_true(is.numeric(result$Simpson))
})


# ---- Test 16: indices parameter — clarke_warwick only ----
test_that("indices = 'clarke_warwick' returns only CW indices", {
  result <- batch_analysis(df_single, indices = "clarke_warwick")

  # Site + N_Species + Delta + Delta_star + AvTD + VarTD = 6 columns
  expect_equal(ncol(result), 6)
  expect_equal(names(result),
               c("Site", "N_Species", "Delta", "Delta_star", "AvTD", "VarTD"))
})


# ---- Test 17: indices parameter — ozkan_pto only ----
test_that("indices = 'ozkan_pto' returns only pTO indices", {
  result <- batch_analysis(df_single, indices = "ozkan_pto", full = FALSE)

  # Site + N_Species + 8 pTO = 10 columns
  expect_equal(ncol(result), 10)
  pto_cols <- c("uTO", "TO", "uTO_plus", "TO_plus",
                "uTO_max", "TO_max", "uTO_plus_max", "TO_plus_max")
  expect_true(all(pto_cols %in% names(result)))
})


# ---- Test 18: indices parameter — classical + clarke_warwick ----
test_that("indices = c('classical', 'clarke_warwick') returns 8 columns", {
  result <- batch_analysis(df_single,
                           indices = c("classical", "clarke_warwick"))

  # Site + N_Species + 2 classical + 4 CW = 8 columns
  expect_equal(ncol(result), 8)
  expect_true("Shannon" %in% names(result))
  expect_true("Simpson" %in% names(result))
  expect_true("Delta" %in% names(result))
  expect_true("VarTD" %in% names(result))
  # No pTO columns
  expect_false("uTO" %in% names(result))
})


# ---- Test 19: indices parameter — unambiguous partial matching ----
test_that("unambiguous abbreviations work for indices", {
  # "clas" -> classical
  result_clas <- batch_analysis(df_single, indices = "clas")
  expect_equal(ncol(result_clas), 4)
  expect_true("Shannon" %in% names(result_clas))

  # "clark" -> clarke_warwick
  result_clark <- batch_analysis(df_single, indices = "clark")
  expect_equal(ncol(result_clark), 6)
  expect_true("Delta" %in% names(result_clark))

  # "oz" -> ozkan_pto
  result_oz <- batch_analysis(df_single, indices = "oz", full = FALSE)
  expect_equal(ncol(result_oz), 10)
  expect_true("uTO" %in% names(result_oz))
})


# ---- Test 20: indices parameter — ambiguous partial match gives error ----
test_that("ambiguous abbreviations give error", {
  # "cl" matches both "classical" and "clarke_warwick"
  expect_error(batch_analysis(df_single, indices = "cl"),
               "Unknown index group")

  # "cla" also ambiguous
  expect_error(batch_analysis(df_single, indices = "cla"),
               "Unknown index group")
})


# ---- Test 21: indices parameter — invalid group gives error ----
test_that("invalid index group gives error", {
  expect_error(batch_analysis(df_single, indices = "nonexistent"),
               "Unknown index group")

  expect_error(batch_analysis(df_single, indices = ""),
               "Unknown index group")
})


# ---- Test 22: indices parameter — default returns all 16 columns ----
test_that("default indices returns all 16 columns", {
  result <- batch_analysis(df_single)
  expect_equal(ncol(result), 16)

  # All groups present
  expect_true("Shannon" %in% names(result))
  expect_true("Delta" %in% names(result))
  expect_true("uTO" %in% names(result))
})


# ---- Test 23: indices parameter — multi-site works with selection ----
test_that("indices selection works with multi-site data", {
  result <- batch_analysis(df_multi, indices = "classical")
  expect_equal(nrow(result), 2)
  expect_equal(ncol(result), 4)
  expect_equal(result$Site, c("A", "B"))
})


# ---- Test 24: Zero abundance species should be filtered ----
test_that("zero abundance species are filtered", {
  df_zeros <- data.frame(
    Species   = c("sp1", "sp2", "sp3", "sp4", "sp5"),
    Genus     = c("G1", "G1", "G2", "G2", "G3"),
    Family    = c("F1", "F1", "F1", "F2", "F2"),
    Order     = c("O1", "O1", "O1", "O1", "O1"),
    Abundance = c(7, 3, 5, 2, 0),
    stringsAsFactors = FALSE
  )
  # Should not throw an error — sp5 (abundance = 0) should be skipped
  result <- batch_analysis(df_zeros)
  expect_equal(nrow(result), 1)
  expect_true(is.numeric(result$Shannon))
})


# ---- Test 25: ozkan_pto accepts any positive abundance values ----
test_that("ozkan_pto accepts abundance values outside Westhoff-Maarel 1-9", {
  df_raw_counts <- data.frame(
    Species   = c("sp1", "sp2", "sp3", "sp4"),
    Genus     = c("G1", "G1", "G2", "G2"),
    Family    = c("F1", "F1", "F1", "F2"),
    Order     = c("O1", "O1", "O1", "O1"),
    Abundance = c(10, 20, 15, 5),
    stringsAsFactors = FALSE
  )

  # Raw counts with pTO should work (no hard W-M restriction)
  result <- batch_analysis(df_raw_counts, indices = "ozkan_pto", full = FALSE)
  expect_equal(nrow(result), 1)
  expect_true(is.numeric(result$uTO))
})


# ---- Test 26: W-M validation — no error without ozkan_pto ----
test_that("raw counts work fine without ozkan_pto", {
  df_raw_counts <- data.frame(
    Species   = c("sp1", "sp2", "sp3", "sp4"),
    Genus     = c("G1", "G1", "G2", "G2"),
    Family    = c("F1", "F1", "F1", "F2"),
    Order     = c("O1", "O1", "O1", "O1"),
    Abundance = c(10, 20, 15, 5),
    stringsAsFactors = FALSE
  )

  # Classical only — should work with raw counts
  result_cl <- batch_analysis(df_raw_counts, indices = "classical")
  expect_equal(nrow(result_cl), 1)
  expect_true(is.numeric(result_cl$Shannon))

  # Clarke & Warwick only — should work with raw counts
  result_cw <- batch_analysis(df_raw_counts, indices = "clarke_warwick")
  expect_equal(nrow(result_cw), 1)
  expect_true(is.numeric(result_cw$Delta))

  # Both classical + CW — should work with raw counts
  result_both <- batch_analysis(df_raw_counts,
                                indices = c("classical", "clarke_warwick"))
  expect_equal(nrow(result_both), 1)
  expect_equal(ncol(result_both), 8)
})


# ---- Test 27: W-M validation — boundary values (1 and 9) are valid ----
test_that("W-M boundary values 1 and 9 are accepted", {
  df_boundary <- data.frame(
    Species   = c("sp1", "sp2", "sp3", "sp4"),
    Genus     = c("G1", "G1", "G2", "G2"),
    Family    = c("F1", "F1", "F1", "F2"),
    Order     = c("O1", "O1", "O1", "O1"),
    Abundance = c(1, 9, 5, 3),
    stringsAsFactors = FALSE
  )

  # Should work — all values within 1-9
  result <- batch_analysis(df_boundary, full = FALSE)
  expect_equal(nrow(result), 1)
  expect_true(is.numeric(result$uTO))
})


# ---- Test 28: Natural sort — Site1..Site11 in correct order ----
test_that("sites are naturally sorted (Site1, Site2, ..., Site10)", {
  # Create 11 sites with shuffled order in data
  sites <- paste0("Site", 1:11)
  rows <- lapply(sites, function(s) {
    data.frame(
      Site      = s,
      Species   = c("sp1", "sp2", "sp3"),
      Genus     = c("G1", "G1", "G2"),
      Family    = c("F1", "F1", "F2"),
      Order     = c("O1", "O1", "O1"),
      Abundance = c(5, 3, 2),
      stringsAsFactors = FALSE
    )
  })
  df_11 <- do.call(rbind, rows)

  result <- batch_analysis(df_11, indices = "classical")

  # Should be Site1, Site2, ..., Site10, Site11 (NOT Site1, Site10, Site11, Site2)
  expect_equal(result$Site, paste0("Site", 1:11))
})


# ---- Test 30: Natural sort — mixed prefixes ----
test_that("natural sort works with mixed prefixes", {
  df_mixed <- data.frame(
    Site      = rep(c("A2", "A10", "A1", "B1", "B10", "B2"), each = 3),
    Species   = rep(c("sp1", "sp2", "sp3"), 6),
    Genus     = rep(c("G1", "G1", "G2"), 6),
    Family    = rep(c("F1", "F1", "F2"), 6),
    Order     = rep(c("O1", "O1", "O1"), 6),
    Abundance = rep(c(5, 3, 2), 6),
    stringsAsFactors = FALSE
  )

  result <- batch_analysis(df_mixed, indices = "classical")
  expect_equal(result$Site, c("A1", "A2", "A10", "B1", "B2", "B10"))
})


# ---- Test 31: Natural sort — no numbers (pure alphabetic) ----
test_that("sites without numbers are sorted alphabetically", {
  df_alpha <- data.frame(
    Site      = rep(c("Cherry", "Apple", "Banana"), each = 3),
    Species   = rep(c("sp1", "sp2", "sp3"), 3),
    Genus     = rep(c("G1", "G1", "G2"), 3),
    Family    = rep(c("F1", "F1", "F2"), 3),
    Order     = rep(c("O1", "O1", "O1"), 3),
    Abundance = rep(c(5, 3, 2), 3),
    stringsAsFactors = FALSE
  )

  result <- batch_analysis(df_alpha, indices = "classical")
  expect_equal(result$Site, c("Apple", "Banana", "Cherry"))
})
