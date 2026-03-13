# ---- test fixtures ----
tax_10sp <- data.frame(
  Species = paste0("sp", 1:10),
  Genus   = rep(c("G1", "G2", "G3", "G4", "G5"), each = 2),
  Family  = rep(c("F1", "F1", "F2", "F2", "F3"), each = 2),
  Order   = rep(c("O1", "O1", "O2", "O2", "O2"), each = 2),
  stringsAsFactors = FALSE
)

tax_4sp <- data.frame(
  Species = paste0("sp", 1:4),
  Genus   = c("G1", "G1", "G2", "G2"),
  Family  = c("F1", "F1", "F2", "F2"),
  stringsAsFactors = FALSE
)


# ---- simulate_td() ----

test_that("simulate_td returns correct structure with both indices", {
  sim <- simulate_td(tax_10sp, n_sim = 50, seed = 1)

  expect_s3_class(sim, "td_simulation")
  expect_s3_class(sim, "data.frame")
  expect_true(all(c("s", "mean_avtd", "lower_avtd", "upper_avtd",
                     "mean_vartd", "lower_vartd", "upper_vartd") %in%
                    names(sim)))
  # Default s_range: 2 to 10

expect_equal(nrow(sim), 9)
  expect_equal(sim$s, 2:10)
})

test_that("simulate_td returns correct structure with avtd only", {
  sim <- simulate_td(tax_10sp, n_sim = 50, index = "avtd", seed = 1)

  expect_true("mean_avtd" %in% names(sim))
  expect_false("mean_vartd" %in% names(sim))
  expect_equal(attr(sim, "index"), "avtd")
})

test_that("simulate_td returns correct structure with vartd only", {
  sim <- simulate_td(tax_10sp, n_sim = 50, index = "vartd", seed = 1)

  expect_true("mean_vartd" %in% names(sim))
  expect_false("mean_avtd" %in% names(sim))
  expect_equal(attr(sim, "index"), "vartd")
})

test_that("simulate_td: lower <= mean <= upper", {
  sim <- simulate_td(tax_10sp, n_sim = 200, seed = 42)

  expect_true(all(sim$lower_avtd <= sim$mean_avtd))
  expect_true(all(sim$mean_avtd <= sim$upper_avtd))
  expect_true(all(sim$lower_vartd <= sim$mean_vartd))
  expect_true(all(sim$mean_vartd <= sim$upper_vartd))
})

test_that("simulate_td: seed reproducibility", {
  sim1 <- simulate_td(tax_10sp, n_sim = 50, seed = 123)
  sim2 <- simulate_td(tax_10sp, n_sim = 50, seed = 123)

  expect_identical(sim1, sim2)
})

test_that("simulate_td: custom s_range", {
  sim <- simulate_td(tax_10sp, s_range = c(3, 5, 7), n_sim = 50, seed = 1)

  expect_equal(nrow(sim), 3)
  expect_equal(sim$s, c(3L, 5L, 7L))
})

test_that("simulate_td: attributes are correct", {
  sim <- simulate_td(tax_10sp, n_sim = 100, ci = 0.90, seed = 1)

  expect_equal(attr(sim, "ci"), 0.90)
  expect_equal(attr(sim, "index"), "both")
  expect_equal(attr(sim, "n_sim"), 100L)
  expect_equal(attr(sim, "pool_size"), 10)
})

test_that("simulate_td: small taxonomy works (S=4, min s=2)", {
  sim <- simulate_td(tax_4sp, n_sim = 50, seed = 1)

  expect_equal(nrow(sim), 3)
  expect_equal(sim$s, 2:4)
})

test_that("simulate_td: S == pool_size gives zero variance for avtd", {
  sim <- simulate_td(tax_4sp, s_range = 4, n_sim = 50, seed = 1)

  # When s == pool_size, every draw is the same, so lower == upper
  expect_equal(sim$lower_avtd, sim$upper_avtd)
  expect_equal(sim$mean_avtd, sim$lower_avtd)
})

test_that("simulate_td: errors on invalid input", {
  expect_error(simulate_td("not a df"), "must be a data frame")
  expect_error(simulate_td(data.frame(x = 1)), "at least 2 columns")
  expect_error(
    simulate_td(data.frame(sp = c("a", "b"), g = c("g1", "g2"))),
    "at least 3 species"
  )
  expect_error(simulate_td(tax_10sp, n_sim = 0), "positive integer")
  expect_error(simulate_td(tax_10sp, ci = 1.5), "between 0 and 1")
  expect_error(simulate_td(tax_10sp, s_range = c(0, 1)), "between 2 and")
})

test_that("simulate_td: print method works", {
  sim <- simulate_td(tax_10sp, n_sim = 50, seed = 1)
  expect_output(print(sim), "Taxonomic Distinctness Simulation")
  expect_output(print(sim), "Species pool: 10 species")
})


# ---- plot_funnel() ----

test_that("plot_funnel returns ggplot object", {
  skip_if_not_installed("ggplot2")
  sim <- simulate_td(tax_10sp, n_sim = 50, seed = 1)
  p <- plot_funnel(sim)
  expect_s3_class(p, "ggplot")
})

test_that("plot_funnel: avtd and vartd index selection", {
  skip_if_not_installed("ggplot2")
  sim <- simulate_td(tax_10sp, n_sim = 50, seed = 1)

  p1 <- plot_funnel(sim, index = "avtd")
  p2 <- plot_funnel(sim, index = "vartd")
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
})

test_that("plot_funnel: with observed points", {
  skip_if_not_installed("ggplot2")
  sim <- simulate_td(tax_10sp, n_sim = 50, seed = 1)

  obs <- data.frame(
    site  = c("Site_A", "Site_B"),
    s     = c(5, 8),
    value = c(2.5, 1.8)
  )
  p <- plot_funnel(sim, observed = obs)
  expect_s3_class(p, "ggplot")
  # Should have geom_point layer for observed
  layer_classes <- vapply(p$layers, function(l) class(l$geom)[1], "")
  expect_true("GeomPoint" %in% layer_classes)
})

test_that("plot_funnel: without labels", {
  skip_if_not_installed("ggplot2")
  sim <- simulate_td(tax_10sp, n_sim = 50, seed = 1)
  obs <- data.frame(site = "A", s = 5, value = 2.0)

  p <- plot_funnel(sim, observed = obs, point_labels = FALSE)
  layer_classes <- vapply(p$layers, function(l) class(l$geom)[1], "")
  expect_false("GeomText" %in% layer_classes)
})

test_that("plot_funnel: errors on wrong sim_result", {
  skip_if_not_installed("ggplot2")
  expect_error(plot_funnel(data.frame(x = 1)), "td_simulation")
})

test_that("plot_funnel: errors on mismatched index", {
  skip_if_not_installed("ggplot2")
  sim <- simulate_td(tax_10sp, n_sim = 50, index = "avtd", seed = 1)
  expect_error(plot_funnel(sim, index = "vartd"), "was computed for")
})

test_that("plot_funnel: errors on malformed observed", {
  skip_if_not_installed("ggplot2")
  sim <- simulate_td(tax_10sp, n_sim = 50, seed = 1)
  expect_error(plot_funnel(sim, observed = data.frame(x = 1)), "must have columns")
  expect_error(plot_funnel(sim, observed = "not a df"), "must be a data frame")
})
