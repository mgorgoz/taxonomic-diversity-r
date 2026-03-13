#' Simulate Expected AvTD/VarTD Under Random Sampling
#'
#' Generates the null distribution of Average Taxonomic Distinctness
#' (AvTD) and/or Variation in Taxonomic Distinctness (VarTD) by
#' randomly drawing species subsets from a regional species pool. Used
#' to construct funnel plots for statistical testing (Clarke & Warwick
#' 1998, 2001).
#'
#' @param tax_tree A data frame representing the full regional species
#'   pool taxonomy. First column is species names, subsequent columns
#'   are taxonomic ranks from lowest to highest.
#' @param s_range Integer vector of species richness values to
#'   simulate. Default `NULL` uses `2:S` where S is the total number
#'   of species in `tax_tree`.
#' @param n_sim Number of random draws per species richness value
#'   (default 999).
#' @param index Which index to simulate: `"avtd"`, `"vartd"`, or
#'   `"both"` (default).
#' @param weights Optional numeric vector of weights for taxonomic
#'   levels. Passed to [avtd()] and [vartd()].
#' @param ci Confidence interval width (default 0.95).
#' @param seed Optional random seed for reproducibility.
#'
#' @return A data frame with class `"td_simulation"` containing
#'   columns:
#'   \describe{
#'     \item{s}{Species richness (number of species drawn)}
#'     \item{mean_avtd}{Mean simulated AvTD (if index includes avtd)}
#'     \item{lower_avtd}{Lower CI bound for AvTD}
#'     \item{upper_avtd}{Upper CI bound for AvTD}
#'     \item{mean_vartd}{Mean simulated VarTD (if index includes vartd)}
#'     \item{lower_vartd}{Lower CI bound for VarTD}
#'     \item{upper_vartd}{Upper CI bound for VarTD}
#'   }
#'
#'   Attributes: `ci`, `index`, `n_sim`, `pool_size`.
#'
#' @details
#' For each value of S in `s_range`, `n_sim` random subsets of S
#' species are drawn (without replacement) from the full species pool
#' in `tax_tree`. AvTD and/or VarTD are computed for each random
#' subset. The mean and percentile-based confidence limits are
#' recorded.
#'
#' The resulting object can be passed to [plot_funnel()] to produce
#' the classic Clarke & Warwick funnel plot.
#'
#' @references
#' Clarke, K.R. & Warwick, R.M. (1998). A taxonomic distinctness index
#' and its statistical properties. Journal of Applied Ecology, 35,
#' 523-531.
#'
#' Clarke, K.R. & Warwick, R.M. (2001). A further biodiversity index
#' applicable to species lists: variation in taxonomic distinctness.
#' Marine Ecology Progress Series, 216, 265-278.
#'
#' @seealso [plot_funnel()] for visualisation, [avtd()] and [vartd()]
#'   for the underlying calculations.
#'
#' @examples
#' tax <- data.frame(
#'   Species = paste0("sp", 1:10),
#'   Genus   = rep(c("G1", "G2", "G3", "G4", "G5"), each = 2),
#'   Family  = rep(c("F1", "F1", "F2", "F2", "F3"), each = 2),
#'   Order   = rep(c("O1", "O1", "O2", "O2", "O2"), each = 2),
#'   stringsAsFactors = FALSE
#' )
#' sim <- simulate_td(tax, n_sim = 99, seed = 42)
#' sim
#'
#' @export
simulate_td <- function(tax_tree,
                         s_range = NULL,
                         n_sim = 999L,
                         index = c("both", "avtd", "vartd"),
                         weights = NULL,
                         ci = 0.95,
                         seed = NULL) {

  index <- match.arg(index)

  if (!is.data.frame(tax_tree)) {
    stop("'tax_tree' must be a data frame.", call. = FALSE)
  }
  if (ncol(tax_tree) < 2) {
    stop("'tax_tree' must have at least 2 columns.", call. = FALSE)
  }

  all_species <- as.character(tax_tree[[1]])
  pool_size <- length(all_species)

  if (pool_size < 3) {
    stop("Species pool must contain at least 3 species.", call. = FALSE)
  }

  if (is.null(s_range)) {
    s_range <- seq(2L, pool_size)
  }
  s_range <- as.integer(s_range)
  s_range <- s_range[s_range >= 2L & s_range <= pool_size]

  if (length(s_range) == 0) {
    stop("'s_range' must contain values between 2 and ", pool_size, ".",
         call. = FALSE)
  }

  n_sim <- as.integer(n_sim)
  if (n_sim < 1L) {
    stop("'n_sim' must be a positive integer.", call. = FALSE)
  }

  if (ci <= 0 || ci >= 1) {
    stop("'ci' must be between 0 and 1 (exclusive).", call. = FALSE)
  }

  if (!is.null(seed)) set.seed(seed)

  do_avtd <- index %in% c("avtd", "both")
  do_vartd <- index %in% c("vartd", "both")

  alpha <- 1 - ci
  q_lo <- alpha / 2
  q_hi <- 1 - alpha / 2

  results <- vector("list", length(s_range))

  for (k in seq_along(s_range)) {
    s <- s_range[k]

    avtd_vals <- if (do_avtd) numeric(n_sim) else NULL
    vartd_vals <- if (do_vartd) numeric(n_sim) else NULL

    for (i in seq_len(n_sim)) {
      spp <- sample(all_species, s, replace = FALSE)

      if (do_avtd) {
        avtd_vals[i] <- avtd(spp, tax_tree, weights = weights)
      }
      if (do_vartd) {
        vartd_vals[i] <- vartd(spp, tax_tree, weights = weights)
      }
    }

    row <- list(s = s)

    if (do_avtd) {
      row$mean_avtd  <- mean(avtd_vals)
      row$lower_avtd <- stats::quantile(avtd_vals, q_lo, names = FALSE)
      row$upper_avtd <- stats::quantile(avtd_vals, q_hi, names = FALSE)
    }
    if (do_vartd) {
      row$mean_vartd  <- mean(vartd_vals)
      row$lower_vartd <- stats::quantile(vartd_vals, q_lo, names = FALSE)
      row$upper_vartd <- stats::quantile(vartd_vals, q_hi, names = FALSE)
    }

    results[[k]] <- as.data.frame(row)
  }

  out <- do.call(rbind, results)
  rownames(out) <- NULL

  attr(out, "ci") <- ci
  attr(out, "index") <- index
  attr(out, "n_sim") <- n_sim
  attr(out, "pool_size") <- pool_size
  class(out) <- c("td_simulation", "data.frame")

  out
}


#' @export
print.td_simulation <- function(x, ...) {
  idx <- attr(x, "index")
  ci_pct <- attr(x, "ci") * 100
  n_sim <- attr(x, "n_sim")
  pool <- attr(x, "pool_size")

  cat("Taxonomic Distinctness Simulation\n")
  cat("  Index:", idx, "\n")
  cat("  Species pool:", pool, "species\n")
  cat("  Simulations per S:", n_sim, "\n")
  cat("  Confidence interval:", ci_pct, "%\n")
  cat("  S range:", min(x$s), "-", max(x$s), "\n\n")

  print.data.frame(utils::head(x, 10), row.names = FALSE)
  if (nrow(x) > 10) {
    cat("... (", nrow(x), " rows total)\n")
  }
  invisible(x)
}
