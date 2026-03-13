#' Funnel Plot for AvTD/VarTD
#'
#' Produces a Clarke & Warwick style funnel plot showing expected
#' confidence limits for Average Taxonomic Distinctness (AvTD) and/or
#' Variation in Taxonomic Distinctness (VarTD) as a function of
#' species richness. Observed site values can be overlaid to assess
#' whether they fall within or outside the expected range.
#'
#' @param sim_result A `td_simulation` object returned by
#'   [simulate_td()].
#' @param observed Optional data frame with columns `site` (character),
#'   `s` (integer, species richness), and `value` (numeric, observed
#'   AvTD or VarTD). Points are plotted on the funnel.
#' @param index Which index to plot when `sim_result` contains both:
#'   `"avtd"` (default) or `"vartd"`.
#' @param title Optional plot title. If `NULL`, generated automatically.
#' @param point_labels Logical; if `TRUE` (default), label observed
#'   points with site names.
#' @param mean_color Color of the mean line (default: `"darkblue"`).
#' @param ci_color Fill color of the confidence band (default:
#'   `"steelblue"`).
#'
#' @return A ggplot object.
#'
#' @details
#' The funnel shape arises because small samples (low S) have greater
#' random variation in AvTD/VarTD, producing wider confidence bands.
#' As S approaches the full species pool, the band narrows.
#'
#' Observed points falling below the lower confidence limit suggest
#' the community has lower taxonomic breadth than expected by chance,
#' potentially indicating environmental stress or biotic
#' homogenisation.
#'
#' Requires the \pkg{ggplot2} package.
#'
#' @references
#' Clarke, K.R. & Warwick, R.M. (1998). A taxonomic distinctness index
#' and its statistical properties. Journal of Applied Ecology, 35,
#' 523-531.
#'
#' @seealso [simulate_td()] for generating the simulation,
#'   [avtd()] and [vartd()] for the underlying calculations.
#'
#' @examples
#' \donttest{
#' tax <- data.frame(
#'   Species = paste0("sp", 1:10),
#'   Genus   = rep(c("G1", "G2", "G3", "G4", "G5"), each = 2),
#'   Family  = rep(c("F1", "F1", "F2", "F2", "F3"), each = 2),
#'   Order   = rep(c("O1", "O1", "O2", "O2", "O2"), each = 2),
#'   stringsAsFactors = FALSE
#' )
#' sim <- simulate_td(tax, n_sim = 99, seed = 42)
#'
#' # Basic funnel plot
#' plot_funnel(sim)
#'
#' # With observed sites
#' obs <- data.frame(
#'   site  = c("Site_A", "Site_B"),
#'   s     = c(5, 8),
#'   value = c(2.5, 1.8)
#' )
#' plot_funnel(sim, observed = obs)
#' }
#'
#' @importFrom rlang .data
#' @export
plot_funnel <- function(sim_result,
                         observed = NULL,
                         index = c("avtd", "vartd"),
                         title = NULL,
                         point_labels = TRUE,
                         mean_color = "darkblue",
                         ci_color = "steelblue") {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plot_funnel().", call. = FALSE)
  }

  if (!inherits(sim_result, "td_simulation")) {
    stop("'sim_result' must be a 'td_simulation' object from simulate_td().",
         call. = FALSE)
  }

  index <- match.arg(index)
  sim_index <- attr(sim_result, "index")
  ci_val <- attr(sim_result, "ci")

  # Determine which columns to use
  if (sim_index == "both") {
    mean_col  <- paste0("mean_", index)
    lower_col <- paste0("lower_", index)
    upper_col <- paste0("upper_", index)
  } else if (sim_index == index) {
    mean_col  <- paste0("mean_", index)
    lower_col <- paste0("lower_", index)
    upper_col <- paste0("upper_", index)
  } else {
    stop("sim_result was computed for '", sim_index,
         "' but index='", index, "' was requested.", call. = FALSE)
  }

  if (!mean_col %in% names(sim_result)) {
    stop("Column '", mean_col, "' not found in sim_result.", call. = FALSE)
  }

  # Prepare plotting data
  plot_df <- data.frame(
    s     = sim_result$s,
    mean  = sim_result[[mean_col]],
    lower = sim_result[[lower_col]],
    upper = sim_result[[upper_col]]
  )

  # Index labels
  index_labels <- c(
    avtd  = "Average Taxonomic Distinctness",
    vartd = "Variation in Taxonomic Distinctness"
  )
  index_symbols <- c(avtd = expression(Delta^"+"), vartd = expression(Lambda^"+"))

  ylab <- index_labels[index]

  if (is.null(title)) {
    ci_pct <- round(ci_val * 100)
    title <- paste0("Funnel Plot - ", index_labels[index],
                    " (", ci_pct, "% CI)")
  }

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$s)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
      fill = ci_color, alpha = 0.3
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = .data$mean),
      color = mean_color, linewidth = 1
    ) +
    ggplot2::labs(
      title = title,
      x = "Number of Species (S)",
      y = ylab
    ) +
    ggplot2::theme_minimal(base_size = 12)

  # Add observed points
  if (!is.null(observed)) {
    if (!is.data.frame(observed)) {
      stop("'observed' must be a data frame.", call. = FALSE)
    }
    required <- c("site", "s", "value")
    missing <- setdiff(required, names(observed))
    if (length(missing) > 0) {
      stop("'observed' must have columns: ",
           paste(required, collapse = ", "), call. = FALSE)
    }

    p <- p + ggplot2::geom_point(
      data = observed,
      ggplot2::aes(x = .data$s, y = .data$value),
      color = "red", size = 3, shape = 16
    )

    if (point_labels) {
      p <- p + ggplot2::geom_text(
        data = observed,
        ggplot2::aes(x = .data$s, y = .data$value, label = .data$site),
        color = "red", size = 3.5, hjust = -0.15, vjust = -0.5
      )
    }
  }

  return(p)
}
