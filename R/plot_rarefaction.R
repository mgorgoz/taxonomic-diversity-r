#' Plot Taxonomic Rarefaction Curve
#'
#' Visualises a rarefaction curve with confidence intervals using ggplot2.
#' Accepts output from [rarefaction_taxonomic()].
#'
#' @param rare_result A data frame returned by [rarefaction_taxonomic()].
#' @param title Optional plot title. If `NULL`, an automatic title is
#'   generated based on the index used.
#' @param xlab Label for the x-axis (default: `"Sample Size (individuals)"`).
#' @param ylab Label for the y-axis. If `NULL`, determined from the index.
#' @param ci_color Fill color for the confidence interval ribbon
#'   (default: `"steelblue"`).
#' @param line_color Color of the mean line (default: `"darkblue"`).
#'
#' @return A ggplot object.
#'
#' @details
#' The plot shows the mean diversity value at each sample size as a solid
#' line, surrounded by a shaded ribbon representing the bootstrap confidence
#' interval. A vertical dashed line marks the total sample size (full data).
#'
#' Requires the \pkg{ggplot2} package.
#'
#' @seealso [rarefaction_taxonomic()] for computing the rarefaction curve.
#'
#' @examples
#' \dontrun{
#' comm <- c(sp1 = 10, sp2 = 5, sp3 = 8, sp4 = 2, sp5 = 3)
#' tax <- data.frame(
#'   Species = paste0("sp", 1:5),
#'   Genus   = c("G1", "G1", "G2", "G2", "G3"),
#'   Family  = c("F1", "F1", "F1", "F2", "F2"),
#'   stringsAsFactors = FALSE
#' )
#' rare <- rarefaction_taxonomic(comm, tax, index = "shannon", n_boot = 50)
#' plot_rarefaction(rare)
#' }
#'
#' @importFrom rlang .data
#' @export
plot_rarefaction <- function(rare_result,
                              title = NULL,
                              xlab = "Sample Size (individuals)",
                              ylab = NULL,
                              ci_color = "steelblue",
                              line_color = "darkblue") {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plot_rarefaction().", call. = FALSE)
  }

  # Extract metadata
  idx <- attr(rare_result, "index")
  total_n <- attr(rare_result, "total_n")
  ci_val <- attr(rare_result, "ci")

  # Index labels
  index_labels <- c(
    species  = "Species Richness",
    shannon  = "Shannon H'",
    simpson  = "Gini-Simpson (1-D)",
    uTO      = "uTO (Unweighted Taxonomic Diversity)",
    TO       = "TO (Weighted Taxonomic Diversity)",
    uTO_plus = "uTO+ (Unweighted Taxonomic Distance)",
    TO_plus  = "TO+ (Weighted Taxonomic Distance)",
    avtd     = "AvTD (Average Taxonomic Distinctness)"
  )

  if (is.null(ylab)) {
    ylab <- if (!is.null(idx) && idx %in% names(index_labels)) {
      index_labels[idx]
    } else {
      "Index Value"
    }
  }

  if (is.null(title)) {
    ci_pct <- if (!is.null(ci_val)) paste0(ci_val * 100, "%") else "95%"
    title <- paste0("Taxonomic Rarefaction Curve - ",
                    if (!is.null(idx)) index_labels[idx] else "Index",
                    " (", ci_pct, " CI)")
  }

  p <- ggplot2::ggplot(rare_result,
                        ggplot2::aes(x = .data$sample_size,
                                     y = .data$mean)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower,
                                       ymax = .data$upper),
                          fill = ci_color, alpha = 0.3) +
    ggplot2::geom_line(color = line_color, linewidth = 1) +
    ggplot2::geom_point(color = line_color, size = 1.5) +
    ggplot2::labs(title = title, x = xlab, y = ylab) +
    ggplot2::theme_minimal(base_size = 12)

  # Add vertical line at total sample size
  if (!is.null(total_n)) {
    p <- p + ggplot2::geom_vline(xintercept = total_n,
                                  linetype = "dashed",
                                  color = "red", alpha = 0.6) +
      ggplot2::annotate("text", x = total_n, y = max(rare_result$mean),
                         label = paste0("N=", total_n),
                         hjust = -0.1, vjust = 1, color = "red", size = 3.5)
  }

  return(p)
}
