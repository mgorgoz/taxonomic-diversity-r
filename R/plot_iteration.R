#' Plot pTO Iteration Results from Run 2 or Run 3
#'
#' Visualizes the iteration-by-iteration pTO values from stochastic
#' resampling (Run 2) or sensitivity analysis (Run 3). Each iteration's
#' value is shown as a point, with the deterministic (Run 1) value
#' displayed as a horizontal reference line.
#'
#' @param resample_result The list returned by
#'   \code{\link{ozkan_pto_resample}} (Run 2) or
#'   \code{\link{ozkan_pto_sensitivity}} (Run 3).
#' @param component Character string specifying which pTO component to plot.
#'   One of \code{"uTO"}, \code{"TO"}, \code{"uTO_plus"}, \code{"TO_plus"}.
#'   Default is \code{"TO"}.
#' @param title Optional character string for the plot title.
#'
#' @return A \code{ggplot} object showing iteration values as points,
#'   the deterministic value as a dashed red line, and the maximum
#'   value as a dashed blue line.
#'
#' @details
#' The plot includes three visual elements:
#' \itemize{
#'   \item \strong{Grey points}: Individual iteration values
#'   \item \strong{Red dashed line}: Deterministic (Run 1) value
#'   \item \strong{Blue dashed line}: Maximum value across all iterations
#' }
#'
#' This helps assess how stochastic species removal affects the pTO
#' index and whether the maximum exceeds the deterministic value.
#'
#' @examples
#' \donttest{
#' comm <- c(sp1 = 10, sp2 = 20, sp3 = 15, sp4 = 5)
#' tax <- build_tax_tree(
#'   species = paste0("sp", 1:4),
#'   Genus = c("G1", "G1", "G2", "G2"),
#'   Family = c("F1", "F1", "F1", "F2")
#' )
#' res <- ozkan_pto_resample(comm, tax, n_iter = 100, seed = 42)
#' plot_iteration(res, component = "TO")
#' }
#'
#' @seealso \code{\link{ozkan_pto_resample}}, \code{\link{ozkan_pto_sensitivity}}
#'
#' @importFrom rlang .data
#' @export
plot_iteration <- function(resample_result,
                           component = "TO",
                           title = NULL) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plot_iteration(). ",
         "Install it with: install.packages(\"ggplot2\")",
         call. = FALSE)
  }

  valid_components <- c("uTO", "TO", "uTO_plus", "TO_plus")
  if (!component %in% valid_components) {
    stop("'component' must be one of: ",
         paste(valid_components, collapse = ", "), ".",
         call. = FALSE)
  }

  if (is.null(resample_result$iteration_results)) {
    stop("'resample_result' must contain 'iteration_results' data frame. ",
         "Use the output of ozkan_pto_resample() or ozkan_pto_sensitivity().",
         call. = FALSE)
  }

  iter_df <- resample_result$iteration_results

  if (!component %in% names(iter_df)) {
    stop("Component '", component, "' not found in iteration_results.",
         call. = FALSE)
  }

  # Get deterministic and max values
  det_key <- paste0(component, "_det")
  max_key <- paste0(component, "_max")

  det_val <- resample_result[[det_key]]
  max_val <- resample_result[[max_key]]

  # If Run 3 result, det values are not directly stored - use iteration 1
  if (is.null(det_val)) {
    det_val <- iter_df[[component]][1]
  }
  if (is.null(max_val)) {
    max_val <- max(iter_df[[component]], na.rm = TRUE)
  }

  # Build iteration index (skip row 1 if it's the deterministic run)
  iter_df$iteration <- seq_len(nrow(iter_df))

  # Component label for axis
  comp_labels <- c(
    uTO = "uTO (Unweighted Taxonomic Diversity)",
    TO = "TO (Weighted Taxonomic Diversity)",
    uTO_plus = "uTO+ (Unweighted Taxonomic Distance)",
    TO_plus = "TO+ (Weighted Taxonomic Distance)"
  )

  p <- ggplot2::ggplot(iter_df,
                       ggplot2::aes(x = .data$iteration,
                                    y = .data[[component]])) +
    ggplot2::geom_point(color = "grey50", alpha = 0.6, size = 1.5) +
    ggplot2::geom_hline(yintercept = det_val,
                        linetype = "dashed", color = "red", linewidth = 0.8) +
    ggplot2::geom_hline(yintercept = max_val,
                        linetype = "dashed", color = "blue", linewidth = 0.8) +
    ggplot2::annotate("text", x = nrow(iter_df) * 0.95, y = det_val,
                      label = paste("Det =", round(det_val, 3)),
                      vjust = -0.5, hjust = 1, color = "red", size = 3.5) +
    ggplot2::annotate("text", x = nrow(iter_df) * 0.95, y = max_val,
                      label = paste("Max =", round(max_val, 3)),
                      vjust = 1.5, hjust = 1, color = "blue", size = 3.5) +
    ggplot2::labs(
      x = "Iteration",
      y = comp_labels[component],
      title = title %||% paste("Stochastic Resampling -", component)
    ) +
    ggplot2::theme_minimal()

  return(p)
}
