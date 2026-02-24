#' Radar (Spider) Chart for Multi-Community Index Comparison
#'
#' Creates a radar chart comparing diversity indices across multiple
#' communities. Each axis represents a different index, and each
#' community is drawn as a colored polygon. Values are normalized
#' to 0-1 scale so that indices with different magnitudes can be
#' compared visually.
#'
#' @param communities A named list of community vectors (named numeric).
#' @param tax_tree A data frame representing the taxonomic hierarchy,
#'   as produced by \code{\link{build_tax_tree}}.
#' @param indices Character vector specifying which indices to include.
#'   Default is all 10 indices. Available: \code{"Shannon"}, \code{"Simpson"},
#'   \code{"Delta"}, \code{"Delta_star"}, \code{"AvTD"}, \code{"VarTD"},
#'   \code{"uTO"}, \code{"TO"}, \code{"uTO_plus"}, \code{"TO_plus"}.
#' @param title Optional character string for the plot title.
#'
#' @return A \code{ggplot} object.
#'
#' @details
#' Each index value is normalized using min-max scaling across the
#' communities being compared:
#'
#' \deqn{x_{norm} = \frac{x - x_{min}}{x_{max} - x_{min}}}
#'
#' If all communities have the same value for an index (i.e.,
#' \eqn{x_{max} = x_{min}}), the normalized value is set to 0.5.
#'
#' The radar chart is built using polar coordinates in ggplot2.
#' Each community appears as a colored polygon overlay, making it
#' easy to spot which community scores higher on which indices.
#'
#' @examples
#' \donttest{
#' tax <- build_tax_tree(
#'   species = c("sp1", "sp2", "sp3", "sp4"),
#'   Genus   = c("G1", "G1", "G2", "G2"),
#'   Family  = c("F1", "F1", "F1", "F2")
#' )
#' comms <- list(
#'   Site_A = c(sp1 = 10, sp2 = 20, sp3 = 15, sp4 = 5),
#'   Site_B = c(sp1 = 5, sp2 = 5, sp3 = 5, sp4 = 5)
#' )
#' plot_radar(comms, tax)
#' }
#'
#' @seealso \code{\link{compare_indices}} for tabular comparison
#'
#' @importFrom rlang .data
#' @export
plot_radar <- function(communities, tax_tree,
                       indices = NULL,
                       title = NULL) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plot_radar(). ",
         "Install it with: install.packages(\"ggplot2\")",
         call. = FALSE)
  }

  # Get comparison table
  comp_df <- compare_indices(communities, tax_tree, plot = FALSE)

  all_indices <- c("Shannon", "Simpson", "Delta", "Delta_star",
                   "AvTD", "VarTD", "uTO", "TO", "uTO_plus", "TO_plus")

  if (is.null(indices)) {
    indices <- all_indices
  }

  bad_idx <- setdiff(indices, all_indices)
  if (length(bad_idx) > 0) {
    stop("Unknown indices: ", paste(bad_idx, collapse = ", "),
         ". Available: ", paste(all_indices, collapse = ", "),
         call. = FALSE)
  }

  # Extract index columns and normalize to 0-1
  n_communities <- nrow(comp_df)
  n_indices <- length(indices)

  norm_vals <- as.data.frame(lapply(indices, function(idx) {
    vals <- comp_df[[idx]]
    rng <- range(vals, na.rm = TRUE)
    if (rng[2] == rng[1]) {
      rep(0.5, length(vals))
    } else {
      (vals - rng[1]) / (rng[2] - rng[1])
    }
  }))
  names(norm_vals) <- indices

  # Build long-format data for radar
  # Each community needs points for all indices + close the polygon
  radar_list <- list()
  for (i in seq_len(n_communities)) {
    comm_name <- comp_df$Community[i]
    vals <- as.numeric(norm_vals[i, ])

    # Close the polygon by repeating the first point
    vals_closed <- c(vals, vals[1])
    idx_closed <- c(indices, indices[1])
    angle_closed <- seq_len(n_indices + 1)

    radar_list[[i]] <- data.frame(
      Community = comm_name,
      Index     = idx_closed,
      Value     = vals_closed,
      Angle     = angle_closed,
      stringsAsFactors = FALSE
    )
  }
  radar_df <- do.call(rbind, radar_list)

  # Factor ordering for indices
  radar_df$Index <- factor(radar_df$Index, levels = indices)

  p <- ggplot2::ggplot(radar_df,
                       ggplot2::aes(x = .data$Angle, y = .data$Value,
                                    color = .data$Community,
                                    fill = .data$Community,
                                    group = .data$Community)) +
    ggplot2::geom_polygon(alpha = 0.15, linewidth = 0.8) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_x_continuous(
      breaks = seq_len(n_indices),
      labels = indices,
      limits = c(0.5, n_indices + 0.5)
    ) +
    ggplot2::scale_y_continuous(
      limits = c(-0.1, 1.05),
      breaks = c(0, 0.25, 0.5, 0.75, 1.0)
    ) +
    ggplot2::coord_polar() +
    ggplot2::labs(
      title = title %||% "Diversity Indices Radar Chart",
      y = "Normalized Value",
      x = NULL,
      color = "Community",
      fill = "Community"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 9),
      panel.grid.major = ggplot2::element_line(color = "grey85"),
      legend.position = "top"
    )

  return(p)
}
