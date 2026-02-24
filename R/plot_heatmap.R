#' Plot Taxonomic Distance Heatmap
#'
#' Visualizes the pairwise taxonomic distance matrix as a colored heatmap
#' using ggplot2. Species are ordered by hierarchical clustering so that
#' taxonomically similar species appear adjacent.
#'
#' @param tax_tree A data frame representing the taxonomic hierarchy,
#'   as produced by \code{\link{build_tax_tree}}.
#' @param label_size Numeric value controlling the size of species labels.
#'   Default is 3.
#' @param title Optional character string for the plot title.
#' @param low_color Color for the lowest distance (most similar).
#'   Default is \code{"white"}.
#' @param high_color Color for the highest distance (most distant).
#'   Default is \code{"#B22222"} (firebrick red).
#'
#' @return A \code{ggplot} object.
#'
#' @details
#' The heatmap displays the full symmetric distance matrix computed by
#' \code{\link{tax_distance_matrix}}. The diagonal (self-distance = 0)
#' appears in the lowest color. Species are reordered using hierarchical
#' clustering (UPGMA) to reveal taxonomic groupings visually.
#'
#' @examples
#' \donttest{
#' tax <- build_tax_tree(
#'   species = c("sp1", "sp2", "sp3", "sp4"),
#'   Genus   = c("G1", "G1", "G2", "G2"),
#'   Family  = c("F1", "F1", "F1", "F2")
#' )
#' plot_heatmap(tax)
#' }
#'
#' @seealso \code{\link{tax_distance_matrix}}, \code{\link{plot_taxonomic_tree}}
#'
#' @importFrom rlang .data
#' @export
plot_heatmap <- function(tax_tree,
                         label_size = 3,
                         title = NULL,
                         low_color = "white",
                         high_color = "#B22222") {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plot_heatmap(). ",
         "Install it with: install.packages(\"ggplot2\")",
         call. = FALSE)
  }

  if (!is.data.frame(tax_tree)) {
    stop("'tax_tree' must be a data frame.", call. = FALSE)
  }

  # Compute distance matrix
  dm <- tax_distance_matrix(tax_tree)

  # Reorder by hierarchical clustering
  hc <- stats::hclust(stats::as.dist(dm), method = "average")
  ord <- hc$labels[hc$order]

  # Reshape to long format
  species <- rownames(dm)
  long_df <- data.frame(
    Species1 = rep(species, each = length(species)),
    Species2 = rep(species, times = length(species)),
    Distance = as.vector(dm),
    stringsAsFactors = FALSE
  )

  # Apply clustering order
  long_df$Species1 <- factor(long_df$Species1, levels = ord)
  long_df$Species2 <- factor(long_df$Species2, levels = rev(ord))

  p <- ggplot2::ggplot(long_df,
                       ggplot2::aes(x = .data$Species1, y = .data$Species2,
                                    fill = .data$Distance)) +
    ggplot2::geom_tile(color = "grey90", linewidth = 0.3) +
    ggplot2::geom_text(
      ggplot2::aes(label = round(.data$Distance, 1)),
      size = label_size, color = "grey30"
    ) +
    ggplot2::scale_fill_gradient(
      low = low_color, high = high_color,
      name = "Taxonomic\nDistance"
    ) +
    ggplot2::labs(x = NULL, y = NULL, title = title) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1,
                                           size = label_size * 3),
      axis.text.y = ggplot2::element_text(size = label_size * 3),
      panel.grid = ggplot2::element_blank()
    ) +
    ggplot2::coord_fixed()

  return(p)
}
