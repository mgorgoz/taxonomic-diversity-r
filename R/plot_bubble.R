#' Bubble Chart of Species Contributions to Diversity
#'
#' Creates a bubble chart showing each species' abundance (x-axis),
#' average taxonomic distance to other species (y-axis), and relative
#' contribution to the community (bubble size). Species that are both
#' abundant and taxonomically distant from others contribute most to
#' overall taxonomic diversity.
#'
#' @param community Named numeric vector of species abundances.
#' @param tax_tree A data frame representing the taxonomic hierarchy,
#'   as produced by \code{\link{build_tax_tree}}.
#' @param color_by Character string specifying which taxonomic rank to use
#'   for coloring bubbles. Must match a column name in \code{tax_tree}.
#'   If \code{NULL} (default), the highest available rank is used.
#' @param title Optional character string for the plot title.
#'
#' @return A \code{ggplot} object.
#'
#' @details
#' For each species \eqn{i}, the average taxonomic distance is calculated as:
#'
#' \deqn{\bar{\omega}_i = \frac{1}{S-1} \sum_{j \neq i} \omega_{ij}}
#'
#' where \eqn{\omega_{ij}} is the pairwise taxonomic distance and \eqn{S}
#' is the number of species. Bubble size represents the product of
#' relative abundance and average distance, indicating each species'
#' contribution to overall taxonomic diversity.
#'
#' @examples
#' \dontrun{
#' comm <- c(sp1 = 25, sp2 = 18, sp3 = 30, sp4 = 12, sp5 = 8)
#' tax <- build_tax_tree(
#'   species = paste0("sp", 1:5),
#'   Genus   = c("G1", "G1", "G2", "G2", "G3"),
#'   Family  = c("F1", "F1", "F1", "F2", "F2"),
#'   Order   = c("O1", "O1", "O1", "O1", "O1")
#' )
#' plot_bubble(comm, tax)
#' }
#'
#' @seealso \code{\link{tax_distance_matrix}}, \code{\link{compare_indices}}
#'
#' @importFrom rlang .data
#' @export
plot_bubble <- function(community, tax_tree,
                        color_by = NULL,
                        title = NULL) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plot_bubble(). ",
         "Install it with: install.packages(\"ggplot2\")",
         call. = FALSE)
  }

  if (!is.numeric(community) || is.null(names(community))) {
    stop("'community' must be a named numeric vector.", call. = FALSE)
  }

  if (!is.data.frame(tax_tree)) {
    stop("'tax_tree' must be a data frame.", call. = FALSE)
  }

  species <- names(community[community > 0])
  n_sp <- length(species)

  if (n_sp < 2) {
    stop("At least 2 species with positive abundance are required.",
         call. = FALSE)
  }

  # Compute distance matrix
  dm <- tax_distance_matrix(tax_tree, species = species)

  # Average taxonomic distance for each species
  avg_dist <- vapply(seq_len(n_sp), function(i) {
    sum(dm[i, -i]) / (n_sp - 1)
  }, numeric(1))

  # Relative abundance
  abund <- community[species]
  rel_abund <- abund / sum(abund)

  # Contribution = relative abundance * average distance
  contribution <- unname(rel_abund * avg_dist)

  # Color grouping
  rank_names <- names(tax_tree)[-1]
  if (is.null(color_by)) {
    color_by <- rank_names[length(rank_names)]
  }
  if (!color_by %in% rank_names) {
    stop("'color_by' must be one of: ",
         paste(rank_names, collapse = ", "), ".",
         call. = FALSE)
  }

  sp_idx <- match(species, as.character(tax_tree[[1]]))
  groups <- as.character(tax_tree[[color_by]][sp_idx])

  # Build data frame
  bubble_df <- data.frame(
    Species      = species,
    Abundance    = unname(abund),
    Avg_Distance = avg_dist,
    Contribution = contribution,
    Group        = groups,
    stringsAsFactors = FALSE
  )

  p <- ggplot2::ggplot(bubble_df,
                       ggplot2::aes(x = .data$Abundance,
                                    y = .data$Avg_Distance,
                                    size = .data$Contribution,
                                    color = .data$Group)) +
    ggplot2::geom_point(alpha = 0.7) +
    ggplot2::geom_text(
      ggplot2::aes(label = .data$Species),
      size = 3, vjust = -1.2, show.legend = FALSE
    ) +
    ggplot2::scale_size_continuous(
      name = "Contribution",
      range = c(3, 15)
    ) +
    ggplot2::labs(
      x = "Abundance",
      y = "Average Taxonomic Distance",
      color = color_by,
      title = title %||% "Species Contributions to Taxonomic Diversity"
    ) +
    ggplot2::theme_minimal()

  return(p)
}
