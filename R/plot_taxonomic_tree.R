#' Plot Taxonomic Tree as a Dendrogram
#'
#' Visualizes a taxonomic hierarchy as a dendrogram (tree diagram) using
#' ggplot2. The function converts the taxonomic distance matrix into a
#' hierarchical clustering object and renders it as a horizontal dendrogram
#' with species labels colored by a chosen taxonomic rank.
#'
#' The dendrogram is constructed from the pairwise taxonomic distance matrix
#' (computed via \code{\link{tax_distance_matrix}}) using hierarchical
#' clustering (\code{\link[stats]{hclust}} with \code{method = "average"}).
#' Branch heights reflect taxonomic distance: species in the same genus
#' merge at the lowest level, while species in different orders merge at
#' the highest level.
#'
#' When a \code{community} vector is provided, species labels are annotated
#' with abundance values in parentheses, e.g., "Quercus_coccifera (25)".
#'
#' @param tax_tree A data frame representing the taxonomic hierarchy,
#'   as produced by \code{\link{build_tax_tree}}. First column must be
#'   species names, subsequent columns are taxonomic ranks from lowest
#'   to highest (e.g., Genus, Family, Order).
#' @param community Optional named numeric vector of species abundances.
#'   Names must match species in \code{tax_tree}. When provided, abundance
#'   values are shown next to species labels.
#' @param color_by Character string specifying which taxonomic rank to use
#'   for coloring species labels. Must match a column name in \code{tax_tree}.
#'   If \code{NULL} (default), the highest available rank is used.
#' @param label_size Numeric value controlling the size of species labels.
#'   Default is 3.
#' @param title Optional character string for the plot title. If \code{NULL}
#'   (default), no title is displayed.
#'
#' @return A \code{ggplot} object that can be further customized with
#'   ggplot2 functions (e.g., \code{+ theme()}, \code{+ labs()}).
#'
#' @details
#' This function requires the \pkg{ggplot2} package. If ggplot2 is not
#' installed, the function will stop with an informative error message.
#'
#' The clustering method used is UPGMA (Unweighted Pair Group Method
#' with Arithmetic Mean), which is standard for taxonomic classification
#' trees where branch lengths represent average distances between groups.
#'
#' @examples
#' \dontrun{
#' # Build a simple taxonomic tree
#' tax <- build_tax_tree(
#'   species = c("Quercus_robur", "Quercus_petraea", "Pinus_nigra",
#'               "Pinus_brutia", "Juniperus_excelsa"),
#'   Genus   = c("Quercus", "Quercus", "Pinus", "Pinus", "Juniperus"),
#'   Family  = c("Fagaceae", "Fagaceae", "Pinaceae", "Pinaceae", "Cupressaceae"),
#'   Order   = c("Fagales", "Fagales", "Pinales", "Pinales", "Pinales")
#' )
#'
#' # Basic dendrogram
#' plot_taxonomic_tree(tax)
#'
#' # Color by Family, with abundance info
#' comm <- c(Quercus_robur = 25, Quercus_petraea = 18,
#'           Pinus_nigra = 30, Pinus_brutia = 12,
#'           Juniperus_excelsa = 8)
#' plot_taxonomic_tree(tax, community = comm, color_by = "Family")
#' }
#'
#' @seealso \code{\link{build_tax_tree}} for creating the taxonomy input,
#'   \code{\link{tax_distance_matrix}} for the underlying distance calculation.
#'
#' @references
#' Clarke, K.R. & Warwick, R.M. (1998). A taxonomic distinctness index
#' and its statistical properties. \emph{Journal of Applied Ecology}, 35,
#' 523--531.
#'
#' @importFrom rlang .data
#' @export
plot_taxonomic_tree <- function(tax_tree,
                                community = NULL,
                                color_by = NULL,
                                label_size = 3,
                                title = NULL) {


  # --- Check ggplot2 availability ---
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plot_taxonomic_tree(). ",
         "Install it with: install.packages(\"ggplot2\")",
         call. = FALSE)
  }

  # --- Input validation ---
  if (!is.data.frame(tax_tree)) {
    stop("'tax_tree' must be a data frame.", call. = FALSE)
  }
  if (ncol(tax_tree) < 2) {
    stop("'tax_tree' must have at least 2 columns (Species + taxonomic ranks).",
         call. = FALSE)
  }

  species <- as.character(tax_tree[[1]])
  n_sp <- length(species)
  rank_names <- names(tax_tree)[-1]

  if (n_sp < 2) {
    stop("At least 2 species are required to plot a dendrogram.",
         call. = FALSE)
  }

  # --- Validate color_by ---
  if (is.null(color_by)) {
    color_by <- rank_names[length(rank_names)]  # highest rank
  }
  if (!color_by %in% rank_names) {
    stop("'color_by' must be one of: ",
         paste(rank_names, collapse = ", "), ".",
         call. = FALSE)
  }

  # --- Validate community ---
  if (!is.null(community)) {
    if (!is.numeric(community) || is.null(names(community))) {
      stop("'community' must be a named numeric vector.", call. = FALSE)
    }
    missing_sp <- setdiff(species, names(community))
    if (length(missing_sp) > 0) {
      stop("Species not found in 'community': ",
           paste(missing_sp, collapse = ", "), call. = FALSE)
    }
  }

  # --- Compute distance matrix and hierarchical clustering ---
  dist_mat <- tax_distance_matrix(tax_tree)
  hc <- stats::hclust(stats::as.dist(dist_mat), method = "average")

  # --- Build dendrogram segment data ---
  dend <- stats::as.dendrogram(hc)
  dend_data <- dendro_to_segments(dend)

  # --- Build label data (leaf positions) ---
  leaf_order <- hc$order
  leaf_labels <- hc$labels[leaf_order]

  # Create label data frame
  label_df <- data.frame(
    x = seq_along(leaf_labels),
    y = 0,
    label = leaf_labels,
    stringsAsFactors = FALSE
  )

  # Add color group from taxonomy
  color_idx <- match(label_df$label, species)
  label_df$group <- as.character(tax_tree[[color_by]][color_idx])

  # Add abundance to labels if community is provided
  if (!is.null(community)) {
    abundances <- community[label_df$label]
    label_df$label <- paste0(label_df$label, " (", abundances, ")")
  }

  # --- Build ggplot ---
  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = dend_data,
      ggplot2::aes(x = .data$x, y = .data$y,
                   xend = .data$xend, yend = .data$yend),
      color = "grey40", linewidth = 0.5
    ) +
    ggplot2::geom_text(
      data = label_df,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$label,
                   color = .data$group),
      hjust = 1, size = label_size, nudge_y = -0.05
    ) +
    ggplot2::coord_flip() +
    ggplot2::scale_y_continuous(
      name = "Taxonomic distance",
      limits = c(min(-0.5, -max(dend_data$y) * 0.3), max(dend_data$y) * 1.05)
    ) +
    ggplot2::scale_x_continuous(breaks = NULL) +
    ggplot2::labs(
      x = NULL,
      color = color_by,
      title = title
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = 9)
    )

  return(p)
}


# =========================================================================
# Internal helper: Convert dendrogram to segment data frame
# =========================================================================
#
# Recursively traverses a stats::dendrogram object and extracts the
# (x, y) -> (xend, yend) line segments needed to draw horizontal and
# vertical branches with ggplot2::geom_segment().
#
# @param dend A dendrogram object from stats::as.dendrogram().
# @return A data frame with columns: x, y, xend, yend.
# @keywords internal
dendro_to_segments <- function(dend) {
  segments <- list()
  counter <- new.env(parent = emptyenv())
  counter$leaf <- 0

  get_pos <- function(node) {
    if (stats::is.leaf(node)) {
      counter$leaf <- counter$leaf + 1
      return(list(x = counter$leaf, y = 0))
    }

    height <- attr(node, "height")
    children <- lapply(node, get_pos)

    # Horizontal segments from each child to parent height
    for (child in children) {
      segments[[length(segments) + 1]] <<- data.frame(
        x = child$x, y = child$y,
        xend = child$x, yend = height
      )
    }

    # Vertical segment connecting children at parent height
    child_xs <- vapply(children, function(ch) ch$x, numeric(1))
    segments[[length(segments) + 1]] <<- data.frame(
      x = min(child_xs), y = height,
      xend = max(child_xs), yend = height
    )

    # Return midpoint position
    return(list(x = mean(child_xs), y = height))
  }

  get_pos(dend)

  do.call(rbind, segments)
}
