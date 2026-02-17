#' Compute Taxonomic Distance Matrix
#'
#' Calculates pairwise taxonomic distances between species based on
#' their positions in a taxonomic hierarchy. Distance is computed as
#' the weighted proportion of taxonomic levels at which two species
#' diverge.
#'
#' @param tax_tree A data frame representing the taxonomic hierarchy.
#'   First column must be species names, subsequent columns are
#'   taxonomic ranks from lowest to highest.
#' @param species Optional character vector of species names to include.
#'   If NULL, all species in `tax_tree` are used.
#' @param weights Optional numeric vector of weights for each taxonomic
#'   level. If NULL, equal weights are assigned.
#'
#' @return A symmetric matrix of taxonomic distances between species.
#'   With default equal step weights (1, 2, 3, ...), values range from
#'   0 (same species) to the number of taxonomic levels (maximum
#'   distance when no common ancestor is found at any level).
#'
#' @examples
#' tax <- data.frame(
#'   Species = c("Quercus_robur", "Pinus_nigra", "Fagus_orientalis"),
#'   Genus = c("Quercus", "Pinus", "Fagus"),
#'   Family = c("Fagaceae", "Pinaceae", "Fagaceae"),
#'   Order = c("Fagales", "Pinales", "Fagales"),
#'   stringsAsFactors = FALSE
#' )
#'
#' tax_distance_matrix(tax)
#'
#' @export
tax_distance_matrix <- function(tax_tree, species = NULL, weights = NULL) {
  if (!is.data.frame(tax_tree)) {
    stop("'tax_tree' must be a data frame.", call. = FALSE)
  }

  if (ncol(tax_tree) < 2) {
    stop("'tax_tree' must have at least 2 columns.", call. = FALSE)
  }

  # Get species list
  all_species <- as.character(tax_tree[[1]])


  if (is.null(species)) {
    species <- all_species
  }

  # Subset taxonomy to requested species
  idx <- match(species, all_species)
  if (any(is.na(idx))) {
    stop("Species not found in tax_tree: ",
         paste(species[is.na(idx)], collapse = ", "), call. = FALSE)
  }
  tax_sub <- tax_tree[idx, , drop = FALSE]

  n_sp <- length(species)
  n_levels <- ncol(tax_tree) - 1

  # Set weights (equal step lengths by default: 1, 2, 3, ...)
  if (is.null(weights)) {
    weights <- seq_len(n_levels)
  }
  if (length(weights) != n_levels) {
    stop("'weights' must have length ", n_levels,
         " (one per taxonomic level), got ", length(weights), ".",
         call. = FALSE)
  }

  # Initialize distance matrix
  dist_mat <- matrix(0, nrow = n_sp, ncol = n_sp,
                     dimnames = list(species, species))

  # Compute pairwise distances using Clarke & Warwick path length:
  # Find the first MATCHING taxonomic level from bottom (lowest rank).
  # Distance = weight of that level. If no match found, use max weight.
  max_weight <- weights[n_levels]  # consistent with delta()/delta_star()
  for (i in seq_len(n_sp - 1)) {
    for (j in (i + 1):n_sp) {
      d <- max_weight  # max distance if nothing matches
      for (k in seq_len(n_levels)) {
        if (as.character(tax_sub[i, k + 1]) ==
            as.character(tax_sub[j, k + 1])) {
          d <- weights[k]
          break
        }
      }
      dist_mat[i, j] <- d
      dist_mat[j, i] <- d
    }
  }

  return(dist_mat)
}


#' Build a Taxonomic Tree from Species Data
#'
#' Creates a taxonomic hierarchy data frame from species classification
#' information. This is a convenience function for constructing the
#' `tax_tree` input required by other functions in the package.
#'
#' @param species Character vector of species names.
#' @param ... Named character vectors for each taxonomic rank, in order
#'   from lowest to highest (e.g., Genus, Family, Order).
#'
#' @return A data frame with species as the first column and taxonomic
#'   ranks as subsequent columns.
#'
#' @examples
#' tree <- build_tax_tree(
#'   species = c("Quercus_robur", "Pinus_nigra", "Fagus_orientalis"),
#'   Genus   = c("Quercus", "Pinus", "Fagus"),
#'   Family  = c("Fagaceae", "Pinaceae", "Fagaceae"),
#'   Order   = c("Fagales", "Pinales", "Fagales")
#' )
#'
#' @export
build_tax_tree <- function(species, ...) {
  ranks <- list(...)

  if (length(ranks) == 0) {
    stop("At least one taxonomic rank must be provided.", call. = FALSE)
  }

  # Check all ranks have same length as species
  n <- length(species)
  for (nm in names(ranks)) {
    if (length(ranks[[nm]]) != n) {
      stop("Length of '", nm, "' (", length(ranks[[nm]]),
           ") does not match number of species (", n, ").",
           call. = FALSE)
    }
  }

  df <- data.frame(Species = species, stringsAsFactors = FALSE)
  for (nm in names(ranks)) {
    df[[nm]] <- as.character(ranks[[nm]])
  }

  return(df)
}
