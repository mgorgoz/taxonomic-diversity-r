#' Calculate Deng Entropy for Taxonomic Diversity
#'
#' Computes the Deng entropy index for a community based on taxonomic
#' hierarchy information. Unlike Shannon entropy which only uses species
#' abundances, Deng entropy incorporates the taxonomic relationships
#' among species through a belief function framework.
#'
#' @param community A named numeric vector of species abundances or a
#'   single-row data frame where columns are species names and values
#'   are abundances.
#' @param tax_tree A data frame representing the taxonomic hierarchy.
#'   Each row is a species, columns represent taxonomic ranks from
#'   lowest (e.g., Species) to highest (e.g., Kingdom). The first
#'   column must match the species names in `community`.
#' @param weights An optional named numeric vector of weights for each

#'   taxonomic level. If NULL (default), equal weights are assigned.
#'
#' @return A numeric value representing the Deng entropy of the community.
#'
#' @details
#' Deng entropy extends Shannon entropy by incorporating taxonomic
#' hierarchy information through Dempster-Shafer evidence theory.
#' The formula accounts for the taxonomic distance between species
#' pairs, giving higher diversity values to communities where species
#' are more taxonomically distinct.
#'
#' The Deng entropy \eqn{E_D} is calculated as:
#' \deqn{E_D = -\sum_{i} m(A_i) \log_2 m(A_i)}
#' where \eqn{m(A_i)} is the basic probability assignment (BPA)
#' derived from taxonomic relationships.
#'
#' @references
#' Deng, Y. (2016). Deng entropy. Chaos, Solitons & Fractals, 91, 549-553.
#'
#' @examples
#' # Create a simple community
#' comm <- c(Quercus_robur = 10, Pinus_nigra = 5, Fagus_orientalis = 8)
#'
#' # Create taxonomic tree
#' tax <- data.frame(
#'   Species = c("Quercus_robur", "Pinus_nigra", "Fagus_orientalis"),
#'   Genus = c("Quercus", "Pinus", "Fagus"),
#'   Family = c("Fagaceae", "Pinaceae", "Fagaceae"),
#'   Order = c("Fagales", "Pinales", "Fagales"),
#'   stringsAsFactors = FALSE
#' )
#'
#' deng_entropy(comm, tax)
#'
#' @export
deng_entropy <- function(community, tax_tree, weights = NULL) {
  # Input validation
  if (!is.numeric(community) || any(community < 0)) {
    stop("'community' must be a numeric vector with non-negative values.",
         call. = FALSE)
  }

  if (sum(community) == 0) {
    stop("'community' must have at least one non-zero abundance.",
         call. = FALSE)
  }

  if (!is.data.frame(tax_tree)) {
    stop("'tax_tree' must be a data frame.", call. = FALSE)
  }

  if (ncol(tax_tree) < 2) {
    stop("'tax_tree' must have at least 2 columns (species + 1 rank).",
         call. = FALSE)
  }

  # Get species names
  species_names <- names(community)
  if (is.null(species_names)) {
    stop("'community' must be a named vector.", call. = FALSE)
  }

  # Filter to non-zero abundances
  community <- community[community > 0]
  species_names <- names(community)

  # Check species are in taxonomy
  tax_species <- as.character(tax_tree[[1]])
  missing_sp <- setdiff(species_names, tax_species)
  if (length(missing_sp) > 0) {
    stop("Species not found in tax_tree: ",
         paste(missing_sp, collapse = ", "), call. = FALSE)
  }

  # If only one species, entropy is 0
  if (length(community) == 1) {
    return(0)
  }

  # Number of taxonomic levels (excluding species column)
  n_levels <- ncol(tax_tree) - 1

  # Set weights
  if (is.null(weights)) {
    weights <- rep(1 / n_levels, n_levels)
  } else {
    if (length(weights) != n_levels) {
      stop("Length of 'weights' must equal number of taxonomic levels (",
           n_levels, ").", call. = FALSE)
    }
    weights <- weights / sum(weights)  # Normalize
  }

  # Compute taxonomic distance matrix
  dist_matrix <- tax_distance_matrix(tax_tree, species_names, weights)

  # Relative abundances (proportions)
  p <- community / sum(community)

  # Compute basic probability assignments (BPA) using Deng framework
  n_sp <- length(community)
  bpa <- numeric(n_sp)

  for (i in seq_len(n_sp)) {
    # BPA for species i considers its proportion weighted by
    # taxonomic similarity to all other species
    sim_sum <- 0
    for (j in seq_len(n_sp)) {
      sim_sum <- sim_sum + p[j] * (1 - dist_matrix[i, j])
    }
    bpa[i] <- p[i] * sim_sum
  }

  # Normalize BPA
  bpa <- bpa / sum(bpa)

  # Calculate Deng entropy
  # Avoid log(0) by filtering zero BPA
  bpa_nonzero <- bpa[bpa > 0]
  entropy <- -sum(bpa_nonzero * log2(bpa_nonzero))

  return(entropy)
}
