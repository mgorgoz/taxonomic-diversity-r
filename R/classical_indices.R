#' Shannon Diversity Index
#'
#' Calculates the Shannon-Wiener diversity index (H') for a community.
#'
#' @param community A numeric vector of species abundances.
#' @param base The logarithm base. Default is `exp(1)` (natural log).
#'   Use `2` for bits.
#'
#' @return A numeric value representing the Shannon diversity index.
#'
#' @details
#' The Shannon index is calculated as:
#' \deqn{H' = -\sum_{i=1}^{S} p_i \ln(p_i)}
#' where \eqn{p_i} is the proportion of species \eqn{i} and \eqn{S}
#' is the total number of species.
#'
#' @examples
#' comm <- c(10, 5, 8, 3, 12)
#' shannon(comm)
#'
#' @export
shannon <- function(community, base = exp(1)) {
  if (!is.numeric(community) || any(community < 0)) {
    stop("'community' must be a numeric vector with non-negative values.",
         call. = FALSE)
  }

  # Remove zeros
  community <- community[community > 0]

  if (length(community) == 0) {
    return(0)
  }

  p <- community / sum(community)
  -sum(p * log(p, base = base))
}


#' Simpson Diversity Index
#'
#' Calculates the Simpson diversity index (1 - D) for a community.
#'
#' @param community A numeric vector of species abundances.
#' @param type One of `"inverse"` (1/D), `"gini_simpson"` (1 - D),
#'   or `"dominance"` (D). Default is `"gini_simpson"`.
#'
#' @return A numeric value representing the Simpson index.
#'
#' @details
#' Simpson's dominance index D is calculated as:
#' \deqn{D = \sum_{i=1}^{S} p_i^2}
#' The Gini-Simpson index is \eqn{1 - D} and the inverse Simpson is
#' \eqn{1/D}.
#'
#' @examples
#' comm <- c(10, 5, 8, 3, 12)
#' simpson(comm)
#' simpson(comm, type = "inverse")
#'
#' @export
simpson <- function(community, type = c("gini_simpson", "inverse",
                                        "dominance")) {
  type <- match.arg(type)

  if (!is.numeric(community) || any(community < 0)) {
    stop("'community' must be a numeric vector with non-negative values.",
         call. = FALSE)
  }

  community <- community[community > 0]

  if (length(community) == 0) {
    return(0)
  }

  p <- community / sum(community)
  D <- sum(p^2)

  switch(type,
    dominance = D,
    gini_simpson = 1 - D,
    inverse = 1 / D
  )
}


#' Average Taxonomic Distinctness (Delta+)
#'
#' Calculates the average taxonomic distinctness (AvTD, Delta+) based
#' on Clarke & Warwick (1998). This is a presence/absence-based measure
#' of the average taxonomic distance between all pairs of species.
#'
#' @param species Character vector of species names present in the
#'   community (presence-only data).
#' @param tax_tree A data frame representing the taxonomic hierarchy.
#' @param weights Optional numeric vector of weights for taxonomic
#'   levels.
#'
#' @return A numeric value representing the average taxonomic
#'   distinctness (Delta+).
#'
#' @references
#' Clarke, K.R. & Warwick, R.M. (1998). A taxonomic distinctness index
#' and its statistical properties. Journal of Applied Ecology, 35,
#' 523-531.
#'
#' @examples
#' tax <- data.frame(
#'   Species = c("Quercus_robur", "Pinus_nigra", "Fagus_orientalis",
#'               "Abies_nordmanniana"),
#'   Genus = c("Quercus", "Pinus", "Fagus", "Abies"),
#'   Family = c("Fagaceae", "Pinaceae", "Fagaceae", "Pinaceae"),
#'   Order = c("Fagales", "Pinales", "Fagales", "Pinales"),
#'   stringsAsFactors = FALSE
#' )
#'
#' spp <- c("Quercus_robur", "Pinus_nigra", "Fagus_orientalis")
#' avtd(spp, tax)
#'
#' @export
avtd <- function(species, tax_tree, weights = NULL) {
  if (length(species) < 2) {
    stop("At least 2 species are required to compute AvTD.",
         call. = FALSE)
  }

  dist_mat <- tax_distance_matrix(tax_tree, species, weights)

  # Average of upper triangle
  n <- length(species)
  sum_dist <- sum(dist_mat[upper.tri(dist_mat)])
  n_pairs <- n * (n - 1) / 2

  sum_dist / n_pairs
}


#' Variation in Taxonomic Distinctness (Lambda+)
#'
#' Calculates the variation in taxonomic distinctness (VarTD, Lambda+)
#' based on Clarke & Warwick (2001).
#'
#' @param species Character vector of species names present in the
#'   community.
#' @param tax_tree A data frame representing the taxonomic hierarchy.
#' @param weights Optional numeric vector of weights for taxonomic
#'   levels.
#'
#' @return A numeric value representing the variation in taxonomic
#'   distinctness (Lambda+).
#'
#' @references
#' Clarke, K.R. & Warwick, R.M. (2001). A further biodiversity index
#' applicable to species lists: variation in taxonomic distinctness.
#' Marine Ecology Progress Series, 216, 265-278.
#'
#' @examples
#' tax <- data.frame(
#'   Species = c("Quercus_robur", "Pinus_nigra", "Fagus_orientalis",
#'               "Abies_nordmanniana"),
#'   Genus = c("Quercus", "Pinus", "Fagus", "Abies"),
#'   Family = c("Fagaceae", "Pinaceae", "Fagaceae", "Pinaceae"),
#'   Order = c("Fagales", "Pinales", "Fagales", "Pinales"),
#'   stringsAsFactors = FALSE
#' )
#'
#' spp <- c("Quercus_robur", "Pinus_nigra", "Fagus_orientalis")
#' vartd(spp, tax)
#'
#' @export
vartd <- function(species, tax_tree, weights = NULL) {
  if (length(species) < 2) {
    stop("At least 2 species are required to compute VarTD.",
         call. = FALSE)
  }

  dist_mat <- tax_distance_matrix(tax_tree, species, weights)

  n <- length(species)
  delta_plus <- avtd(species, tax_tree, weights)

  # VarTD = variance of omega_ij around Delta+
  upper_vals <- dist_mat[upper.tri(dist_mat)]
  n_pairs <- n * (n - 1) / 2

  sum((upper_vals - delta_plus)^2) / n_pairs
}
