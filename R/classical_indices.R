#' Shannon Diversity Index
#'
#' Calculates the Shannon-Wiener diversity index (H') for a community,
#' optionally applying a bias correction for small samples.
#'
#' @param community A numeric vector of species abundances (counts).
#' @param base The logarithm base. Default is `exp(1)` (natural log).
#'   Use `2` for bits.
#' @param correction Bias correction method. One of `"none"` (default,
#'   naive MLE), `"miller_madow"`, `"grassberger"`, or `"chao_shen"`.
#'   See Details.
#'
#' @return A numeric value representing the Shannon diversity index.
#'
#' @details
#' The naive (MLE) Shannon index is calculated as:
#' \deqn{H' = -\sum_{i=1}^{S} p_i \ln(p_i)}
#' where \eqn{p_i = n_i / N} is the proportion of species \eqn{i},
#' \eqn{N} is the total number of individuals, and \eqn{S} is the
#' number of species observed.
#'
#' The MLE estimator has a known negative bias that is significant for
#' small samples. Three bias-correction methods are available:
#'
#' **Miller-Madow** (1955): Adds a first-order bias correction term:
#' \deqn{H_{MM} = H_{MLE} + \frac{S_{obs} - 1}{2N}}
#'
#' **Grassberger** (2003): Uses the digamma function instead of the
#' logarithm:
#' \deqn{H_G = \ln(N) - \frac{1}{N} \sum_i n_i \psi(n_i)}
#' where \eqn{\psi} is the digamma function.
#'
#' **Chao-Shen** (2003): Applies a Good-Turing coverage correction
#' with Horvitz-Thompson weighting:
#' \deqn{\hat{C} = 1 - f_1 / N}
#' \deqn{H_{CS} = -\sum_i \frac{\hat{p}_i \ln \hat{p}_i}{1 -
#'   (1 - \hat{p}_i)^N}}
#' where \eqn{\hat{p}_i = \hat{C} \cdot n_i / N} and \eqn{f_1} is
#' the number of singletons.
#'
#' Bias corrections require integer abundance counts. A warning is
#' issued if non-integer values are detected with `correction != "none"`.
#'
#' @references
#' Miller, G.A. & Madow, W.G. (1954). On the maximum likelihood
#' estimate of the Shannon-Wiener index of diversity. AFCRC-TR-54-75.
#'
#' Grassberger, P. (2003). Entropy estimates from insufficient
#' samplings. arXiv:physics/0307138.
#'
#' Chao, A. & Shen, T.-J. (2003). Nonparametric estimation of
#' Shannon's index of diversity when there are unseen species in
#' sample. Environmental and Ecological Statistics, 10, 429-443.
#'
#' @seealso [simpson()] for Simpson diversity, [deng_entropy_level()] for
#'   Deng entropy (a generalization of Shannon).
#'
#' @examples
#' comm <- c(10, 5, 8, 3, 12)
#' shannon(comm)
#' shannon(comm, correction = "miller_madow")
#' shannon(comm, correction = "grassberger")
#' shannon(comm, correction = "chao_shen")
#'
#' @export
shannon <- function(community, base = exp(1),
                    correction = c("none", "miller_madow",
                                   "grassberger", "chao_shen")) {
  correction <- match.arg(correction)

  if (!is.numeric(community) || any(community < 0)) {
    stop("'community' must be a numeric vector with non-negative values.",
         call. = FALSE)
  }

  # Remove zeros
  community <- community[community > 0]

  if (length(community) == 0) return(0)
  if (length(community) == 1) return(0)

  # Warn if non-integer counts are used with correction
  if (correction != "none" && any(abs(community - round(community)) > 1e-9)) {
    warning("Bias correction requires integer abundance counts. ",
            "Non-integer values detected.", call. = FALSE)
  }

  N <- sum(community)
  S_obs <- length(community)

  if (correction == "none") {
    p <- community / N
    return(-sum(p * log(p, base = base)))
  }

  if (correction == "miller_madow") {
    p <- community / N
    H_naive <- -sum(p * log(p, base = base))
    # Bias term is in nats; convert to requested base
    return(H_naive + (S_obs - 1) / (2 * N * log(base)))
  }

  if (correction == "grassberger") {
    # H_G = log(N) - (1/N) * sum(n_i * digamma(n_i))
    n <- community
    H_G <- log(N) - (1 / N) * sum(n * digamma(n))
    # Convert from nats to requested base
    return(H_G / log(base))
  }

  if (correction == "chao_shen") {
    n <- community
    f1 <- sum(n == 1)  # number of singletons
    C_hat <- 1 - f1 / N
    # Guard: if all species are singletons, C_hat = 0
    if (C_hat == 0) C_hat <- 1 / N
    p_hat <- C_hat * n / N
    # Horvitz-Thompson correction: avoid division by zero
    ht_denom <- 1 - (1 - p_hat)^N
    ht_denom[ht_denom == 0] <- 1e-10
    H_CS <- -sum(p_hat * log(p_hat, base = base) / ht_denom)
    return(H_CS)
  }
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
#' @seealso [shannon()] for Shannon diversity.
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


#' Taxonomic Diversity Index (Delta)
#'
#' Calculates the taxonomic diversity index (Delta) from Warwick &
#' Clarke (1995). This is the average weighted path length between
#' every pair of individuals, including same-species pairs (weighted 0).
#'
#' @param community A named numeric vector of species abundances.
#' @param tax_tree A data frame with taxonomic hierarchy.
#' @param weights Optional numeric vector of path weights for each
#'   taxonomic level. If NULL, a linear scale is used (1, 2, 3, ...).
#'
#' @return A numeric value representing taxonomic diversity (Delta).
#'
#' @details
#' \deqn{\Delta = \frac{\sum\sum_{i<j} w_{ij} x_i x_j + \sum_i 0
#'   \cdot x_i(x_i-1)/2}{\sum\sum_{i<j} x_i x_j + \sum_i
#'   x_i(x_i-1)/2}}
#'
#' @references
#' Warwick, R.M. & Clarke, K.R. (1995). New 'biodiversity' measures
#' reveal a decrease in taxonomic distinctness with increasing stress.
#' Marine Ecology Progress Series, 129, 301-305.
#'
#' @seealso [delta_star()] for taxonomic distinctness (excluding same-species),
#'   [avtd()] for presence/absence-based AvTD,
#'   [ozkan_pto()] for Deng entropy-based alternative.
#'
#' @examples
#' comm <- c(sp1 = 5, sp2 = 3, sp3 = 3, sp4 = 1, sp5 = 3)
#' tax <- data.frame(
#'   Species = paste0("sp", 1:5),
#'   Genus = c("G1", "G1", "G2", "G2", "G2"),
#'   Family = c("F1", "F1", "F1", "F2", "F2"),
#'   stringsAsFactors = FALSE
#' )
#' delta(comm, tax)
#'
#' @export
delta <- function(community, tax_tree, weights = NULL) {
  if (!is.numeric(community) || any(community < 0)) {
    stop("'community' must be a numeric vector with non-negative values.",
         call. = FALSE)
  }

  community <- community[community > 0]
  if (length(community) < 2) return(0)

  species_names <- names(community)
  if (is.null(species_names)) {
    stop("'community' must be a named vector.", call. = FALSE)
  }

  n_levels <- ncol(tax_tree) - 1
  if (is.null(weights)) {
    weights <- seq_len(n_levels)
  }
  if (length(weights) != n_levels) {
    stop("'weights' must have length ", n_levels,
         " (one per taxonomic level), got ", length(weights), ".",
         call. = FALSE)
  }

  # Get path weights between all species pairs
  n_sp <- length(community)
  x <- as.numeric(community)

  # Find the taxonomic level at which each pair diverges
  tax_species <- as.character(tax_tree[[1]])
  idx <- match(species_names, tax_species)
  if (any(is.na(idx))) {
    missing <- species_names[is.na(idx)]
    stop("Species not found in tax_tree: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  tax_sub <- tax_tree[idx, , drop = FALSE]

  # Numerator: sum of w_ij * xi * xj for i < j
  numerator <- 0
  # Denominator: sum of xi*xj for i<j + sum of xi*(xi-1)/2
  denom_cross <- 0

  for (i in seq_len(n_sp - 1)) {
    for (j in (i + 1):n_sp) {
      # Clarke & Warwick path length: number of steps from species
      # to their lowest common ancestor. With equal step lengths,
      # w_ij = index of the first MATCHING level from bottom.
      # Genus=1, Family=2, Order=3, ... If same genus -> w=1,
      # same family but different genus -> w=2, etc.
      w_ij <- n_levels  # max distance if nothing matches
      for (lev in seq_len(n_levels)) {
        if (as.character(tax_sub[i, lev + 1]) ==
            as.character(tax_sub[j, lev + 1])) {
          w_ij <- weights[lev]
          break
        }
      }

      numerator <- numerator + w_ij * x[i] * x[j]
      denom_cross <- denom_cross + x[i] * x[j]
    }
  }

  # Same-species pairs contribute 0 to numerator
  denom_same <- sum(x * (x - 1) / 2)

  numerator / (denom_cross + denom_same)
}


#' Taxonomic Distinctness (Delta*)
#'
#' Calculates the taxonomic distinctness (Delta*) from Warwick &
#' Clarke (1995). This is the average weighted path length between
#' individuals of different species only.
#'
#' @inheritParams delta
#'
#' @return A numeric value representing taxonomic distinctness (Delta*).
#'
#' @details
#' \deqn{\Delta^* = \frac{\sum\sum_{i<j} w_{ij} x_i x_j}
#'   {\sum\sum_{i<j} x_i x_j}}
#'
#' @references
#' Warwick, R.M. & Clarke, K.R. (1995). New 'biodiversity' measures
#' reveal a decrease in taxonomic distinctness with increasing stress.
#' Marine Ecology Progress Series, 129, 301-305.
#'
#' @seealso [delta()] for taxonomic diversity (including same-species),
#'   [avtd()] and [vartd()] for presence/absence measures.
#'
#' @examples
#' comm <- c(sp1 = 5, sp2 = 3, sp3 = 3, sp4 = 1, sp5 = 3)
#' tax <- data.frame(
#'   Species = paste0("sp", 1:5),
#'   Genus = c("G1", "G1", "G2", "G2", "G2"),
#'   Family = c("F1", "F1", "F1", "F2", "F2"),
#'   stringsAsFactors = FALSE
#' )
#' delta_star(comm, tax)
#'
#' @export
delta_star <- function(community, tax_tree, weights = NULL) {
  if (!is.numeric(community) || any(community < 0)) {
    stop("'community' must be a numeric vector with non-negative values.",
         call. = FALSE)
  }

  community <- community[community > 0]
  if (length(community) < 2) return(0)

  species_names <- names(community)
  if (is.null(species_names)) {
    stop("'community' must be a named vector.", call. = FALSE)
  }

  n_levels <- ncol(tax_tree) - 1
  if (is.null(weights)) {
    weights <- seq_len(n_levels)
  }
  if (length(weights) != n_levels) {
    stop("'weights' must have length ", n_levels,
         " (one per taxonomic level), got ", length(weights), ".",
         call. = FALSE)
  }

  n_sp <- length(community)
  x <- as.numeric(community)

  tax_species <- as.character(tax_tree[[1]])
  idx <- match(species_names, tax_species)
  if (any(is.na(idx))) {
    missing <- species_names[is.na(idx)]
    stop("Species not found in tax_tree: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  tax_sub <- tax_tree[idx, , drop = FALSE]

  numerator <- 0
  denominator <- 0

  for (i in seq_len(n_sp - 1)) {
    for (j in (i + 1):n_sp) {
      # Clarke & Warwick path length: first matching level from bottom
      w_ij <- n_levels  # max distance if nothing matches
      for (lev in seq_len(n_levels)) {
        if (as.character(tax_sub[i, lev + 1]) ==
            as.character(tax_sub[j, lev + 1])) {
          w_ij <- weights[lev]
          break
        }
      }

      numerator <- numerator + w_ij * x[i] * x[j]
      denominator <- denominator + x[i] * x[j]
    }
  }

  if (denominator == 0) return(0)
  numerator / denominator
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
  upper_vals <- dist_mat[upper.tri(dist_mat)]
  n_pairs <- n * (n - 1) / 2

  # Delta+ = mean of upper triangle (same as avtd)
  delta_plus <- sum(upper_vals) / n_pairs

  # VarTD = variance of omega_ij around Delta+
  sum((upper_vals - delta_plus)^2) / n_pairs
}
