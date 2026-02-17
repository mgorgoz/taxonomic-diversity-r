#' Calculate Deng Entropy at a Single Taxonomic Level
#'
#' Computes the Deng entropy (Ed) for a given set of group proportions
#' at a specific taxonomic level. This is the core entropy calculation
#' from Deng (2016), which generalizes Shannon entropy through the
#' Dempster-Shafer evidence theory framework.
#'
#' @param abundances A numeric vector of abundances for each group
#'   (node) at the given taxonomic level.
#' @param group_sizes Optional integer vector of focal element sizes
#'   (`|Fi|`) for each group. At species level this is NULL (all sizes
#'   are 1, reducing to Shannon entropy). At higher taxonomic levels,
#'   each value represents the number of species within that group.
#'
#' @return A numeric value representing the Deng entropy at that level.
#'
#' @details
#' The Deng entropy is calculated as:
#' \deqn{E_d = -\sum_{i} m(F_i) \ln \frac{m(F_i)}{2^{|F_i|} - 1}}
#'
#' At species level, each focal element has cardinality 1, so Deng
#' entropy reduces to Shannon entropy:
#' \deqn{E_d^S = H = -\sum_i p_i \ln p_i}
#'
#' At higher levels (genus, family, etc.), \eqn{|F_i|} equals the
#' number of species within each group, and the mass function is
#' the normalized proportion of total abundance in each group.
#'
#' @references
#' Deng, Y. (2016). Deng entropy. Chaos, Solitons & Fractals, 91, 549-553.
#'
#' @seealso [ozkan_pto()] which uses this function internally,
#'   [shannon()] for classical Shannon entropy.
#'
#' @examples
#' # Shannon entropy (species level, |Fi| = 1 for all)
#' deng_entropy_level(c(4, 2, 3, 1, 2, 3, 2, 2))
#'
#' # Deng entropy at genus level with group sizes
#' deng_entropy_level(c(9, 3, 7), group_sizes = c(3, 2, 3))
#'
#' @export
deng_entropy_level <- function(abundances, group_sizes = NULL) {
  if (!is.numeric(abundances) || any(abundances < 0)) {
    stop("'abundances' must be a numeric vector with non-negative values.",
         call. = FALSE)
  }

  # Remove zeros
  if (!is.null(group_sizes)) {
    keep <- abundances > 0
    abundances <- abundances[keep]
    group_sizes <- group_sizes[keep]
  } else {
    abundances <- abundances[abundances > 0]
  }

  if (length(abundances) == 0) return(0)
  if (length(abundances) == 1) return(0)

  # Mass function: proportions
  m <- abundances / sum(abundances)

  if (is.null(group_sizes)) {
    # Species level: |Fi| = 1 for all, reduces to Shannon
    # Ed = -sum(m * ln(m / (2^1 - 1))) = -sum(m * ln(m))
    return(-sum(m * log(m)))
  }

  # Higher levels: |Fi| = group_sizes[i]
  # Ed = -sum(m(Fi) * ln(m(Fi) / (2^|Fi| - 1)))
  ed <- -sum(m * log(m / (2^group_sizes - 1)))
  return(ed)
}
