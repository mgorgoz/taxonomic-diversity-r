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
#' @param correction Bias correction method for Shannon entropy
#'   estimation. Only applied at species level (`group_sizes = NULL`).
#'   One of `"none"` (default), `"miller_madow"`, `"grassberger"`,
#'   or `"chao_shen"`. See [shannon()] for details. A warning is
#'   issued if correction is requested with non-NULL `group_sizes`.
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
#' Bias correction is only meaningful at the species level where Deng
#' entropy equals Shannon entropy. At higher taxonomic levels the mass
#' function has a different structure and bias-correction formulas do
#' not apply.
#'
#' @references
#' Deng, Y. (2016). Deng entropy. Chaos, Solitons & Fractals, 91, 549-553.
#'
#' @seealso [ozkan_pto()] which uses this function internally,
#'   [shannon()] for classical Shannon entropy and bias corrections.
#'
#' @examples
#' # Shannon entropy (species level, |Fi| = 1 for all)
#' deng_entropy_level(c(4, 2, 3, 1, 2, 3, 2, 2))
#'
#' # With bias correction at species level
#' deng_entropy_level(c(4, 2, 3, 1, 2), correction = "chao_shen")
#'
#' # Deng entropy at genus level with group sizes
#' deng_entropy_level(c(9, 3, 7), group_sizes = c(3, 2, 3))
#'
#' @export
deng_entropy_level <- function(abundances, group_sizes = NULL,
                                correction = c("none", "miller_madow",
                                               "grassberger",
                                               "chao_shen")) {
  correction <- match.arg(correction)

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

  if (is.null(group_sizes)) {
    # Species level: |Fi| = 1 for all, reduces to Shannon
    # Delegate to shannon() which handles bias correction
    return(shannon(abundances, base = exp(1), correction = correction))
  }

  # Higher levels: bias correction not applicable
  if (correction != "none") {
    warning("Bias correction is only supported at species level ",
            "(group_sizes = NULL). Falling back to correction = 'none'.",
            call. = FALSE)
  }

  # Mass function: proportions
  m <- abundances / sum(abundances)

  # Higher levels: |Fi| = group_sizes[i]
  # Ed = -sum(m(Fi) * ln(m(Fi) / (2^|Fi| - 1)))
  ed <- -sum(m * log(m / (2^group_sizes - 1)))
  return(ed)
}
