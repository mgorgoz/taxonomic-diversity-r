#' Calculate Ozkan's Taxonomic Diversity Index (pTO)
#'
#' Computes the four components of the Deng entropy-based taxonomic
#' diversity measure proposed by Ozkan (2018): weighted/unweighted
#' taxonomic diversity (TO, uTO) and weighted/unweighted taxonomic
#' distance (TO+, uTO+).
#'
#' @param community A named numeric vector of species abundances.
#'   Names must match the first column of `tax_tree`.
#' @param tax_tree A data frame with taxonomic hierarchy. First column
#'   is species names, subsequent columns are taxonomic ranks from
#'   lowest to highest (e.g., Species, Genus, Family, Order, Class,
#'   Phylum, Kingdom).
#' @param max_level Integer or \code{NULL}. Maximum number of taxonomic
#'   levels (above Species) to include in the product formula. When
#'   \code{NULL} (default), all available levels are used. When set to
#'   \code{"auto"}, the function automatically detects the highest
#'   informative level (where Deng entropy > 0 at nk=0) and truncates
#'   the product there. A positive integer limits to that many levels
#'   (e.g., \code{max_level = 2} uses only Genus and Family).
#'
#' @return A named list with components:
#'   \describe{
#'     \item{uTO}{Unweighted taxonomic diversity (all levels)}
#'     \item{TO}{Weighted taxonomic diversity (all levels)}
#'     \item{uTO_plus}{Unweighted taxonomic distance (all levels)}
#'     \item{TO_plus}{Weighted taxonomic distance (all levels)}
#'     \item{uTO_max}{Unweighted taxonomic diversity (max informative level)}
#'     \item{TO_max}{Weighted taxonomic diversity (max informative level)}
#'     \item{uTO_plus_max}{Unweighted taxonomic distance (max informative level)}
#'     \item{TO_plus_max}{Weighted taxonomic distance (max informative level)}
#'     \item{Ed_levels}{Deng entropy at each taxonomic level (nk=0 slice)}
#'     \item{max_informative_level}{Integer: highest level with Ed > 0}
#'   }
#'
#' @details
#' The method uses the slicing procedure from Ozkan (2018). At each
#' slice (nk = 0, 1, ..., n_s), species with abundance <= nk are removed.
#' The surviving species receive EQUAL weight (1/count) — abundance
#' information enters indirectly through which species survive each slice.
#'
#' Deng entropy at each taxonomic level is computed using these equal
#' proportions, where the mass function m(Fi) = count_in_group / total_count
#' and |Fi| = number of species in that taxonomic group.
#'
#' The core product formula at each slice is:
#'
#' \deqn{\prod_{i=1}^{L} \left( w_i \left( \frac{(e^{E_d^S})^2}
#'   {e^{E_d^i}} + 1 \right) \right)}
#'
#' where \eqn{E_d^S} is the Deng entropy at species level and
#' \eqn{E_d^i} is the Deng entropy at level i, computed using
#' presence/absence (equal weight) proportions.
#'
#' pTO+ (taxonomic distance) uses only the nk=0 slice:
#' \deqn{pT_O^+ = \ln \prod_{i=1}^{L} \left( w_i \left(
#'   \frac{(e^{E_d^S})^2}{e^{E_d^i}} + 1 \right) \right)}
#'
#' pTO (taxonomic diversity) aggregates across all slices:
#' \deqn{pT_O = \ln \left( \frac{\sum_{k=0}^{n_s} (n_s - n_k)
#'   \prod_{i=1}^{L} \left( w_i \left( \frac{(e^{E_d^S})^2}
#'   {e^{E_d^i}} + 1 \right) \right)}{n_s + \sum n_k} \right)}
#'
#' @references
#' Ozkan, K. (2018). A new proposed measure for estimating taxonomic
#' diversity. Turkish Journal of Forestry, 19(4), 336-346.
#' DOI: 10.18182/tjf.441061
#'
#' Deng, Y. (2016). Deng entropy. Chaos, Solitons & Fractals, 91,
#' 549-553.
#'
#' @seealso [deng_entropy_level()] for the core Deng entropy calculation,
#'   [pto_components()] for a convenience wrapper returning a named vector,
#'   [delta()] and [avtd()] for Clarke & Warwick alternatives.
#'
#' @examples
#' # Simple example with 5 species
#' comm <- c(sp1 = 4, sp2 = 2, sp3 = 3, sp4 = 1, sp5 = 2)
#' tax <- data.frame(
#'   Species = paste0("sp", 1:5),
#'   Genus   = c("G1", "G1", "G1", "G2", "G2"),
#'   Family  = c("F1", "F1", "F1", "F1", "F1"),
#'   stringsAsFactors = FALSE
#' )
#' ozkan_pto(comm, tax)
#'
#' # With auto max-level detection
#' ozkan_pto(comm, tax, max_level = "auto")
#'
#' @export
ozkan_pto <- function(community, tax_tree, max_level = NULL) {
  # --- Input validation ---
  if (!is.numeric(community) || any(community < 0)) {
    stop("'community' must be a numeric vector with non-negative values.",
         call. = FALSE)
  }
  if (is.null(names(community))) {
    stop("'community' must be a named vector.", call. = FALSE)
  }
  if (!is.data.frame(tax_tree)) {
    stop("'tax_tree' must be a data frame.", call. = FALSE)
  }
  if (ncol(tax_tree) < 2) {
    stop("'tax_tree' must have at least 2 columns.", call. = FALSE)
  }

  # Filter to non-zero abundances
  community <- community[community > 0]
  species_names <- names(community)

  # Check species in taxonomy
  tax_species <- as.character(tax_tree[[1]])
  missing_sp <- setdiff(species_names, tax_species)
  if (length(missing_sp) > 0) {
    stop("Species not found in tax_tree: ",
         paste(missing_sp, collapse = ", "), call. = FALSE)
  }

  # Subset taxonomy to community species
  idx <- match(species_names, tax_species)
  tax_sub <- tax_tree[idx, , drop = FALSE]

  n_species <- length(community)
  n_levels <- ncol(tax_tree) - 1  # Excluding species column

  # --- Validate max_level ---
  use_auto_max <- FALSE
  if (!is.null(max_level)) {
    if (is.character(max_level) && tolower(max_level) == "auto") {
      use_auto_max <- TRUE
    } else if (is.numeric(max_level)) {
      max_level <- as.integer(max_level)
      if (max_level < 1 || max_level > n_levels) {
        stop("'max_level' must be between 1 and ", n_levels,
             " (number of taxonomic ranks above Species).", call. = FALSE)
      }
    } else {
      stop("'max_level' must be NULL, \"auto\", or a positive integer.",
           call. = FALSE)
    }
  }

  # --- If only one species, pTO = 0 ---
  if (n_species == 1) {
    result <- list(
      uTO = 0, TO = 0, uTO_plus = 0, TO_plus = 0,
      uTO_max = 0, TO_max = 0, uTO_plus_max = 0, TO_plus_max = 0,
      Ed_levels = setNames(rep(0, n_levels + 1),
                           c("Species", names(tax_tree)[-1])),
      max_informative_level = 0L
    )
    class(result) <- "ozkan_pto"
    return(result)
  }

  # --- Helper: compute Deng entropy at all levels for a given set of
  #     present species (using EQUAL weights / presence-absence) ---
  compute_deng_all_levels <- function(present_mask, tax_data, n_lev) {
    tax_present <- tax_data[present_mask, , drop = FALSE]
    n_present <- sum(present_mask)

    if (n_present <= 1) {
      ed <- setNames(rep(0, n_lev + 1),
                     c("Species", names(tax_data)[-1]))
      return(list(Ed = ed, active = integer(0)))
    }

    ed <- numeric(n_lev + 1)
    names(ed) <- c("Species", names(tax_data)[-1])

    # Species level: equal weights -> Ed_S = ln(S)
    # Each species gets m_i = 1/S, |Fi| = 1
    # Ed = -sum((1/S) * ln((1/S) / (2^1 - 1))) = -sum((1/S) * ln(1/S)) = ln(S)
    ed[1] <- log(n_present)

    active <- 1L  # Species level is always active if > 1 species

    for (lev in seq_len(n_lev)) {
      col_idx <- lev + 1
      groups <- as.character(tax_present[[col_idx]])
      unique_groups <- unique(groups)
      n_nodes <- length(unique_groups)

      if (n_nodes <= 1) {
        ed[lev + 1] <- 0
        break
      }

      active <- c(active, as.integer(lev + 1))

      # Count species per group (presence-based)
      group_counts <- numeric(n_nodes)
      group_sizes <- integer(n_nodes)  # |Fi| = species in each group

      for (g in seq_along(unique_groups)) {
        mask <- groups == unique_groups[g]
        group_counts[g] <- sum(mask)   # Number of present species
        group_sizes[g] <- sum(mask)    # Same as count for presence-based
      }

      # Deng entropy: m_i = count_i / total, |Fi| = group_sizes[i]
      ed[lev + 1] <- deng_entropy_level(group_counts,
                                         group_sizes = group_sizes)
    }

    return(list(Ed = ed, active = active))
  }

  # --- Step 1: Compute Deng entropy at nk=0 (all species present) ---
  all_present <- rep(TRUE, n_species)
  result_nk0 <- compute_deng_all_levels(all_present, tax_sub, n_levels)
  Ed <- result_nk0$Ed
  active_levels <- result_nk0$active

  # --- Determine max informative level (highest level with Ed > 0) ---
  # active_levels contains 1-based indices of levels with Ed > 0
  # Species level (index 1) is always active; we want the highest rank level
  max_info_level <- if (length(active_levels) > 0) max(active_levels) else 1L

  # Determine which levels to use for "max" calculation
  if (use_auto_max) {
    max_active <- active_levels  # only levels with Ed > 0
  } else if (is.integer(max_level)) {
    # max_level = number of rank levels above species to include
    # active_levels indices: 1=Species, 2=Genus, 3=Family, ...
    # So we keep active_levels where index <= max_level + 1
    max_active <- active_levels[active_levels <= (max_level + 1L)]
  } else {
    # NULL: max version uses auto-detected informative levels
    max_active <- active_levels
  }

  # --- Helper: compute core product for a set of active levels ---
  compute_product <- function(Ed_vec, levels_to_use) {
    Ed_S_local <- Ed_vec[1]
    e_Ed_S_sq_local <- (exp(Ed_S_local))^2
    core_u <- 1
    core_w <- 1
    for (lv in seq_along(levels_to_use)) {
      lev_idx <- levels_to_use[lv]
      Ed_i <- Ed_vec[lev_idx]
      ratio <- e_Ed_S_sq_local / exp(Ed_i) + 1
      core_u <- core_u * (1 * ratio)
      w_i <- lev_idx
      core_w <- core_w * (w_i * ratio)
    }
    list(unweighted = core_u, weighted = core_w)
  }

  # --- Step 2: Compute the core product (pTO+) using nk=0 entropies ---
  # Full version (all active levels)
  prod_full <- compute_product(Ed, active_levels)
  uTO_plus <- unname(log(prod_full$unweighted))
  TO_plus <- unname(log(prod_full$weighted))

  # Max version (informative levels only)
  prod_max <- compute_product(Ed, max_active)
  uTO_plus_max <- unname(log(prod_max$unweighted))
  TO_plus_max <- unname(log(prod_max$weighted))

  # --- Step 3: Slicing procedure for pTO ---
  max_abundance <- max(community)
  abundances <- as.numeric(community)
  n_s <- max_abundance

  # Accumulators for full (all levels) and max (informative levels)
  sum_u_full <- 0
  sum_w_full <- 0
  sum_u_max <- 0
  sum_w_max <- 0

  for (step in seq_len(n_s)) {
    nk <- step - 1  # Amount subtracted (0-indexed: step 1 = nk 0)
    factor_val <- n_s - nk

    # Species survive if abundance > nk
    present_mask <- abundances > nk

    if (sum(present_mask) == 0) break

    # Compute Deng entropy for this slice using presence-based weights
    tax_step <- tax_sub[present_mask, , drop = FALSE]
    n_present <- sum(present_mask)

    if (n_present <= 1) {
      # Single species: Ed = 0 at all levels, product = 1+1 = 2
      ratio_species <- 0 + 1  # (e^0)^2 / e^0 + 1 = 1/1 + 1 = 2
      core_u_full <- 1 * ratio_species
      core_w_full <- 1 * ratio_species
      core_u_max <- core_u_full
      core_w_max <- core_w_full
    } else {
      # Compute full Deng entropy at all levels for this slice
      result_step <- compute_deng_all_levels(
        rep(TRUE, n_present), tax_step, n_levels)
      Ed_step <- result_step$Ed
      active_step <- result_step$active

      # Determine max-active levels for this slice
      if (use_auto_max) {
        max_active_step <- active_step
      } else if (is.integer(max_level)) {
        max_active_step <- active_step[active_step <= (max_level + 1L)]
      } else {
        max_active_step <- active_step
      }

      # Full product (all active levels)
      prod_step_full <- compute_product(Ed_step, active_step)
      core_u_full <- prod_step_full$unweighted
      core_w_full <- prod_step_full$weighted

      # Max product (informative levels only)
      prod_step_max <- compute_product(Ed_step, max_active_step)
      core_u_max <- prod_step_max$unweighted
      core_w_max <- prod_step_max$weighted
    }

    sum_u_full <- sum_u_full + factor_val * core_u_full
    sum_w_full <- sum_w_full + factor_val * core_w_full
    sum_u_max <- sum_u_max + factor_val * core_u_max
    sum_w_max <- sum_w_max + factor_val * core_w_max
  }

  # Denominator: ns + sum(nk) = ns + sum(0:(ns-1)) = ns + ns*(ns-1)/2
  denom <- n_s + sum(0:(n_s - 1))

  # pTO = ln(numerator / denominator)
  uTO <- unname(log(sum_u_full / denom))
  TO <- unname(log(sum_w_full / denom))
  uTO_max <- unname(log(sum_u_max / denom))
  TO_max <- unname(log(sum_w_max / denom))

  result <- list(
    uTO = uTO,
    TO = TO,
    uTO_plus = uTO_plus,
    TO_plus = TO_plus,
    uTO_max = uTO_max,
    TO_max = TO_max,
    uTO_plus_max = uTO_plus_max,
    TO_plus_max = TO_plus_max,
    Ed_levels = Ed,
    max_informative_level = as.integer(max_info_level - 1L)
  )
  class(result) <- "ozkan_pto"
  return(result)
}


#' Calculate All Eight pTO Components (Convenience Wrapper)
#'
#' Returns a named numeric vector with all eight Ozkan (2018) components:
#' four using all taxonomic levels and four using only the informative
#' levels (max version), matching the Excel macro's Run 1+2+3 output.
#'
#' @inheritParams ozkan_pto
#'
#' @return A named numeric vector with eight elements:
#'   \code{uTO}, \code{TO}, \code{uTO_plus}, \code{TO_plus},
#'   \code{uTO_max}, \code{TO_max}, \code{uTO_plus_max}, \code{TO_plus_max}.
#'
#' @seealso [ozkan_pto()] for the full result including per-level entropy.
#'
#' @examples
#' comm <- c(sp1 = 4, sp2 = 2, sp3 = 3, sp4 = 1, sp5 = 2)
#' tax <- data.frame(
#'   Species = paste0("sp", 1:5),
#'   Genus   = c("G1", "G1", "G1", "G2", "G2"),
#'   Family  = c("F1", "F1", "F1", "F1", "F1"),
#'   stringsAsFactors = FALSE
#' )
#' pto_components(comm, tax)
#'
#' @export
pto_components <- function(community, tax_tree) {
  result <- ozkan_pto(community, tax_tree, max_level = "auto")
  c(uTO          = unname(result$uTO),
    TO           = unname(result$TO),
    uTO_plus     = unname(result$uTO_plus),
    TO_plus      = unname(result$TO_plus),
    uTO_max      = unname(result$uTO_max),
    TO_max       = unname(result$TO_max),
    uTO_plus_max = unname(result$uTO_plus_max),
    TO_plus_max  = unname(result$TO_plus_max))
}
