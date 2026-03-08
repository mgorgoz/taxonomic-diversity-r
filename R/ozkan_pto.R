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
#'
#' @return A named list with components:
#'   \describe{
#'     \item{uTO}{Unweighted taxonomic diversity}
#'     \item{TO}{Weighted taxonomic diversity}
#'     \item{uTO_plus}{Unweighted taxonomic distance}
#'     \item{TO_plus}{Weighted taxonomic distance}
#'     \item{Ed_levels}{Deng entropy at each taxonomic level (nk=0 slice)}
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
#' pTO (taxonomic diversity) aggregates across slices (max 9 steps,
#' matching the original Ozkan macro implementation):
#' \deqn{pT_O = \ln \left( \frac{\sum_{k=0}^{n_a-1} (n_a - k)
#'   \prod_{i=1}^{L} \left( w_i \left( \frac{(e^{E_d^S})^2}
#'   {e^{E_d^i}} + 1 \right) \right)}{n_a (n_a + 1) / 2} \right)}
#'
#' where \eqn{n_a} is the number of active slicing steps (steps where
#' at least one species survives, capped at 9).
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
#' @export
ozkan_pto <- function(community, tax_tree) {
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

  # --- If only one species, pTO = 0 ---
  if (n_species == 1) {
    return(list(
      uTO = 0, TO = 0, uTO_plus = 0, TO_plus = 0,
      Ed_levels = setNames(rep(0, n_levels + 1),
                           c("Species", names(tax_tree)[-1]))
    ))
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

  # --- Step 2: Compute the core product (pTO+) using nk=0 entropies ---
  Ed_S <- Ed[1]  # Species-level entropy = ln(S)
  e_Ed_S_sq <- (exp(Ed_S))^2  # S^2

  core_unweighted <- 1
  core_weighted <- 1

  for (lv in seq_along(active_levels)) {
    lev_idx <- active_levels[lv]
    Ed_i <- Ed[lev_idx]
    ratio <- e_Ed_S_sq / exp(Ed_i) + 1

    # Unweighted: wi = 1
    core_unweighted <- core_unweighted * (1 * ratio)

    # Weighted: wi = rank index (species=1, genus=2, ..., kingdom=7)
    w_i <- lev_idx
    core_weighted <- core_weighted * (w_i * ratio)
  }

  # pTO+ (taxonomic distance, presence/absence based)
  uTO_plus <- unname(log(core_unweighted))
  TO_plus <- unname(log(core_weighted))

  # --- Step 3: Slicing procedure for pTO ---
  # Following the original Ozkan macro: max 9 slicing steps (nk = 0..8).
  # Each step removes species with abundance <= nk.
  # Only steps where at least 2 species survive are considered "active"
  # (single-species steps have trivial entropy and are excluded).
  # Weights are based on n_active (count of multi-species steps).
  max_abundance <- max(community)
  abundances <- as.numeric(community)
  max_steps <- 9L  # Macro uses 9 TE sheets (1TE..9TE)
  n_steps <- min(max_abundance, max_steps)

  # Compute products at each slice (only for steps with >= 2 species)
  products_u <- numeric(0)  # unweighted products (multi-species only)
  products_w <- numeric(0)  # weighted products (multi-species only)

  for (step in seq_len(n_steps)) {
    nk <- step - 1  # Threshold (0-indexed: step 1 = nk 0)

    # Species survive if abundance > nk
    present_mask <- abundances > nk
    n_present <- sum(present_mask)

    if (n_present <= 1) next  # Skip single-species or empty steps

    # Compute Deng entropy for this slice using presence-based weights
    tax_step <- tax_sub[present_mask, , drop = FALSE]

    result_step <- compute_deng_all_levels(
      rep(TRUE, n_present), tax_step, n_levels)
    Ed_step <- result_step$Ed
    active_step <- result_step$active

    Ed_S_step <- Ed_step[1]
    e_Ed_S_sq_step <- (exp(Ed_S_step))^2

    core_u_step <- 1
    core_w_step <- 1

    for (a in seq_along(active_step)) {
      lev_idx <- active_step[a]
      Ed_i_step <- Ed_step[lev_idx]
      ratio_step <- e_Ed_S_sq_step / exp(Ed_i_step) + 1

      core_u_step <- core_u_step * (1 * ratio_step)

      w_i <- lev_idx
      core_w_step <- core_w_step * (w_i * ratio_step)
    }

    products_u <- c(products_u, core_u_step)
    products_w <- c(products_w, core_w_step)
  }

  # Count active (multi-species) steps and compute weighted average
  n_active <- length(products_u)

  if (n_active == 0) {
    uTO <- 0
    TO <- 0
  } else {
    # Weights: n_active, n_active-1, ..., 1 (descending for ordered steps)
    weights_vec <- n_active:1
    denom <- sum(weights_vec)  # = n_active * (n_active + 1) / 2

    sum_unweighted <- sum(products_u * weights_vec)
    sum_weighted <- sum(products_w * weights_vec)

    # pTO = ln(weighted_average)
    uTO <- log(sum_unweighted / denom)
    TO <- log(sum_weighted / denom)
  }

  return(list(
    uTO = uTO,
    TO = TO,
    uTO_plus = uTO_plus,
    TO_plus = TO_plus,
    Ed_levels = Ed
  ))
}


#' Calculate All Four pTO Components (Convenience Wrapper)
#'
#' Returns a named numeric vector with all four Ozkan (2018) components.
#'
#' @inheritParams ozkan_pto
#'
#' @return A named numeric vector with uTO, TO, uTO_plus, TO_plus.
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
  result <- ozkan_pto(community, tax_tree)
  c(uTO = unname(result$uTO),
    TO = unname(result$TO),
    uTO_plus = unname(result$uTO_plus),
    TO_plus = unname(result$TO_plus))
}
