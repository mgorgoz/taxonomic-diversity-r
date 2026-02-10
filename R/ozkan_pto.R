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
#' @param weighted Logical. If TRUE (default), use weighted calculation
#'   where w_i increases from 1 (species) to max level. If FALSE, all
#'   w_i = 1 (unweighted).
#'
#' @return A named list with components:
#'   \describe{
#'     \item{uTO}{Unweighted taxonomic diversity}
#'     \item{TO}{Weighted taxonomic diversity}
#'     \item{uTO_plus}{Unweighted taxonomic distance}
#'     \item{TO_plus}{Weighted taxonomic distance}
#'     \item{Ed_levels}{Deng entropy at each taxonomic level}
#'   }
#'
#' @details
#' The method applies the slicing procedure to abundance data and
#' computes Deng entropy at each taxonomic level of the Linnean
#' hierarchy. The formula is:
#'
#' \deqn{pT_O = \ln \left( \frac{\sum_{k=0}^{n_s} (n_s - n_k)
#'   \prod_{i=1}^{L} \left( w_i \left( \frac{(e^{E_d^S})^2}
#'   {e^{E_d^i}} + 1 \right) \right)}{n_s + \sum n_k} \right)}
#'
#' where \eqn{E_d^S} is the Deng entropy at species level (Shannon H),
#' \eqn{E_d^i} is the Deng entropy at level i, and the product runs
#' over levels where the number of nodes > 1.
#'
#' The taxonomic distance (pTO+) uses only presence/absence data
#' (equivalent to n_s = 1, n_k = 0):
#'
#' \deqn{pT_O^+ = \prod_{i=1}^{L} \left( w_i \left(
#'   \frac{(e^{E_d^S})^2}{e^{E_d^i}} + 1 \right) \right)}
#'
#' @references
#' Ozkan, K. (2018). A new proposed measure for estimating taxonomic
#' diversity. Turkish Journal of Forestry, 19(4), 336-346.
#' DOI: 10.18182/tjf.441061
#'
#' Deng, Y. (2016). Deng entropy. Chaos, Solitons & Fractals, 91,
#' 549-553.
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

  # --- Step 1: Compute Deng entropy at each taxonomic level ---
  Ed <- numeric(n_levels + 1)
  level_names <- c("Species", names(tax_tree)[-1])
  names(Ed) <- level_names

  # Species level: Ed_S = Shannon entropy (|Fi| = 1)
  Ed[1] <- deng_entropy_level(as.numeric(community))

  # Higher levels: aggregate abundances by group, compute Deng entropy
  active_levels <- 1  # Track which levels have > 1 node

  for (lev in seq_len(n_levels)) {
    col_idx <- lev + 1  # Column in tax_tree
    groups <- as.character(tax_sub[[col_idx]])
    unique_groups <- unique(groups)
    n_nodes <- length(unique_groups)

    # Termination rule: if only 1 node, Ed = 0, stop here
    if (n_nodes <= 1) {
      Ed[lev + 1] <- 0
      break
    }

    active_levels <- c(active_levels, lev + 1)

    # Aggregate abundance by group at this level
    group_abundances <- numeric(n_nodes)
    group_sizes <- integer(n_nodes)  # Number of species in each group

    for (g in seq_along(unique_groups)) {
      mask <- groups == unique_groups[g]
      group_abundances[g] <- sum(as.numeric(community[mask]))
      group_sizes[g] <- sum(mask)
    }

    # Deng entropy at this level
    Ed[lev + 1] <- deng_entropy_level(group_abundances,
                                       group_sizes = group_sizes)
  }

  # --- Step 2: Compute the core product (pTO+) ---
  # For both weighted and unweighted versions

  Ed_S <- Ed[1]  # Species-level entropy
  e_Ed_S_sq <- (exp(Ed_S))^2

  # Core product for unweighted (all wi = 1)
  core_unweighted <- 1
  # Core product for weighted (wi = i)
  core_weighted <- 1

  for (k in seq_along(active_levels)) {
    lev_idx <- active_levels[k]
    if (lev_idx == 1) next  # Skip species level in product

    Ed_i <- Ed[lev_idx]
    ratio <- e_Ed_S_sq / exp(Ed_i) + 1

    # Unweighted: wi = 1
    core_unweighted <- core_unweighted * (1 * ratio)

    # Weighted: wi = rank index (species=1, genus=2, ..., kingdom=7)
    w_i <- lev_idx
    core_weighted <- core_weighted * (w_i * ratio)
  }

  # pTO+ (taxonomic distance, presence/absence based)
  uTO_plus <- core_unweighted
  TO_plus <- core_weighted

  # --- Step 3: Slicing procedure for pTO ---
  # Determine slicing steps
  max_abundance <- max(community)
  abundances <- as.numeric(community)

  # Build slicing steps: for each step k, subtract k from abundances
  # Record which species remain (abundance - k > 0)
  n_s <- max_abundance  # Total number of steps

  # Compute sum of (ns - nk) * core_product for each step
  # At each step, we recompute Deng entropy with the reduced community

  sum_weighted <- 0
  sum_unweighted <- 0

  for (step in seq_len(n_s)) {
    nk <- step - 1  # Amount subtracted (0-indexed: step 1 = nk 0)
    factor <- n_s - nk

    # Reduce abundances: subtract nk, keep positive
    reduced <- abundances - nk
    reduced[reduced < 0] <- 0
    # Convert to presence/absence at this step
    present <- as.integer(reduced > 0)

    if (sum(present) == 0) break

    # Recompute Deng entropy at each level with reduced community
    # Only include species that are still present
    present_mask <- present > 0
    comm_step <- reduced[present_mask]
    tax_step <- tax_sub[present_mask, , drop = FALSE]

    # Species-level entropy for this step
    if (length(comm_step) <= 1) {
      Ed_S_step <- 0
    } else {
      Ed_S_step <- deng_entropy_level(comm_step)
    }

    e_Ed_S_sq_step <- (exp(Ed_S_step))^2

    # Core product for this step
    core_u_step <- 1
    core_w_step <- 1

    for (lev in seq_len(n_levels)) {
      col_idx <- lev + 1
      groups_step <- as.character(tax_step[[col_idx]])
      unique_groups_step <- unique(groups_step)
      n_nodes_step <- length(unique_groups_step)

      if (n_nodes_step <= 1) break

      # Aggregate
      ga <- numeric(n_nodes_step)
      gs <- integer(n_nodes_step)
      for (g in seq_along(unique_groups_step)) {
        mask <- groups_step == unique_groups_step[g]
        ga[g] <- sum(comm_step[mask])
        gs[g] <- sum(mask)
      }

      Ed_i_step <- deng_entropy_level(ga, group_sizes = gs)
      ratio_step <- e_Ed_S_sq_step / exp(Ed_i_step) + 1

      core_u_step <- core_u_step * (1 * ratio_step)
      w_i <- lev + 1
      core_w_step <- core_w_step * (w_i * ratio_step)
    }

    sum_unweighted <- sum_unweighted + factor * core_u_step
    sum_weighted <- sum_weighted + factor * core_w_step
  }

  # Denominator: ns + sum(nk) = ns + sum(0:(ns-1)) = ns + ns*(ns-1)/2
  denom <- n_s + sum(0:(n_s - 1))

  # pTO = ln(numerator / denominator)
  uTO <- log(sum_unweighted / denom)
  TO <- log(sum_weighted / denom)

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
