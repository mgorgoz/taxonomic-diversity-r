#' Jackknife Analysis for Ozkan's pTO Index (Islem 1 / Run 1)
#'
#' Implements the leave-one-out jackknife procedure from the Ozkan Excel
#' macro (Islem 1). Removes each species one at a time, recalculates pTO,
#' and identifies "happy" (contributing) and "unhappy" (non-contributing)
#' species. A species is "happy" if its removal decreases the pTO index,
#' indicating it positively contributes to the community's taxonomic diversity.
#'
#' @param community A named numeric vector of species abundances.
#'   Names must match the first column of `tax_tree`.
#' @param tax_tree A data frame with taxonomic hierarchy. First column
#'   is species names, subsequent columns are taxonomic ranks.
#' @param component Character string specifying which pTO component to use
#'   for the happy/unhappy classification. One of `"uTO_plus"` (default),
#'   `"TO_plus"`, `"uTO"`, or `"TO"`.
#'
#' @return A named list with components:
#'   \describe{
#'     \item{full_result}{The [ozkan_pto()] result for the full community}
#'     \item{jackknife_results}{Data frame with leave-one-out results
#'       per species}
#'     \item{species_status}{Named logical vector: `TRUE` = happy
#'       (contributing), `FALSE` = unhappy (non-contributing)}
#'     \item{n_happy}{Number of happy species}
#'     \item{n_unhappy}{Number of unhappy species}
#'   }
#'
#' @details
#' The jackknife procedure follows the Excel macro's Islem 1 logic:
#' \enumerate{
#'   \item Compute pTO for the full community.
#'   \item For each species i, remove it and compute pTO for the remaining
#'     community (leave-one-out).
#'   \item Compare each leave-one-out result against the full-community value.
#'   \item If removing species i DECREASES the specified component (pTO becomes
#'     smaller), species i is classified as "happy" (contributing).
#'   \item If removing species i does NOT decrease the component, species i is
#'     classified as "unhappy" (non-contributing).
#' }
#'
#' The happy/unhappy classification is used by [ozkan_pto_resample()]
#' (Islem 2) and [ozkan_pto_sensitivity()] (Islem 3) to apply
#' different resampling probabilities to each species category.
#'
#' @references
#' Ozkan, K. (2018). A new proposed measure for estimating taxonomic
#' diversity. Turkish Journal of Forestry, 19(4), 336-346.
#' DOI: 10.18182/tjf.441061
#'
#' @seealso [ozkan_pto()] for the core calculation,
#'   [ozkan_pto_resample()] for Run 2,
#'   [ozkan_pto_full()] for the full 3-run pipeline.
#'
#' @examples
#' comm <- c(sp1 = 4, sp2 = 2, sp3 = 3, sp4 = 1, sp5 = 2)
#' tax <- data.frame(
#'   Species = paste0("sp", 1:5),
#'   Genus   = c("G1", "G1", "G1", "G2", "G2"),
#'   Family  = c("F1", "F1", "F1", "F1", "F1"),
#'   stringsAsFactors = FALSE
#' )
#' jk <- ozkan_pto_jackknife(comm, tax)
#' jk$species_status   # Which species are happy (contributing)?
#' jk$n_happy           # How many happy species?
#'
#' @export
ozkan_pto_jackknife <- function(community, tax_tree,
                                 component = "uTO_plus") {
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
  valid_components <- c("uTO_plus", "TO_plus", "uTO", "TO")
  if (!component %in% valid_components) {
    stop("'component' must be one of: ",
         paste(valid_components, collapse = ", "), call. = FALSE)
  }

  # Filter to non-zero abundances
  community <- community[community > 0]
  if (length(community) < 2) {
    stop("Need at least 2 species with positive abundance for jackknife.",
         call. = FALSE)
  }

  species_names <- names(community)
  n_species <- length(community)

  # --- Step 1: Compute pTO for full community ---
  full_result <- ozkan_pto(community, tax_tree)
  full_value <- full_result[[component]]

  # --- Step 2: Leave-one-out for each species ---
  loo_results <- matrix(NA_real_, nrow = n_species, ncol = 4)
  colnames(loo_results) <- c("uTO_plus", "TO_plus", "uTO", "TO")

  for (i in seq_len(n_species)) {
    loo_comm <- community[-i]

    if (length(loo_comm) < 2) {
      loo_results[i, ] <- c(0, 0, 0, 0)
    } else {
      loo <- ozkan_pto(loo_comm, tax_tree)
      loo_results[i, ] <- c(unname(loo$uTO_plus), unname(loo$TO_plus),
                             unname(loo$uTO), unname(loo$TO))
    }
  }

  # --- Step 3: Classify species ---
  # Happy  = removing it DECREASES the component (species contributes)
  # Unhappy = removing it does NOT decrease (species doesn't contribute)
  is_happy <- loo_results[, component] < full_value
  names(is_happy) <- species_names

  # Build results data frame
  jk_df <- as.data.frame(loo_results)
  jk_df$species <- species_names
  jk_df$is_happy <- is_happy
  jk_df <- jk_df[, c("species", "uTO_plus", "TO_plus", "uTO", "TO",
                       "is_happy")]

  list(
    full_result       = full_result,
    jackknife_results = jk_df,
    species_status    = is_happy,
    n_happy           = sum(is_happy),
    n_unhappy         = sum(!is_happy)
  )
}


#' Stochastic Resampling of Ozkan's pTO Index (Islem 2 / Run 2)
#'
#' Implements the stochastic resampling procedure from Ozkan's Excel macro
#' (Islem 2). First performs a jackknife (Islem 1) to identify "happy"
#' (contributing) and "unhappy" (non-contributing) species, then runs
#' stochastic resampling where unhappy species are always included and
#' happy species are randomly included with 50\% probability.
#'
#' @param community A named numeric vector of species abundances.
#'   Names must match the first column of `tax_tree`.
#' @param tax_tree A data frame with taxonomic hierarchy. First column
#'   is species names, subsequent columns are taxonomic ranks.
#' @param n_iter Number of stochastic iterations to run (default: 101).
#'   Must be >= 101.
#' @param seed Optional random seed for reproducibility (default: NULL).
#'
#' @return A named list with components:
#'   \describe{
#'     \item{uTO_plus_max}{Maximum unweighted taxonomic distance across
#'       iterations}
#'     \item{TO_plus_max}{Maximum weighted taxonomic distance across
#'       iterations}
#'     \item{uTO_max}{Maximum unweighted taxonomic diversity across
#'       iterations}
#'     \item{TO_max}{Maximum weighted taxonomic diversity across iterations}
#'     \item{uTO_plus_det}{Deterministic uTO+ (first iteration, same as
#'       [ozkan_pto()])}
#'     \item{TO_plus_det}{Deterministic TO+ (first iteration)}
#'     \item{uTO_det}{Deterministic uTO (first iteration)}
#'     \item{TO_det}{Deterministic TO (first iteration)}
#'     \item{n_iter}{Number of iterations performed}
#'     \item{species_status}{Named logical vector from jackknife
#'       (`TRUE` = happy)}
#'     \item{jackknife}{Full jackknife result from
#'       [ozkan_pto_jackknife()]}
#'     \item{n_positive}{Number of iterations with positive uTO+}
#'     \item{iteration_results}{Data frame with all iteration results}
#'   }
#'
#' @details
#' The algorithm follows the Excel macro's Islem 1 + Islem 2 logic:
#' \enumerate{
#'   \item Run jackknife ([ozkan_pto_jackknife()]) to classify each
#'     species as happy or unhappy.
#'   \item Iteration 1: Use the original community (deterministic).
#'   \item Iterations 2..n_iter: For each species:
#'     \itemize{
#'       \item Unhappy species (AA = 0): always included with original
#'         abundance.
#'       \item Happy species (AA > 0): randomly included (50\% probability)
#'         or excluded. Uses `RANDBETWEEN(0,1) * abundance`.
#'     }
#'   \item Return the maximum of each component across all iterations.
#' }
#'
#' @references
#' Ozkan, K. (2018). A new proposed measure for estimating taxonomic
#' diversity. Turkish Journal of Forestry, 19(4), 336-346.
#' DOI: 10.18182/tjf.441061
#'
#' @seealso [ozkan_pto_jackknife()] for the jackknife step,
#'   [ozkan_pto_sensitivity()] for Run 3,
#'   [ozkan_pto_full()] for the full pipeline.
#'
#' @examples
#' comm <- c(sp1 = 4, sp2 = 2, sp3 = 3, sp4 = 1, sp5 = 2)
#' tax <- data.frame(
#'   Species = paste0("sp", 1:5),
#'   Genus   = c("G1", "G1", "G1", "G2", "G2"),
#'   Family  = c("F1", "F1", "F1", "F1", "F1"),
#'   stringsAsFactors = FALSE
#' )
#' set.seed(42)
#' result <- ozkan_pto_resample(comm, tax, n_iter = 101)
#' result$species_status  # Happy/unhappy classification
#'
#' @export
ozkan_pto_resample <- function(community, tax_tree, n_iter = 101L,
                                seed = NULL) {
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
  if (!is.numeric(n_iter) || length(n_iter) != 1 || n_iter < 101) {
    stop("'n_iter' must be a single integer >= 101.", call. = FALSE)
  }
  n_iter <- as.integer(n_iter)

  # Filter to non-zero abundances
  community <- community[community > 0]
  if (length(community) < 2) {
    stop("Need at least 2 species with positive abundance for resampling.",
         call. = FALSE)
  }

  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }

  species_names <- names(community)
  n_species <- length(community)

  # --- Step 1: Jackknife (Islem 1) ---
  jk <- ozkan_pto_jackknife(community, tax_tree)
  is_happy <- jk$species_status

  # --- Step 2: Stochastic resampling (Islem 2) ---
  results <- matrix(NA_real_, nrow = n_iter, ncol = 4)
  colnames(results) <- c("uTO_plus", "TO_plus", "uTO", "TO")

  # Iteration 1: deterministic (same as jackknife full result)
  det <- jk$full_result
  results[1, ] <- c(unname(det$uTO_plus), unname(det$TO_plus),
                     unname(det$uTO), unname(det$TO))

  # Iterations 2..n_iter: stochastic resampling
  n_happy <- sum(is_happy)

  for (iter in 2:n_iter) {
    resampled <- community

    # Happy species: 50% random inclusion (RANDBETWEEN(0,1) * abundance)
    # Unhappy species: always included (no change needed)
    if (n_happy > 0) {
      include_happy <- sample(c(0L, 1L), n_happy, replace = TRUE)
      resampled[is_happy] <- community[is_happy] * include_happy
    }

    # Keep only non-zero species
    resampled <- resampled[resampled > 0]

    if (length(resampled) < 2) {
      results[iter, ] <- c(0, 0, 0, 0)
      next
    }

    iter_result <- ozkan_pto(resampled, tax_tree)
    results[iter, ] <- c(unname(iter_result$uTO_plus),
                          unname(iter_result$TO_plus),
                          unname(iter_result$uTO),
                          unname(iter_result$TO))
  }

  # --- Compute maximums across all iterations ---
  max_vals <- apply(results, 2, max, na.rm = TRUE)

  # Build iteration results data frame
  iter_df <- as.data.frame(results)
  iter_df$iteration <- seq_len(n_iter)
  iter_df <- iter_df[, c("iteration", "uTO_plus", "TO_plus", "uTO", "TO")]

  # Count iterations with positive uTO+ (needed for Run 3 probability)
  n_positive <- sum(results[, "uTO_plus"] > 0, na.rm = TRUE)

  result <- list(
    uTO_plus_max      = unname(max_vals["uTO_plus"]),
    TO_plus_max       = unname(max_vals["TO_plus"]),
    uTO_max           = unname(max_vals["uTO"]),
    TO_max            = unname(max_vals["TO"]),
    uTO_plus_det      = unname(results[1, "uTO_plus"]),
    TO_plus_det       = unname(results[1, "TO_plus"]),
    uTO_det           = unname(results[1, "uTO"]),
    TO_det            = unname(results[1, "TO"]),
    n_iter            = n_iter,
    species_status    = is_happy,
    jackknife         = jk,
    n_positive        = n_positive,
    iteration_results = iter_df
  )
  class(result) <- "ozkan_pto_resample"
  return(result)
}


#' Sensitivity Analysis of Ozkan's pTO Index (Islem 3 / Run 3)
#'
#' Implements the sensitivity analysis procedure from Ozkan's Excel macro
#' (Islem 3). Uses the jackknife results from Run 2 to apply species-specific
#' inclusion probabilities: unhappy species get \eqn{P = (S-1)/S}, happy
#' species get a data-driven probability derived from Run 2 iteration results.
#'
#' @param community A named numeric vector of species abundances.
#' @param tax_tree A data frame with taxonomic hierarchy.
#' @param run2_result The result from [ozkan_pto_resample()].
#' @param n_iter Number of iterations (default: same as Run 2).
#' @param seed Optional random seed for reproducibility.
#'
#' @return A named list with components:
#'   \describe{
#'     \item{uTO_plus_max}{Maximum uTO+ across Run 1, 2, and 3}
#'     \item{TO_plus_max}{Maximum TO+ across all runs}
#'     \item{uTO_max}{Maximum uTO across all runs}
#'     \item{TO_max}{Maximum TO across all runs}
#'     \item{run3_uTO_plus_max}{Maximum uTO+ from Run 3 only}
#'     \item{run3_TO_plus_max}{Maximum TO+ from Run 3 only}
#'     \item{run3_uTO_max}{Maximum uTO from Run 3 only}
#'     \item{run3_TO_max}{Maximum TO from Run 3 only}
#'     \item{n_iter}{Number of iterations performed}
#'     \item{species_probs}{Named numeric vector of inclusion probabilities}
#'     \item{prob_happy}{Probability used for happy species}
#'     \item{prob_unhappy}{Probability used for unhappy species}
#'     \item{iteration_results}{Data frame with all Run 3 iteration results}
#'   }
#'
#' @details
#' The algorithm follows the Excel macro's Islem 3 logic:
#'
#' For each species, the inclusion probability depends on its jackknife
#' classification from Islem 1:
#'
#' \itemize{
#'   \item \strong{Unhappy species} (AA = 0, non-contributing): Included with
#'     probability \eqn{(S-1)/S}, where S is total species count. In the Excel
#'     formula: `IF(RANDBETWEEN(1, S) > 1, H2, 0)`.
#'   \item \strong{Happy species} (AA > 0, contributing): Included with
#'     probability derived from Run 2 results. In the Excel formula:
#'     `IF(L25 >= RANDBETWEEN(0, K22), H2, 0)`, where L25 is a
#'     summary score from Run 2 and K22 is the iteration count.
#' }
#'
#' The happy species probability is computed as:
#' \deqn{P_{happy} = \frac{\max(0, N_{positive} - S) + 1}{N_{iter} + 1}}
#'
#' where \eqn{N_{positive}} is the number of Run 2 iterations that produced
#' a positive uTO+ value and S is the species count.
#'
#' The maximum pTO across all three runs (Run 1, 2, 3) is the final result.
#'
#' @references
#' Ozkan, K. (2018). A new proposed measure for estimating taxonomic
#' diversity. Turkish Journal of Forestry, 19(4), 336-346.
#' DOI: 10.18182/tjf.441061
#'
#' @seealso [ozkan_pto_resample()] for Run 2,
#'   [ozkan_pto_full()] for the full pipeline.
#'
#' @examples
#' comm <- c(sp1 = 4, sp2 = 2, sp3 = 3, sp4 = 1, sp5 = 2)
#' tax <- data.frame(
#'   Species = paste0("sp", 1:5),
#'   Genus   = c("G1", "G1", "G1", "G2", "G2"),
#'   Family  = c("F1", "F1", "F1", "F1", "F1"),
#'   stringsAsFactors = FALSE
#' )
#' set.seed(42)
#' run2 <- ozkan_pto_resample(comm, tax, n_iter = 101)
#' ozkan_pto_sensitivity(comm, tax, run2, n_iter = 101)
#'
#' @export
ozkan_pto_sensitivity <- function(community, tax_tree, run2_result,
                                   n_iter = NULL, seed = NULL) {
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
  if (!is.list(run2_result) || is.null(run2_result$iteration_results)) {
    stop("'run2_result' must be the output of ozkan_pto_resample().",
         call. = FALSE)
  }
  if (is.null(run2_result$species_status)) {
    stop("'run2_result' must contain 'species_status' from jackknife. ",
         "Please re-run ozkan_pto_resample().", call. = FALSE)
  }

  # Use same n_iter as Run 2 if not specified
  if (is.null(n_iter)) {
    n_iter <- run2_result$n_iter
  }
  n_iter <- as.integer(n_iter)

  # Filter to non-zero
  community <- community[community > 0]
  if (length(community) < 2) {
    stop("Need at least 2 species with positive abundance.", call. = FALSE)
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  species_names <- names(community)
  n_species <- length(community)

  # --- Species-specific inclusion probabilities ---
  is_happy <- run2_result$species_status

  # Unhappy species: P = (S-1)/S
  # Excel formula: IF(RANDBETWEEN(1, S) > 1, H2, 0) -> P = (S-1)/S
  prob_unhappy <- (n_species - 1) / n_species

  # Happy species: P based on Run 2 summary
  # Excel formula: IF(L25 >= RANDBETWEEN(0, K22), H2, 0)
  # L25 = n_positive - S (frozen at end of Run 2)
  # K22 = n_iter (Run 2 iteration count)
  # P = (L25 + 1) / (K22 + 1) when L25 >= 0, else P = 0
  n_positive <- run2_result$n_positive
  if (is.null(n_positive)) {
    # Fallback: count from iteration results
    n_positive <- sum(run2_result$iteration_results$uTO_plus > 0,
                      na.rm = TRUE)
  }

  L25_equiv <- n_positive - n_species
  K22_equiv <- run2_result$n_iter

  if (L25_equiv >= 0) {
    prob_happy <- (L25_equiv + 1) / (K22_equiv + 1)
  } else {
    prob_happy <- 0
  }

  # Build species probability vector
  species_probs <- ifelse(is_happy, prob_happy, prob_unhappy)
  names(species_probs) <- species_names

  # --- Resampling loop ---
  results <- matrix(NA_real_, nrow = n_iter, ncol = 4)
  colnames(results) <- c("uTO_plus", "TO_plus", "uTO", "TO")

  # Iteration 1: deterministic
  det_result <- ozkan_pto(community, tax_tree)
  results[1, ] <- c(unname(det_result$uTO_plus), unname(det_result$TO_plus),
                     unname(det_result$uTO), unname(det_result$TO))

  # Iterations 2..n_iter: probability-weighted resampling
  for (iter in 2:n_iter) {
    include <- stats::runif(n_species) < species_probs
    resampled <- community * as.integer(include)
    resampled <- resampled[resampled > 0]

    if (length(resampled) < 2) {
      results[iter, ] <- c(0, 0, 0, 0)
      next
    }

    iter_result <- ozkan_pto(resampled, tax_tree)
    results[iter, ] <- c(unname(iter_result$uTO_plus),
                          unname(iter_result$TO_plus),
                          unname(iter_result$uTO),
                          unname(iter_result$TO))
  }

  # --- Compute Run 3 maximums ---
  run3_max <- apply(results, 2, max, na.rm = TRUE)

  # --- Overall maximums across Run 1, 2, 3 ---
  overall_max <- pmax(
    c(run2_result$uTO_plus_max, run2_result$TO_plus_max,
      run2_result$uTO_max, run2_result$TO_max),
    run3_max
  )
  names(overall_max) <- c("uTO_plus", "TO_plus", "uTO", "TO")

  # Build iteration results data frame
  iter_df <- as.data.frame(results)
  iter_df$iteration <- seq_len(n_iter)
  iter_df <- iter_df[, c("iteration", "uTO_plus", "TO_plus", "uTO", "TO")]

  result <- list(
    uTO_plus_max      = unname(overall_max["uTO_plus"]),
    TO_plus_max       = unname(overall_max["TO_plus"]),
    uTO_max           = unname(overall_max["uTO"]),
    TO_max            = unname(overall_max["TO"]),
    run3_uTO_plus_max = unname(run3_max["uTO_plus"]),
    run3_TO_plus_max  = unname(run3_max["TO_plus"]),
    run3_uTO_max      = unname(run3_max["uTO"]),
    run3_TO_max       = unname(run3_max["TO"]),
    n_iter            = n_iter,
    species_probs     = species_probs,
    prob_happy        = prob_happy,
    prob_unhappy      = prob_unhappy,
    iteration_results = iter_df
  )
  class(result) <- "ozkan_pto_sensitivity"
  return(result)
}


#' Full Ozkan pTO Pipeline (Islem 1 + 2 + 3)
#'
#' Runs the complete Ozkan taxonomic diversity analysis pipeline:
#' jackknife (Islem 1), stochastic resampling (Islem 2), and sensitivity
#' analysis (Islem 3), returning the maximum values across all three runs.
#' This is equivalent to running all three steps in the Excel macro
#' sequentially.
#'
#' @param community A named numeric vector of species abundances.
#'   Names must match the first column of `tax_tree`.
#' @param tax_tree A data frame with taxonomic hierarchy. First column
#'   is species names, subsequent columns are taxonomic ranks.
#' @param n_iter Number of stochastic iterations for Run 2 and Run 3
#'   (default: 101, minimum: 101).
#' @param seed Optional random seed for reproducibility. If provided,
#'   Run 2 uses this seed and Run 3 uses `seed + 1` to ensure
#'   independent randomness.
#'
#' @return A named list with components:
#'   \describe{
#'     \item{uTO_plus}{Final maximum uTO+ across all 3 runs}
#'     \item{TO_plus}{Final maximum TO+ across all 3 runs}
#'     \item{uTO}{Final maximum uTO across all 3 runs}
#'     \item{TO}{Final maximum TO across all 3 runs}
#'     \item{run1}{Deterministic pTO result (from [ozkan_pto()])}
#'     \item{run2}{Full Run 2 result (from [ozkan_pto_resample()])}
#'     \item{run3}{Full Run 3 result (from [ozkan_pto_sensitivity()])}
#'     \item{jackknife}{Jackknife result with species classifications}
#'   }
#'
#' @details
#' This function implements the full Excel macro pipeline in a single call:
#' \enumerate{
#'   \item \strong{Islem 1}: Leave-one-out jackknife to identify
#'     contributing (happy) vs non-contributing (unhappy) species, plus
#'     deterministic pTO calculation.
#'   \item \strong{Islem 2}: Stochastic resampling -- unhappy species are
#'     always included, happy species get 50\% random inclusion.
#'   \item \strong{Islem 3}: Sensitivity analysis -- unhappy species get
#'     \eqn{P = (S-1)/S}, happy species get a data-driven probability.
#'   \item \strong{Final result}: Maximum values across all three runs.
#' }
#'
#' @references
#' Ozkan, K. (2018). A new proposed measure for estimating taxonomic
#' diversity. Turkish Journal of Forestry, 19(4), 336-346.
#' DOI: 10.18182/tjf.441061
#'
#' @seealso [ozkan_pto()] for deterministic calculation only,
#'   [ozkan_pto_resample()] for Run 2 only,
#'   [ozkan_pto_sensitivity()] for Run 3 only,
#'   [ozkan_pto_jackknife()] for jackknife only.
#'
#' @examples
#' comm <- c(sp1 = 4, sp2 = 2, sp3 = 3, sp4 = 1, sp5 = 2)
#' tax <- data.frame(
#'   Species = paste0("sp", 1:5),
#'   Genus   = c("G1", "G1", "G1", "G2", "G2"),
#'   Family  = c("F1", "F1", "F1", "F1", "F1"),
#'   stringsAsFactors = FALSE
#' )
#' set.seed(42)
#' result <- ozkan_pto_full(comm, tax, n_iter = 101)
#' result$uTO_plus  # Final maximum uTO+
#' result$TO_plus   # Final maximum TO+
#'
#' @export
ozkan_pto_full <- function(community, tax_tree, n_iter = 101L,
                            seed = NULL) {
  # --- Input validation (basic; ozkan_pto_resample does full validation) ---
  if (!is.numeric(n_iter) || length(n_iter) != 1 || n_iter < 101) {
    stop("'n_iter' must be a single integer >= 101.", call. = FALSE)
  }
  n_iter <- as.integer(n_iter)

  # --- Run 2 (includes Run 1 jackknife internally) ---
  run2 <- ozkan_pto_resample(community, tax_tree,
                              n_iter = n_iter, seed = seed)

  # --- Run 3 (uses Run 2 results) ---
  # Use seed + 1 for Run 3 to ensure independent randomness
  seed3 <- if (!is.null(seed)) seed + 1L else NULL
  run3 <- ozkan_pto_sensitivity(community, tax_tree, run2,
                                 n_iter = n_iter, seed = seed3)

  # --- Final result: MAX across all 3 runs ---
  list(
    uTO_plus  = run3$uTO_plus_max,
    TO_plus   = run3$TO_plus_max,
    uTO       = run3$uTO_max,
    TO        = run3$TO_max,
    run1      = run2$jackknife$full_result,
    run2      = run2,
    run3      = run3,
    jackknife = run2$jackknife
  )
}
