#' Taxonomic Diversity Rarefaction
#'
#' Computes rarefaction curves for taxonomic diversity indices by
#' subsampling individuals from the community at increasing sample sizes.
#' Uses bootstrap resampling to estimate expected diversity and
#' confidence intervals at each sample size.
#'
#' @param community A named numeric vector of species abundances (integers).
#'   Names must match the first column of `tax_tree`.
#' @param tax_tree A data frame with taxonomic hierarchy. First column
#'   is species names, subsequent columns are taxonomic ranks.
#' @param index Which index to rarefy. One of `"shannon"`, `"simpson"`,
#'   `"uTO"`, `"TO"`, `"uTO_plus"`, `"TO_plus"`, `"avtd"`, `"species"`
#'   (default: `"shannon"`).
#' @param steps Number of sample-size steps along the curve (default: 20).
#' @param n_boot Number of bootstrap replicates per step (default: 100).
#' @param ci Confidence interval width (default: 0.95).
#' @param seed Optional random seed for reproducibility (default: NULL).
#' @param correction Bias correction for the Shannon index. One of
#'   `"none"` (default), `"miller_madow"`, `"grassberger"`, or
#'   `"chao_shen"`. Only used when `index = "shannon"`. Passed to
#'   [shannon()]. See [shannon()] for details.
#' @param parallel Logical. If `TRUE`, use parallel processing to speed up
#'   bootstrap resampling across sample sizes. Default `FALSE`.
#' @param n_cores Number of CPU cores to use when `parallel = TRUE`. Default
#'   `NULL` uses up to 2 cores (CRAN policy limit).
#'
#' @return A data frame with columns:
#'   \describe{
#'     \item{sample_size}{Number of individuals in the subsample}
#'     \item{mean}{Mean index value across bootstrap replicates}
#'     \item{lower}{Lower bound of the confidence interval}
#'     \item{upper}{Upper bound of the confidence interval}
#'     \item{sd}{Standard deviation across replicates}
#'   }
#'
#' @details
#' The algorithm works as follows:
#' \enumerate{
#'   \item Expand the abundance vector into an individual-level vector
#'     (e.g., c(sp1=3, sp2=2) becomes c("sp1","sp1","sp1","sp2","sp2")).
#'   \item For each sample size (from min to total N), draw `n_boot`
#'     random subsamples without replacement.
#'   \item For each subsample, count species abundances and compute
#'     the chosen diversity index.
#'   \item Return the mean and confidence interval at each step.
#' }
#'
#' @references
#' Gotelli, N.J. & Colwell, R.K. (2001). Quantifying biodiversity:
#' procedures and pitfalls in the measurement and comparison of species
#' richness. Ecology Letters, 4, 379-391.
#'
#' Ozkan, K. (2018). A new proposed measure for estimating taxonomic
#' diversity. Turkish Journal of Forestry, 19(4), 336-346.
#'
#' Ozkan, K. & Mert, A. (2022). Comparisons of Deng entropy-based
#' taxonomic diversity measures with the other diversity measures and
#' introduction to the new proposed (reinforced) estimators. FORESTIST,
#' 72(2). DOI: 10.5152/forestist.2021.21025
#'
#' @seealso [plot_rarefaction()] for visualising the rarefaction curve,
#'   [ozkan_pto()] for full pTO calculation, [shannon()] and [simpson()]
#'   for classical indices, [avtd()] for average taxonomic distinctness.
#'
#' @examples
#' comm <- c(sp1 = 10, sp2 = 5, sp3 = 8, sp4 = 2, sp5 = 3)
#' tax <- data.frame(
#'   Species = paste0("sp", 1:5),
#'   Genus   = c("G1", "G1", "G2", "G2", "G3"),
#'   Family  = c("F1", "F1", "F1", "F2", "F2"),
#'   stringsAsFactors = FALSE
#' )
#' rarefaction_taxonomic(comm, tax, index = "shannon", n_boot = 50)
#'
#' @export
rarefaction_taxonomic <- function(community, tax_tree,
                                  index = c("shannon", "simpson", "species",
                                            "uTO", "TO", "uTO_plus",
                                            "TO_plus", "avtd"),
                                  steps = 20, n_boot = 100,
                                  ci = 0.95, seed = NULL,
                                  correction = c("none", "miller_madow",
                                                 "grassberger",
                                                 "chao_shen"),
                                  parallel = FALSE,
                                  n_cores = NULL) {

  index <- match.arg(index)
  correction <- match.arg(correction)

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
  if (steps < 2) {
    stop("'steps' must be at least 2.", call. = FALSE)
  }
  if (n_boot < 2) {
    stop("'n_boot' must be at least 2.", call. = FALSE)
  }
  if (ci <= 0 || ci >= 1) {
    stop("'ci' must be between 0 and 1.", call. = FALSE)
  }

  # Filter zero abundances
  community <- community[community > 0]

  # Round to integers for individual-level resampling
  community <- round(community)
  total_n <- sum(community)

  if (total_n < 2) {
    stop("Total abundance must be at least 2 for rarefaction.", call. = FALSE)
  }

  # Set seed if provided
  if (!is.null(seed)) set.seed(seed)

  # Expand to individual-level vector
  individuals <- rep(names(community), times = community)

  # Determine sample sizes
  min_n <- max(2, min(community[community > 0]))
  sample_sizes <- unique(round(seq(min_n, total_n, length.out = steps)))

  # CI quantiles
  alpha <- (1 - ci) / 2
  q_lower <- alpha
  q_upper <- 1 - alpha

  # --- Helper: compute index from a subsample abundance vector ---
  compute_index <- function(sub_comm) {
    if (length(sub_comm) < 2) return(0)

    switch(index,
      species = length(sub_comm),
      shannon = shannon(sub_comm, correction = correction),
      simpson = simpson(sub_comm),
      uTO = {
        result <- ozkan_pto(sub_comm, tax_tree)
        result$uTO
      },
      TO = {
        result <- ozkan_pto(sub_comm, tax_tree)
        result$TO
      },
      uTO_plus = {
        result <- ozkan_pto(sub_comm, tax_tree)
        result$uTO_plus
      },
      TO_plus = {
        result <- ozkan_pto(sub_comm, tax_tree)
        result$TO_plus
      },
      avtd = {
        sp_names <- names(sub_comm)
        if (length(unique(sp_names)) < 2) return(0)
        avtd(sp_names, tax_tree)
      }
    )
  }

  # --- Worker: bootstrap at one sample size ---
  boot_one_step <- function(n) {
    boot_values <- numeric(n_boot)
    for (b in seq_len(n_boot)) {
      sub_individuals <- sample(individuals, size = n, replace = FALSE)
      sub_comm <- table(sub_individuals)
      sub_comm <- setNames(as.numeric(sub_comm), names(sub_comm))
      val <- tryCatch(compute_index(sub_comm), error = function(e) NA_real_)
      boot_values[b] <- val
    }
    boot_clean <- boot_values[!is.na(boot_values)]
    if (length(boot_clean) == 0) boot_clean <- 0
    data.frame(
      sample_size = n,
      mean = mean(boot_clean),
      lower = stats::quantile(boot_clean, q_lower, names = FALSE),
      upper = stats::quantile(boot_clean, q_upper, names = FALSE),
      sd = stats::sd(boot_clean)
    )
  }

  # --- Run: parallel or sequential ---
  if (isTRUE(parallel)) {
    n_cores_use <- if (is.null(n_cores)) min(2L, parallel::detectCores()) else max(1L, as.integer(n_cores))

    if (.Platform$OS.type == "windows") {
      cl <- parallel::makeCluster(n_cores_use)
      on.exit(parallel::stopCluster(cl), add = TRUE)
      parallel::clusterExport(cl, varlist = c(
        "individuals", "n_boot", "q_lower", "q_upper",
        "compute_index"
      ), envir = environment())
      ns <- asNamespace(utils::packageName())
      fns_to_export <- c("shannon", "simpson", "ozkan_pto", "avtd")
      fns_to_export <- fns_to_export[vapply(fns_to_export, exists, logical(1), where = ns)]
      if (length(fns_to_export) > 0) {
        parallel::clusterExport(cl, varlist = fns_to_export, envir = ns)
      }
      if (!is.null(seed)) parallel::clusterSetRNGStream(cl, seed)
      result_list <- parallel::parLapply(cl, sample_sizes, boot_one_step)
    } else {
      result_list <- parallel::mclapply(sample_sizes, boot_one_step,
                                         mc.cores = n_cores_use,
                                         mc.set.seed = TRUE)
    }
  } else {
    result_list <- lapply(sample_sizes, boot_one_step)
  }

  results <- do.call(rbind, result_list)
  rownames(results) <- NULL

  # Add metadata as attributes
  attr(results, "index") <- index
  attr(results, "total_n") <- total_n
  attr(results, "n_boot") <- n_boot
  attr(results, "ci") <- ci
  class(results) <- c("rarefaction_taxonomic", "data.frame")

  return(results)
}
