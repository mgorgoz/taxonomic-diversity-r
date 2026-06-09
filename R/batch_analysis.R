#' Batch Analysis from a Single Data Frame
#'
#' Computes selected diversity indices for one or more sample sites from a
#' single data frame (e.g., imported from Excel). The function automatically
#' detects the site column, taxonomic columns, and abundance column, splits
#' the data by site, and returns a summary data frame with species count and
#' diversity indices per site.
#'
#' When no site column is present (or all values are identical), the entire
#' data set is treated as a single community.
#'
#' Three groups of indices are available via the `indices` parameter:
#' \describe{
#'   \item{\code{"classical"}}{Shannon-Wiener entropy and Gini-Simpson index.
#'     These are species-level diversity measures that do not use taxonomic
#'     hierarchy.}
#'   \item{\code{"clarke_warwick"}}{Delta, Delta*, AvTD, and VarTD.
#'     Taxonomy-aware distinctness measures from Clarke & Warwick (1998).}
#'   \item{\code{"ozkan_pto"}}{Deng entropy-based taxonomic diversity and
#'     distance (uTO, TO, uTO+, TO+) from Ozkan (2018). When \code{full = TRUE},
#'     also computes max values across Run 1+2+3. The Westhoff-Maarel
#'     cover-abundance scale (integer values 1-9) is recommended for
#'     compatibility with the original paper, but any positive numeric
#'     abundance values are accepted.}
#' }
#'
#' @param data A data frame containing species data. Must include at minimum
#'   a species column, at least one taxonomic rank column, and an abundance
#'   column. Optionally includes a site/plot column for multi-site analysis.
#' @param indices Character vector specifying which index groups to compute.
#'   One or more of \code{"classical"} (Shannon, Simpson),
#'   \code{"clarke_warwick"} (Delta, Delta*, AvTD, VarTD), and
#'   \code{"ozkan_pto"} (uTO, TO, uTO+, TO+ with optional Run 2+3).
#'   Default is all three groups. Unambiguous abbreviations are allowed (e.g.,
#'   \code{"clas"} for classical, \code{"clark"} for clarke_warwick,
#'   \code{"oz"} for ozkan_pto). Note that \code{"cl"} and \code{"cla"} are
#'   ambiguous and will produce an error.
#' @param site_column Character string specifying the name of the site column.
#'   If \code{NULL} (default), the function searches for columns named
#'   \code{"Site"}, \code{"site"}, \code{"Plot"}, or \code{"plot"}.
#'   If no such column is found, all data is treated as a single site.
#' @param tax_columns Character vector specifying the names of the taxonomic
#'   columns (from Species to highest rank). If \code{NULL} (default), the
#'   function auto-detects columns named \code{"Species"}, \code{"Genus"},
#'   \code{"Family"}, \code{"Order"}, \code{"Class"}, \code{"Phylum"}, and
#'   \code{"Kingdom"} (case-insensitive).
#' @param abundance_column Character string specifying the name of the
#'   abundance column. Default is \code{"Abundance"} (case-insensitive match).
#' @param correction Bias correction for the Shannon index. One of
#'   `"none"` (default), `"miller_madow"`, `"grassberger"`, or
#'   `"chao_shen"`. Passed to [shannon()]. See [shannon()] for details.
#'   Ignored when `"classical"` is not in `indices`.
#' @param full Logical. If `TRUE` (default), run the full Ozkan pipeline
#'   (Run 1+2+3) using [ozkan_pto_full()] instead of deterministic-only
#'   [pto_components()]. This produces max values across all three runs,
#'   matching the Excel macro output. Set to `FALSE` for deterministic
#'   Run 1 only (faster but incomplete). Ignored when `"ozkan_pto"` is not
#'   in `indices`.
#' @param n_iter Number of stochastic iterations for Run 2 and Run 3 when
#'   `full = TRUE`. Default `101`. Ignored when `full = FALSE` or when
#'   `"ozkan_pto"` is not in `indices`.
#' @param seed Random seed for reproducibility when `full = TRUE`.
#'   Default `42`. Set to `NULL` for non-deterministic results.
#'   Ignored when `full = FALSE` or when `"ozkan_pto"` is not in `indices`.
#' @param parallel Logical. If `TRUE`, use parallel processing to compute
#'   indices for multiple sites concurrently. Default `FALSE`.
#' @param n_cores Number of CPU cores to use when `parallel = TRUE`. Default
#'   `NULL` uses up to 2 cores (CRAN policy limit).
#' @param progress Logical. If `TRUE` (default), display a progress bar
#'   during sequential computation. Ignored when `parallel = TRUE`.
#'   Set to `FALSE` to suppress progress output.
#' @param progress_fn Optional callback function for custom progress reporting
#'   (e.g. Shiny). When provided, it is called after each site completes with
#'   named arguments `i` (current site index), `n` (total sites), and
#'   `site` (site name). Useful for integrating with `shiny::withProgress()`.
#'
#' @return A data frame with one row per site. Columns always include
#'   \code{Site} and \code{N_Species}. Additional columns depend on the
#'   `indices` parameter:
#'   \itemize{
#'     \item \code{"classical"}: \code{Shannon}, \code{Simpson}
#'     \item \code{"clarke_warwick"}: \code{Delta}, \code{Delta_star},
#'       \code{AvTD}, \code{VarTD}
#'     \item \code{"ozkan_pto"}: \code{uTO}, \code{TO}, \code{uTO_plus},
#'       \code{TO_plus}, \code{uTO_max}, \code{TO_max},
#'       \code{uTO_plus_max}, \code{TO_plus_max}
#'   }
#'
#' @examples
#' # All indices (default)
#' \donttest{
#' df <- data.frame(
#'   Species   = c("sp1", "sp2", "sp3", "sp4"),
#'   Genus     = c("G1", "G1", "G2", "G2"),
#'   Family    = c("F1", "F1", "F1", "F2"),
#'   Order     = c("O1", "O1", "O1", "O1"),
#'   Abundance = c(4, 2, 3, 1),
#'   stringsAsFactors = FALSE
#' )
#' batch_analysis(df)
#'
#' # Only classical indices (fast)
#' batch_analysis(df, indices = "classical")
#'
#' # Classical + Clarke & Warwick (no pTO)
#' batch_analysis(df, indices = c("classical", "clarke_warwick"))
#'
#' # Only Ozkan pTO, deterministic Run 1
#' batch_analysis(df, indices = "ozkan_pto", full = FALSE)
#' }
#'
#' @seealso \code{\link{compare_indices}} for analysis with pre-built community
#'   vectors, \code{\link{build_tax_tree}} for building taxonomic trees manually,
#'   \code{\link{ozkan_pto_full}} for the full 3-run pipeline on a single community.
#'
#' @export
batch_analysis <- function(data,
                           indices = c("classical", "clarke_warwick",
                                       "ozkan_pto"),
                           site_column = NULL,
                           tax_columns = NULL,
                           abundance_column = "Abundance",
                           correction = c("none", "miller_madow",
                                          "grassberger", "chao_shen"),
                           full = TRUE,
                           n_iter = 101L,
                           seed = 42L,
                           parallel = FALSE,
                           n_cores = NULL,
                           progress = TRUE,
                           progress_fn = NULL) {
  correction <- match.arg(correction)

  # --- Validate indices parameter ---
  valid_indices <- c("classical", "clarke_warwick", "ozkan_pto")
  indices <- unique(tolower(indices))
  # Allow partial matching
  indices <- vapply(indices, function(x) {
    m <- pmatch(x, valid_indices)
    if (is.na(m)) {
      stop("Unknown index group: '", x, "'. Valid groups: ",
           paste(valid_indices, collapse = ", "), call. = FALSE)
    }
    valid_indices[m]
  }, character(1), USE.NAMES = FALSE)
  indices <- unique(indices)

  do_classical     <- "classical" %in% indices
  do_clarke        <- "clarke_warwick" %in% indices
  do_pto           <- "ozkan_pto" %in% indices

  # --- Input validation ---
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame.", call. = FALSE)
  }

  if (nrow(data) == 0) {
    stop("'data' has no rows.", call. = FALSE)
  }

  col_names <- names(data)
  col_lower <- tolower(col_names)

  # --- Find abundance column ---
  if (is.null(abundance_column)) {
    stop("'abundance_column' cannot be NULL.", call. = FALSE)
  }

  abd_idx <- which(col_lower == tolower(abundance_column))
  if (length(abd_idx) == 0) {
    stop("Abundance column '", abundance_column,
         "' not found in data. Available columns: ",
         paste(col_names, collapse = ", "), call. = FALSE)
  }
  abd_col <- col_names[abd_idx[1]]

  if (!is.numeric(data[[abd_col]])) {
    stop("Abundance column '", abd_col, "' must be numeric.", call. = FALSE)
  }

  # --- Find site column ---
  site_col <- NULL
  site_auto_names <- c("site", "plot")


  if (!is.null(site_column)) {
    # User specified
    site_idx <- which(col_lower == tolower(site_column))
    if (length(site_idx) == 0) {
      stop("Site column '", site_column,
           "' not found in data. Available columns: ",
           paste(col_names, collapse = ", "), call. = FALSE)
    }
    site_col <- col_names[site_idx[1]]
  } else {
    # Auto-detect
    for (sn in site_auto_names) {
      idx <- which(col_lower == sn)
      if (length(idx) > 0) {
        site_col <- col_names[idx[1]]
        break
      }
    }
  }

  # Check if site column has meaningful values
  if (!is.null(site_col)) {
    site_vals <- data[[site_col]]
    # If all NA or all empty string, treat as single site
    if (all(is.na(site_vals)) || all(trimws(as.character(site_vals)) == "")) {
      site_col <- NULL
    }
  }

  # --- Find taxonomic columns ---
  known_tax_ranks <- c("species", "genus", "family", "order",
                       "class", "phylum", "kingdom")

  if (!is.null(tax_columns)) {
    # User specified
    missing <- tax_columns[!tolower(tax_columns) %in% col_lower]
    if (length(missing) > 0) {
      stop("Taxonomic columns not found: ",
           paste(missing, collapse = ", "), call. = FALSE)
    }
    tax_cols <- col_names[match(tolower(tax_columns), col_lower)]
  } else {
    # Auto-detect in canonical order
    tax_cols <- character(0)
    for (rank in known_tax_ranks) {
      idx <- which(col_lower == rank)
      if (length(idx) > 0) {
        tax_cols <- c(tax_cols, col_names[idx[1]])
      }
    }
    if (length(tax_cols) < 2) {
      stop("Could not auto-detect taxonomic columns. ",
           "Need at least Species + one rank column. ",
           "Found: ", paste(tax_cols, collapse = ", "),
           ". Use 'tax_columns' parameter to specify manually.", call. = FALSE)
    }
  }

  # Species column is the first taxonomic column
  species_col <- tax_cols[1]
  rank_cols <- tax_cols[-1]  # Genus, Family, Order, ...

  # --- Split by site ---
  if (is.null(site_col)) {
    # Single site
    site_list <- list(All = data)
  } else {
    site_list <- split(data, data[[site_col]])
    # Natural sort: Site1, Site2, ..., Site9, Site10 (not Site1, Site10, Site2)
    sn       <- names(site_list)
    txt_part <- sub("(\\d+)$", "", sn)
    num_part <- suppressWarnings(as.numeric(sub("^\\D*(\\d+)$", "\\1", sn)))
    num_part[is.na(num_part)] <- 0
    site_list <- site_list[order(txt_part, num_part)]
  }

  # --- Worker function: compute indices for one site ---
  compute_site <- function(site_name) {
    site_data <- site_list[[site_name]]

    # Build community vector
    species_names <- as.character(site_data[[species_col]])
    abundances <- site_data[[abd_col]]
    community <- stats::setNames(abundances, species_names)

    # Remove zero/NA abundances
    community <- community[!is.na(community) & community > 0]

    if (length(community) < 2) {
      warning("Site '", site_name, "' has fewer than 2 species with ",
              "positive abundance. Skipping.", call. = FALSE)
      return(NULL)
    }

    # Build tax_tree for this site's species
    sp_present <- names(community)
    tax_rows <- site_data[match(sp_present, as.character(site_data[[species_col]])), ]
    tax_tree_site <- tax_rows[, tax_cols, drop = FALSE]

    for (j in seq_along(tax_tree_site)) {
      tax_tree_site[[j]] <- as.character(tax_tree_site[[j]])
    }
    rownames(tax_tree_site) <- NULL

    # Build result starting with Site and N_Species
    result_row <- data.frame(
      Site      = site_name,
      N_Species = length(community),
      row.names = NULL,
      stringsAsFactors = FALSE
    )

    # Classical indices: Shannon, Simpson
    if (do_classical) {
      result_row$Shannon <- round(shannon(community, correction = correction), 6)
      result_row$Simpson <- round(simpson(community), 6)
    }

    # Clarke & Warwick indices: Delta, Delta*, AvTD, VarTD
    if (do_clarke) {
      result_row$Delta      <- round(delta(community, tax_tree_site), 6)
      result_row$Delta_star <- round(delta_star(community, tax_tree_site), 6)
      sp_names <- names(community)
      result_row$AvTD       <- round(avtd(sp_names, tax_tree_site), 6)
      result_row$VarTD      <- round(vartd(sp_names, tax_tree_site), 6)
    }

    # Ozkan pTO indices
    if (do_pto) {
      if (isTRUE(full)) {
        # Full pipeline: Run 1 + 2 + 3 (matches Excel macro)
        pto_run1 <- pto_components(community, tax_tree_site)
        full_res <- ozkan_pto_full(community, tax_tree_site,
                                    n_iter = n_iter, seed = seed)
        result_row$uTO          <- round(pto_run1["uTO"], 6)
        result_row$TO           <- round(pto_run1["TO"], 6)
        result_row$uTO_plus     <- round(pto_run1["uTO_plus"], 6)
        result_row$TO_plus      <- round(pto_run1["TO_plus"], 6)
        result_row$uTO_max      <- round(full_res$uTO, 6)
        result_row$TO_max       <- round(full_res$TO, 6)
        result_row$uTO_plus_max <- round(full_res$uTO_plus, 6)
        result_row$TO_plus_max  <- round(full_res$TO_plus, 6)
      } else {
        # Deterministic only: Run 1 (fast)
        pto <- pto_components(community, tax_tree_site)
        result_row$uTO          <- round(pto["uTO"], 6)
        result_row$TO           <- round(pto["TO"], 6)
        result_row$uTO_plus     <- round(pto["uTO_plus"], 6)
        result_row$TO_plus      <- round(pto["TO_plus"], 6)
        result_row$uTO_max      <- round(pto["uTO_max"], 6)
        result_row$TO_max       <- round(pto["TO_max"], 6)
        result_row$uTO_plus_max <- round(pto["uTO_plus_max"], 6)
        result_row$TO_plus_max  <- round(pto["TO_plus_max"], 6)
      }
    }

    result_row
  }

  # --- Run: parallel or sequential ---
  site_names <- names(site_list)

  if (isTRUE(parallel) && length(site_names) > 1) {
    n_cores_use <- if (is.null(n_cores)) min(2L, parallel::detectCores()) else max(1L, as.integer(n_cores))
    if (isTRUE(progress) && interactive()) {
      mode_label <- if (isTRUE(full)) "full pipeline" else "Run 1 only"
      message(sprintf("Processing %d sites in parallel (%d cores, %s)...",
                      length(site_names), n_cores_use, mode_label))
    }

    if (.Platform$OS.type == "windows") {
      cl <- parallel::makeCluster(n_cores_use)
      on.exit(parallel::stopCluster(cl), add = TRUE)
      parallel::clusterExport(cl, varlist = c(
        "site_list", "species_col", "abd_col", "tax_cols", "correction",
        "full", "n_iter", "seed"
      ), envir = environment())
      ns <- asNamespace(utils::packageName())
      parallel::clusterExport(cl, varlist = c(
        "shannon", "simpson", "delta", "delta_star",
        "avtd", "vartd", "pto_components", "ozkan_pto",
        "deng_entropy_level", "ozkan_pto_full",
        "ozkan_pto_resample", "ozkan_pto_sensitivity",
        "ozkan_pto_jackknife"
      ), envir = ns)
      results <- parallel::parLapply(cl, site_names, compute_site)
    } else {
      results <- parallel::mclapply(site_names, compute_site,
                                     mc.cores = n_cores_use)
    }
  } else {
    n_sites <- length(site_names)
    show_progress <- isTRUE(progress) && n_sites > 1 && interactive()
    use_callback  <- is.function(progress_fn) && n_sites > 1

    if (show_progress || use_callback) {
      if (show_progress) {
        mode_label <- if (isTRUE(full)) "full pipeline" else "Run 1 only"
        message(sprintf("Processing %d sites (%s)...", n_sites, mode_label))
        pb <- utils::txtProgressBar(min = 0, max = n_sites, style = 3)
        on.exit(close(pb), add = TRUE)
      }
      results <- vector("list", n_sites)
      for (i in seq_along(site_names)) {
        results[[i]] <- compute_site(site_names[i])
        if (show_progress) utils::setTxtProgressBar(pb, i)
        if (use_callback) {
          tryCatch(
            progress_fn(i = i, n = n_sites, site = site_names[i]),
            error = function(e) NULL
          )
        }
      }
    } else {
      results <- lapply(site_names, compute_site)
    }
  }

  # Remove NULLs (skipped sites)
  results <- results[!vapply(results, is.null, logical(1))]

  if (length(results) == 0) {
    stop("No sites had enough species (>= 2) for analysis.", call. = FALSE)
  }

  result_df <- do.call(rbind, results)
  rownames(result_df) <- NULL
  class(result_df) <- c("batch_analysis", "data.frame")

  return(result_df)
}
