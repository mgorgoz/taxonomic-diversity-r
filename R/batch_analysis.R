#' Batch Analysis from a Single Data Frame
#'
#' Computes all diversity indices for one or more sample sites from a single
#' data frame (e.g., imported from Excel). The function automatically detects
#' the site column, taxonomic columns, and abundance column, splits the data
#' by site, and returns a summary data frame with 14 diversity indices per site.
#'
#' When no site column is present (or all values are identical), the entire
#' data set is treated as a single community.
#'
#' The function calculates the following indices per site:
#' \itemize{
#'   \item \strong{Shannon}: Shannon-Wiener entropy (\code{\link{shannon}})
#'   \item \strong{Simpson}: Gini-Simpson index (\code{\link{simpson}})
#'   \item \strong{Delta}: Clarke & Warwick taxonomic diversity (\code{\link{delta}})
#'   \item \strong{Delta_star}: Clarke & Warwick taxonomic distinctness (\code{\link{delta_star}})
#'   \item \strong{AvTD}: Average taxonomic distinctness (\code{\link{avtd}})
#'   \item \strong{VarTD}: Variation in taxonomic distinctness (\code{\link{vartd}})
#'   \item \strong{uTO}: Unweighted taxonomic diversity (Ozkan pTO, all levels)
#'   \item \strong{TO}: Weighted taxonomic diversity (Ozkan pTO, all levels)
#'   \item \strong{uTO_plus}: Unweighted taxonomic distance (Ozkan pTO, all levels)
#'   \item \strong{TO_plus}: Weighted taxonomic distance (Ozkan pTO, all levels)
#'   \item \strong{uTO_max}: Unweighted taxonomic diversity (informative levels only)
#'   \item \strong{TO_max}: Weighted taxonomic diversity (informative levels only)
#'   \item \strong{uTO_plus_max}: Unweighted taxonomic distance (informative levels only)
#'   \item \strong{TO_plus_max}: Weighted taxonomic distance (informative levels only)
#' }
#'
#' @param data A data frame containing species data. Must include at minimum
#'   a species column, at least one taxonomic rank column, and an abundance
#'   column. Optionally includes a site/plot column for multi-site analysis.
#' @param site_column Character string specifying the name of the site column.
#'   If \code{NULL} (default), the function searches for columns named
#'   \code{"Site"}, \code{"site"}, \code{"Alan"}, \code{"alan"}, \code{"Plot"},
#'   or \code{"plot"}. If no such column is found, all data is treated as a
#'   single site.
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
#' @param parallel Logical. If `TRUE`, use parallel processing to compute
#'   indices for multiple sites concurrently. Default `FALSE`.
#' @param n_cores Number of CPU cores to use when `parallel = TRUE`. Default
#'   `NULL` uses `parallel::detectCores() - 1`.
#'
#' @return A data frame with one row per site and columns:
#'   \code{Site}, \code{Shannon}, \code{Simpson}, \code{Delta},
#'   \code{Delta_star}, \code{AvTD}, \code{VarTD}, \code{uTO}, \code{TO},
#'   \code{uTO_plus}, \code{TO_plus}, \code{uTO_max}, \code{TO_max},
#'   \code{uTO_plus_max}, \code{TO_plus_max}.
#'
#' @examples
#' # Single-site data (no Site column)
#' df <- data.frame(
#'   Species   = c("sp1", "sp2", "sp3", "sp4"),
#'   Genus     = c("G1", "G1", "G2", "G2"),
#'   Family    = c("F1", "F1", "F1", "F2"),
#'   Order     = c("O1", "O1", "O1", "O1"),
#'   Abundance = c(10, 20, 15, 5),
#'   stringsAsFactors = FALSE
#' )
#' batch_analysis(df)
#'
#' # Multi-site data (with Site column)
#' df2 <- data.frame(
#'   Site      = c("A", "A", "A", "B", "B", "B"),
#'   Species   = c("sp1", "sp2", "sp3", "sp1", "sp3", "sp4"),
#'   Genus     = c("G1", "G1", "G2", "G1", "G2", "G2"),
#'   Family    = c("F1", "F1", "F1", "F1", "F1", "F2"),
#'   Order     = c("O1", "O1", "O1", "O1", "O1", "O1"),
#'   Abundance = c(10, 20, 15, 5, 25, 10),
#'   stringsAsFactors = FALSE
#' )
#' batch_analysis(df2)
#'
#' @seealso \code{\link{compare_indices}} for analysis with pre-built community
#'   vectors, \code{\link{build_tax_tree}} for building taxonomic trees manually.
#'
#' @export
batch_analysis <- function(data,
                           site_column = NULL,
                           tax_columns = NULL,
                           abundance_column = "Abundance",
                           correction = c("none", "miller_madow",
                                          "grassberger", "chao_shen"),
                           parallel = FALSE,
                           n_cores = NULL) {
  correction <- match.arg(correction)

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
  site_auto_names <- c("site", "alan", "plot")


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

    # Compute all indices
    H <- shannon(community, correction = correction)
    D <- simpson(community)
    Delta_val <- delta(community, tax_tree_site)
    Delta_s <- delta_star(community, tax_tree_site)

    sp_names <- names(community)
    AvTD_val <- avtd(sp_names, tax_tree_site)
    VarTD_val <- vartd(sp_names, tax_tree_site)

    pto <- pto_components(community, tax_tree_site)

    data.frame(
      Site         = site_name,
      Shannon      = round(H, 6),
      Simpson      = round(D, 6),
      Delta        = round(Delta_val, 6),
      Delta_star   = round(Delta_s, 6),
      AvTD         = round(AvTD_val, 6),
      VarTD        = round(VarTD_val, 6),
      uTO          = round(pto["uTO"], 6),
      TO           = round(pto["TO"], 6),
      uTO_plus     = round(pto["uTO_plus"], 6),
      TO_plus      = round(pto["TO_plus"], 6),
      uTO_max      = round(pto["uTO_max"], 6),
      TO_max       = round(pto["TO_max"], 6),
      uTO_plus_max = round(pto["uTO_plus_max"], 6),
      TO_plus_max  = round(pto["TO_plus_max"], 6),
      row.names    = NULL,
      stringsAsFactors = FALSE
    )
  }

  # --- Run: parallel or sequential ---
  site_names <- names(site_list)

  if (isTRUE(parallel) && length(site_names) > 1) {
    n_cores_use <- if (is.null(n_cores)) max(1L, parallel::detectCores() - 1L) else max(1L, as.integer(n_cores))

    if (.Platform$OS.type == "windows") {
      cl <- parallel::makeCluster(n_cores_use)
      on.exit(parallel::stopCluster(cl), add = TRUE)
      parallel::clusterExport(cl, varlist = c(
        "site_list", "species_col", "abd_col", "tax_cols", "correction"
      ), envir = environment())
      ns <- asNamespace(utils::packageName())
      parallel::clusterExport(cl, varlist = c(
        "shannon", "simpson", "delta", "delta_star",
        "avtd", "vartd", "pto_components", "ozkan_pto",
        "deng_entropy_level"
      ), envir = ns)
      results <- parallel::parLapply(cl, site_names, compute_site)
    } else {
      results <- parallel::mclapply(site_names, compute_site,
                                     mc.cores = n_cores_use)
    }
  } else {
    results <- lapply(site_names, compute_site)
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
