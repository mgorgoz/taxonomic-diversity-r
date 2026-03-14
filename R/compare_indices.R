#' Compare All Diversity Indices Side by Side
#'
#' Computes all available diversity indices for one or more communities and
#' returns them in a single data frame. Optionally produces a grouped bar
#' plot for visual comparison.
#'
#' The function calculates the following indices:
#' \itemize{
#'   \item \strong{Shannon}: Shannon-Wiener entropy (\code{\link{shannon}})
#'   \item \strong{Simpson}: Gini-Simpson index (\code{\link{simpson}})
#'   \item \strong{Delta}: Clarke & Warwick taxonomic diversity (\code{\link{delta}})
#'   \item \strong{Delta_star}: Clarke & Warwick taxonomic distinctness (\code{\link{delta_star}})
#'   \item \strong{AvTD}: Average taxonomic distinctness (\code{\link{avtd}})
#'   \item \strong{VarTD}: Variation in taxonomic distinctness (\code{\link{vartd}})
#'   \item \strong{uTO}: Unweighted taxonomic diversity (Ozkan pTO)
#'   \item \strong{TO}: Weighted taxonomic diversity (Ozkan pTO)
#'   \item \strong{uTO_plus}: Unweighted taxonomic distance (Ozkan pTO)
#'   \item \strong{TO_plus}: Weighted taxonomic distance (Ozkan pTO)
#' }
#'
#' @param communities A named list of community vectors (named numeric),
#'   or a single named numeric vector. When a single vector is provided,
#'   it is wrapped in a list with name "Community".
#' @param tax_tree A data frame representing the taxonomic hierarchy,
#'   as produced by \code{\link{build_tax_tree}}.
#' @param correction Bias correction for the Shannon index. One of
#'   `"none"` (default), `"miller_madow"`, `"grassberger"`, or
#'   `"chao_shen"`. Passed to [shannon()]. See [shannon()] for details.
#' @param plot Logical. If \code{TRUE} and \pkg{ggplot2} is available,
#'   returns a list with both the data frame and a ggplot object.
#'   Default is \code{FALSE}.
#'
#' @return If \code{plot = FALSE}, a data frame with communities as rows
#'   and indices as columns. If \code{plot = TRUE}, a list with two elements:
#'   \describe{
#'     \item{table}{The data frame of index values.}
#'     \item{plot}{A \code{ggplot} object showing a grouped bar chart.}
#'   }
#'
#' @examples
#' tax <- build_tax_tree(
#'   species = c("sp1", "sp2", "sp3", "sp4"),
#'   Genus   = c("G1", "G1", "G2", "G2"),
#'   Family  = c("F1", "F1", "F1", "F2"),
#'   Order   = c("O1", "O1", "O1", "O1")
#' )
#'
#' # Single community
#' comm <- c(sp1 = 10, sp2 = 20, sp3 = 15, sp4 = 5)
#' compare_indices(comm, tax)
#'
#' # Multiple communities
#' comm_list <- list(
#'   Site_A = c(sp1 = 10, sp2 = 20, sp3 = 15, sp4 = 5),
#'   Site_B = c(sp1 = 5, sp2 = 5, sp3 = 5, sp4 = 5)
#' )
#' compare_indices(comm_list, tax)
#'
#' @seealso \code{\link{shannon}}, \code{\link{simpson}},
#'   \code{\link{delta}}, \code{\link{delta_star}},
#'   \code{\link{avtd}}, \code{\link{vartd}},
#'   \code{\link{ozkan_pto}}, \code{\link{pto_components}}
#'
#' @importFrom rlang .data
#' @export
compare_indices <- function(communities, tax_tree,
                            correction = c("none", "miller_madow",
                                           "grassberger", "chao_shen"),
                            plot = FALSE) {
  correction <- match.arg(correction)

  # --- Input validation ---
  if (!is.data.frame(tax_tree)) {
    stop("'tax_tree' must be a data frame.", call. = FALSE)
  }

  # If single community vector, wrap in a list

  if (is.numeric(communities) && !is.null(names(communities))) {
    communities <- list(Community = communities)
  }

  if (!is.list(communities)) {
    stop("'communities' must be a named numeric vector or a named list of numeric vectors.",
         call. = FALSE)
  }

  if (is.null(names(communities))) {
    names(communities) <- paste0("Community_", seq_along(communities))
  }

  # --- Compute indices for each community ---
  species_in_tree <- as.character(tax_tree[[1]])

  results <- lapply(names(communities), function(nm) {
    comm <- communities[[nm]]

    if (!is.numeric(comm) || is.null(names(comm))) {
      stop("Community '", nm, "' must be a named numeric vector.", call. = FALSE)
    }

    # Species present in this community
    sp_present <- names(comm[comm > 0])

    # Classical indices
    H <- shannon(comm, correction = correction)
    D <- simpson(comm)

    # Clarke & Warwick (need species in tax_tree)
    Delta     <- delta(comm, tax_tree)
    Delta_s   <- delta_star(comm, tax_tree)
    AvTD_val  <- avtd(sp_present, tax_tree)
    VarTD_val <- vartd(sp_present, tax_tree)

    # Ozkan pTO
    pto <- pto_components(comm, tax_tree)

    data.frame(
      Community  = nm,
      Shannon    = round(H, 6),
      Simpson    = round(D, 6),
      Delta      = round(Delta, 6),
      Delta_star = round(Delta_s, 6),
      AvTD       = round(AvTD_val, 6),
      VarTD      = round(VarTD_val, 6),
      uTO        = round(pto["uTO"], 6),
      TO         = round(pto["TO"], 6),
      uTO_plus   = round(pto["uTO_plus"], 6),
      TO_plus    = round(pto["TO_plus"], 6),
      row.names  = NULL,
      stringsAsFactors = FALSE
    )
  })

  result_df <- do.call(rbind, results)
  class(result_df) <- c("compare_indices", "data.frame")

  # --- Return table only if plot = FALSE ---
  if (!plot) {
    return(result_df)
  }

  # --- Build bar plot ---
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("Package 'ggplot2' is required for plotting. Returning table only.",
            call. = FALSE)
    return(result_df)
  }

  # Reshape to long format for ggplot
  index_cols <- setdiff(names(result_df), "Community")
  long_df <- data.frame(
    Community = rep(result_df$Community, each = length(index_cols)),
    Index     = rep(index_cols, times = nrow(result_df)),
    Value     = as.numeric(t(as.matrix(result_df[, index_cols]))),
    stringsAsFactors = FALSE
  )

  # Preserve index order
  long_df$Index <- factor(long_df$Index, levels = index_cols)

  p <- ggplot2::ggplot(long_df,
                       ggplot2::aes(x = .data$Index, y = .data$Value,
                                    fill = .data$Community)) +
    ggplot2::geom_col(position = "dodge", width = 0.7) +
    ggplot2::labs(
      title = "Diversity Indices Comparison",
      x = "Index",
      y = "Value",
      fill = "Community"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 10),
      legend.position = "top"
    )

  return(list(table = result_df, plot = p))
}
