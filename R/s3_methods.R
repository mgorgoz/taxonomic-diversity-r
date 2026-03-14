# ---- compare_indices S3 methods ----

#' @method print compare_indices
#' @export
print.compare_indices <- function(x, ...) {
  n_comm <- nrow(x)
  index_cols <- setdiff(names(x), "Community")

  cat("taxdiv -- Index Comparison\n")
  cat("  Communities:", n_comm, "\n")
  cat("  Indices:", paste(index_cols, collapse = ", "), "\n\n")

  print.data.frame(x, row.names = FALSE, ...)

  invisible(x)
}


#' @method summary compare_indices
#' @export
summary.compare_indices <- function(object, ...) {
  index_cols <- setdiff(names(object), "Community")
  n_comm <- nrow(object)

  cat("taxdiv -- Index Comparison Summary\n")
  cat("  Communities:", n_comm, "\n\n")

  if (n_comm == 1) {
    cat("  Single community -- no min/max comparison possible.\n")
    print.data.frame(object, row.names = FALSE)
  } else {
    stats_df <- data.frame(
      Index = index_cols,
      Min   = vapply(index_cols, function(i) min(object[[i]]), numeric(1)),
      Mean  = vapply(index_cols, function(i) mean(object[[i]]), numeric(1)),
      Max   = vapply(index_cols, function(i) max(object[[i]]), numeric(1)),
      row.names = NULL
    )
    stats_df$Min  <- round(stats_df$Min, 6)
    stats_df$Mean <- round(stats_df$Mean, 6)
    stats_df$Max  <- round(stats_df$Max, 6)
    print.data.frame(stats_df, row.names = FALSE)
  }

  invisible(object)
}


#' @method plot compare_indices
#' @export
plot.compare_indices <- function(x, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting.", call. = FALSE)
  }

  index_cols <- setdiff(names(x), "Community")
  long_df <- data.frame(
    Community = rep(x$Community, each = length(index_cols)),
    Index     = rep(index_cols, times = nrow(x)),
    Value     = as.numeric(t(as.matrix(x[, index_cols]))),
    stringsAsFactors = FALSE
  )
  long_df$Index <- factor(long_df$Index, levels = index_cols)

  p <- ggplot2::ggplot(long_df,
                       ggplot2::aes(x = .data$Index, y = .data$Value,
                                    fill = .data$Community)) +
    ggplot2::geom_col(position = "dodge", width = 0.7) +
    ggplot2::labs(title = "Diversity Indices Comparison",
                  x = "Index", y = "Value", fill = "Community") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 10),
      legend.position = "top"
    )

  print(p)
  invisible(p)
}


# ---- batch_analysis S3 methods ----

#' @method print batch_analysis
#' @export
print.batch_analysis <- function(x, ...) {
  n_sites <- nrow(x)
  index_cols <- setdiff(names(x), "Site")

  cat("taxdiv -- Batch Analysis\n")
  cat("  Sites:", n_sites, "\n")
  cat("  Indices:", length(index_cols), "\n\n")

  if (n_sites <= 20) {
    print.data.frame(x, row.names = FALSE, ...)
  } else {
    print.data.frame(utils::head(x, 10), row.names = FALSE, ...)
    cat("... (", n_sites, " sites total)\n")
  }

  invisible(x)
}


#' @method summary batch_analysis
#' @export
summary.batch_analysis <- function(object, ...) {
  index_cols <- setdiff(names(object), "Site")
  n_sites <- nrow(object)

  cat("taxdiv -- Batch Analysis Summary\n")
  cat("  Sites:", n_sites, "\n\n")

  if (n_sites == 1) {
    print.data.frame(object, row.names = FALSE)
  } else {
    stats_df <- data.frame(
      Index = index_cols,
      Min   = round(vapply(index_cols, function(i) min(object[[i]]), numeric(1)), 6),
      Mean  = round(vapply(index_cols, function(i) mean(object[[i]]), numeric(1)), 6),
      Max   = round(vapply(index_cols, function(i) max(object[[i]]), numeric(1)), 6),
      SD    = round(vapply(index_cols, function(i) stats::sd(object[[i]]), numeric(1)), 6),
      row.names = NULL
    )
    print.data.frame(stats_df, row.names = FALSE)
  }

  invisible(object)
}


# ---- ozkan_pto S3 methods ----

#' @method print ozkan_pto
#' @export
print.ozkan_pto <- function(x, ...) {
  cat("taxdiv -- Ozkan pTO Result\n\n")
  cat("  uTO  :", round(x$uTO, 6), "\n")
  cat("  TO   :", round(x$TO, 6), "\n")
  cat("  uTO+ :", round(x$uTO_plus, 6), "\n")
  cat("  TO+  :", round(x$TO_plus, 6), "\n")

  ed <- x$Ed_levels
  if (!is.null(ed) && any(ed != 0)) {
    cat("\n  Deng entropy by level:\n")
    for (i in seq_along(ed)) {
      cat("    ", names(ed)[i], ": ", round(ed[i], 6), "\n", sep = "")
    }
  }

  invisible(x)
}


# ---- ozkan_pto_resample S3 methods ----

#' @method print ozkan_pto_resample
#' @export
print.ozkan_pto_resample <- function(x, ...) {
  n_happy <- sum(x$species_status)
  n_unhappy <- sum(!x$species_status)

  cat("taxdiv -- Stochastic Resampling (Run 2)\n\n")
  cat("  Iterations:", x$n_iter, "\n")
  cat("  Species: ", n_happy, " happy, ", n_unhappy, " unhappy\n\n", sep = "")
  cat("  Maximum values:\n")
  cat("    uTO+ :", round(x$uTO_plus_max, 6), "\n")
  cat("    TO+  :", round(x$TO_plus_max, 6), "\n")
  cat("    uTO  :", round(x$uTO_max, 6), "\n")
  cat("    TO   :", round(x$TO_max, 6), "\n")

  invisible(x)
}


#' @method summary ozkan_pto_resample
#' @export
summary.ozkan_pto_resample <- function(object, ...) {
  n_happy <- sum(object$species_status)
  n_unhappy <- sum(!object$species_status)

  cat("taxdiv -- Stochastic Resampling Summary (Run 2)\n\n")
  cat("  Iterations:", object$n_iter, "\n")
  cat("  Positive iterations (uTO+ > 0):", object$n_positive, "\n\n")

  cat("  Deterministic vs Maximum:\n")
  cat("           Deterministic   Maximum\n")
  cat("    uTO+ : ", sprintf("%-14.6f  %.6f", object$uTO_plus_det, object$uTO_plus_max), "\n")
  cat("    TO+  : ", sprintf("%-14.6f  %.6f", object$TO_plus_det, object$TO_plus_max), "\n")
  cat("    uTO  : ", sprintf("%-14.6f  %.6f", object$uTO_det, object$uTO_max), "\n")
  cat("    TO   : ", sprintf("%-14.6f  %.6f", object$TO_det, object$TO_max), "\n")

  cat("\n  Species status:\n")
  status <- object$species_status
  for (sp in names(status)) {
    label <- if (status[sp]) "happy" else "unhappy"
    cat("    ", sp, ": ", label, "\n", sep = "")
  }

  invisible(object)
}


# ---- ozkan_pto_sensitivity S3 methods ----

#' @method print ozkan_pto_sensitivity
#' @export
print.ozkan_pto_sensitivity <- function(x, ...) {
  cat("taxdiv -- Sensitivity Analysis (Run 3)\n\n")
  cat("  Iterations:", x$n_iter, "\n")
  cat("  P(happy):  ", round(x$prob_happy, 4), "\n")
  cat("  P(unhappy):", round(x$prob_unhappy, 4), "\n\n")
  cat("  Overall maximum (across Run 1-2-3):\n")
  cat("    uTO+ :", round(x$uTO_plus_max, 6), "\n")
  cat("    TO+  :", round(x$TO_plus_max, 6), "\n")
  cat("    uTO  :", round(x$uTO_max, 6), "\n")
  cat("    TO   :", round(x$TO_max, 6), "\n")

  invisible(x)
}


#' @method summary ozkan_pto_sensitivity
#' @export
summary.ozkan_pto_sensitivity <- function(object, ...) {
  cat("taxdiv -- Sensitivity Analysis Summary (Run 3)\n\n")
  cat("  Iterations:", object$n_iter, "\n")
  cat("  P(happy):  ", round(object$prob_happy, 4), "\n")
  cat("  P(unhappy):", round(object$prob_unhappy, 4), "\n\n")

  cat("  Run 3 only vs Overall maximum:\n")
  cat("           Run 3 only     Overall (R1+R2+R3)\n")
  cat("    uTO+ : ", sprintf("%-14.6f  %.6f", object$run3_uTO_plus_max, object$uTO_plus_max), "\n")
  cat("    TO+  : ", sprintf("%-14.6f  %.6f", object$run3_TO_plus_max, object$TO_plus_max), "\n")
  cat("    uTO  : ", sprintf("%-14.6f  %.6f", object$run3_uTO_max, object$uTO_max), "\n")
  cat("    TO   : ", sprintf("%-14.6f  %.6f", object$run3_TO_max, object$TO_max), "\n")

  cat("\n  Species inclusion probabilities:\n")
  probs <- object$species_probs
  for (sp in names(probs)) {
    cat("    ", sp, ": ", round(probs[sp], 4), "\n", sep = "")
  }

  invisible(object)
}


# ---- rarefaction_taxonomic S3 methods ----

#' @method print rarefaction_taxonomic
#' @export
print.rarefaction_taxonomic <- function(x, ...) {
  idx <- attr(x, "index")
  total_n <- attr(x, "total_n")
  n_boot <- attr(x, "n_boot")
  ci_pct <- attr(x, "ci") * 100

  cat("taxdiv -- Rarefaction Curve\n")
  cat("  Index:", idx, "\n")
  cat("  Total N:", total_n, "\n")
  cat("  Bootstrap:", n_boot, "replicates\n")
  cat("  CI:", ci_pct, "%\n")
  cat("  Steps:", nrow(x), "\n\n")

  print.data.frame(utils::head(x, 10), row.names = FALSE, ...)
  if (nrow(x) > 10) {
    cat("... (", nrow(x), " rows total)\n")
  }

  invisible(x)
}


#' @method summary rarefaction_taxonomic
#' @export
summary.rarefaction_taxonomic <- function(object, ...) {
  idx <- attr(object, "index")
  total_n <- attr(object, "total_n")

  cat("taxdiv -- Rarefaction Summary\n")
  cat("  Index:", idx, "\n")
  cat("  Total N:", total_n, "\n")
  cat("  Steps:", nrow(object), "\n\n")

  cat("  Mean range: ", round(min(object$mean), 6), " -- ",
      round(max(object$mean), 6), "\n", sep = "")
  cat("  Final mean (N=", max(object$sample_size), "): ",
      round(object$mean[nrow(object)], 6), "\n", sep = "")

  # Plateau detection: where does the curve reach 95% of max?
  max_mean <- max(object$mean)
  if (max_mean > 0) {
    threshold <- 0.95 * max_mean
    plateau_idx <- which(object$mean >= threshold)[1]
    if (!is.na(plateau_idx)) {
      cat("  95% of max reached at N=", object$sample_size[plateau_idx], "\n",
          sep = "")
    }
  }

  invisible(object)
}
