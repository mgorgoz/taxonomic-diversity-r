#' Launch the taxdiv Explorer Shiny Application
#'
#' Opens an interactive Shiny dashboard for taxonomic diversity analysis.
#' The application allows users to upload species data (Excel or CSV),
#' preview the data, run `batch_analysis()` with a single button click,
#' and download the results as Excel or CSV.
#'
#' @param launch.browser Logical. If `TRUE` (default), opens the app in the
#'   default web browser. Set to `FALSE` to use the RStudio Viewer pane.
#'
#' @return Does not return a value. Launches the Shiny application.
#'
#' @details
#' The following packages must be installed to use this function:
#' `shiny`, `bslib`, `DT`, `openxlsx`, `readxl`.
#' These are listed under `Suggests` in the package `DESCRIPTION` and are
#' not installed automatically with `taxdiv`. Install them with:
#'
#' ```r
#' install.packages(c("shiny", "bslib", "DT", "openxlsx", "readxl"))
#' ```
#'
#' @seealso [batch_analysis()] for the underlying analysis function.
#'
#' @examples
#' \dontrun{
#' taxdiv_explorer()
#' }
#'
#' @export
taxdiv_explorer <- function(launch.browser = TRUE) {

  required <- c("shiny", "bslib", "DT", "openxlsx", "readxl")
  missing  <- required[!vapply(required, requireNamespace, quietly = TRUE,
                               FUN.VALUE = logical(1))]

  if (length(missing) > 0) {
    stop(
      "The following packages are required but not installed:\n",
      "  ", paste(missing, collapse = ", "), "\n\n",
      "Install them with:\n",
      "  install.packages(c(",
      paste0('"', missing, '"', collapse = ", "),
      "))",
      call. = FALSE
    )
  }

  app_dir <- system.file("shinyapp", package = "taxdiv")

  if (!nzchar(app_dir)) {
    stop("Shiny app files not found. Please reinstall taxdiv.", call. = FALSE)
  }

  shiny::runApp(app_dir, launch.browser = launch.browser)
}
