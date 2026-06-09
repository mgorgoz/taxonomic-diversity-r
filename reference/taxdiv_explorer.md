# Launch the taxdiv Explorer Shiny Application

Opens an interactive Shiny dashboard for taxonomic diversity analysis.
The application allows users to upload species data (Excel or CSV),
preview the data, run
[`batch_analysis()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/batch_analysis.md)
with a single button click, and download the results as Excel or CSV.

## Usage

``` r
taxdiv_explorer(launch.browser = TRUE)
```

## Arguments

- launch.browser:

  Logical. If `TRUE` (default), opens the app in the default web
  browser. Set to `FALSE` to use the RStudio Viewer pane.

## Value

Does not return a value. Launches the Shiny application.

## Details

The following packages must be installed to use this function: `shiny`,
`bslib`, `DT`, `openxlsx`, `readxl`. These are listed under `Suggests`
in the package `DESCRIPTION` and are not installed automatically with
`taxdiv`. Install them with:

    install.packages(c("shiny", "bslib", "DT", "openxlsx", "readxl"))

## See also

[`batch_analysis()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/batch_analysis.md)
for the underlying analysis function.

## Examples

``` r
if (FALSE) { # \dontrun{
taxdiv_explorer()
} # }
```
