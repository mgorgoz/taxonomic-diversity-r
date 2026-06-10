## Submission type

This is an update of an existing CRAN package (taxdiv 0.1.0 -> 1.0.0).

Main changes in this release (see NEWS.md for the full list):

* New interactive Shiny dashboard, `taxdiv_explorer()`.
* The Westhoff-Maarel 1-9 abundance scale is recommended but no longer
  enforced; any positive numeric abundances are accepted.
* The minimum for `n_iter` in the stochastic pTO pipeline was lowered from
  101 to 1 (101 was an inherited spreadsheet convention, not a requirement).
* Documentation and README refresh.

## R CMD check results

0 errors | 0 warnings | 1 note

* The note is the standard CRAN incoming feasibility note (maintainer
  identification).

## Test environments

* Local: macOS (aarch64), R 4.5.x
* GitHub Actions: ubuntu-latest, macOS-latest, windows-latest (R release)
* Win-builder: R release and R-devel (run before submission)

## Notes

* Words flagged as possibly misspelled in DESCRIPTION (AvTD, VarTD) are
  standard abbreviations for Average Taxonomic Distinctness and Variation
  in Taxonomic Distinctness (Clarke & Warwick, 1998, 2001).

* A few DOI links in the documentation (e.g. Wiley, Inter-Research MEPS)
  return 401/403 to automated requests but resolve correctly in a browser;
  the DOIs are valid.

## Package description

taxdiv computes taxonomic diversity indices for ecological community data.
It implements the Ozkan (2018) Deng entropy-based method alongside classical
Shannon/Simpson indices and Clarke & Warwick's taxonomic distinctness family.
The package provides 27 exported functions, 13 S3 methods, 7 plot types,
7 vignettes (6 English + 1 Turkish), and 3 example datasets. All 668 unit
tests pass on all tested platforms.

## Downstream dependencies

There are currently no downstream dependencies for this package.
