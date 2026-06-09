# Changelog

## taxdiv 1.0.0

### New Features

- [`batch_analysis()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/batch_analysis.md)
  gains an `indices` parameter for selective index computation. Three
  groups available: `"classical"` (Shannon, Simpson), `"clarke_warwick"`
  (Delta, Delta\*, AvTD, VarTD), and `"ozkan_pto"` (8 pTO indices).
  Default computes all. Unambiguous abbreviations accepted (e.g.,
  `"clas"`, `"clark"`, `"oz"`).

- **Westhoff-Maarel scale**: recommended for compatibility with Ozkan

  2018. and the original Excel macro, but no longer enforced as a hard
        requirement. The 1-9 limit in the original macro was a practical
        constraint of VBA (the macro froze on larger values), not a
        mathematical requirement of the algorithm. In R,
        [`ozkan_pto()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/ozkan_pto.md)
        and related functions accept any positive numeric abundance
        values.

- **`n_iter` minimum lowered from 101 to 1**:
  [`ozkan_pto_resample()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/ozkan_pto_resample.md),
  [`ozkan_pto_sensitivity()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/ozkan_pto_sensitivity.md),
  and
  [`ozkan_pto_full()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/ozkan_pto_full.md)
  now accept any positive integer for `n_iter`. The old `>= 101` floor
  was an inherited Excel-macro convention with no mathematical basis.
  With `n_iter = 1`, only the deterministic first iteration runs, so the
  result equals
  [`ozkan_pto()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/ozkan_pto.md)
  (Run 1). Note: because the procedure returns the maximum across
  iterations, the result is non-decreasing in `n_iter` — keep `n_iter`
  fixed when comparing sites and use a high value (e.g. 500-1000) for
  final estimates.

- **Natural sort** for site names: sites are now ordered as Site1,
  Site2, …, Site9, Site10 instead of the previous alphabetical order
  (Site1, Site10, Site11, Site2, …). Works with any naming pattern
  (A1-A10, Plot1-Plot20, etc.).

### Shiny Dashboard (taxdiv Explorer)

- **Index group selection**: three individual checkboxes with
  descriptions replace the previous “compute everything” approach. Each
  checkbox shows what the index group measures and any requirements.
- **Conditional settings panel**: n_iter and seed fields only appear
  when Ozkan pTO is selected, reducing UI clutter.
- **Adaptive bar chart and summary boxes**: when pTO is not selected,
  the chart and summary boxes automatically switch to Shannon or Delta
  instead of showing nothing.
- **Natural sort in Site Filter**: sidebar site checkboxes now follow
  natural sort order matching the results table.
- **Template file rebuilt**: 4-sheet Excel template (ENTER_DATA,
  EXAMPLE_text, EXAMPLE_numeric, INSTRUCTIONS) with Westhoff-Maarel
  abundances, styled headers (Calibri 12 bold / 11 normal), and English
  instructions with em dash formatting.
- **Download buttons fixed**: the Excel/CSV download buttons on the
  Results card now sit in a `card_footer()` so they remain visible
  regardless of how tall the results table grows (previously the table’s
  fill behaviour could push the buttons out of view).
- **`n_iter` default raised to 500** and its minimum lowered to 1,
  matching the relaxed function-level limit.

### Documentation & Internalization

- Translated all Turkish comments in `test-ozkan_pto_resample.R` to
  English.
- Roxygen documentation for
  [`batch_analysis()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/batch_analysis.md)
  updated with `indices` parameter details, Westhoff-Maarel
  requirements, and correct abbreviation examples.
- All vignette examples updated to Westhoff-Maarel scale abundances
  (`giris_rehberi.Rmd`, `introduction.Rmd`, `ozkan-pto.Rmd`,
  `visualization.Rmd`, `workflow.Rmd`).
- All Roxygen `@examples` updated to Westhoff-Maarel scale across 7 R
  source files.
- `inst/extdata/mediterranean_forest.csv` converted from raw counts to
  Westhoff-Maarel scale.
- Spelling wordlist updated (+11 terms: clarke, warwick, RStudio, etc.).

### Bug Fixes

- Fixed [`pmatch()`](https://rdrr.io/r/base/pmatch.html) ambiguity:
  documented that `"cl"` and `"cla"` are ambiguous between “classical”
  and “clarke_warwick”; minimum unambiguous prefixes are `"clas"` and
  `"clark"`.

### Testing

- 668 assertions across 12 test files (up from 610).
- New tests: index selection (9 tests), Westhoff-Maarel validation (4
  tests), natural sort (3 tests).
- All test data updated to Westhoff-Maarel scale across 5 test files
  (`test-batch_analysis.R`, `test-s3_methods.R`, `test-parallel.R`,
  `test-ozkan_pto.R`, `test-compare_indices.R`).
- R CMD check –as-cran: 0 errors, 0 warnings, 0 notes.

------------------------------------------------------------------------

## taxdiv 0.1.0

CRAN release: 2026-04-01

Initial release of taxdiv — taxonomic diversity indices using Deng
entropy.

### Core Functions

- Deng entropy at any taxonomic level
  ([`deng_entropy_level()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/deng_entropy_level.md))
- Ozkan (2018) pTO with slicing procedure
  ([`ozkan_pto()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/ozkan_pto.md),
  [`pto_components()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/pto_components.md))
  - 8 indices: uTO, TO, uTO+, TO+, uTO_max, TO_max, uTO+\_max, TO+\_max
  - `max_level` parameter for controlling taxonomic depth (NULL, “auto”,
    integer)
  - Presence-based entropy: equal weight (1/S) at each slice
- Ozkan & Mert (2022) reinforced estimators:
- Stochastic resampling — Run 2
  ([`ozkan_pto_resample()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/ozkan_pto_resample.md))
  - Random species inclusion/exclusion (50% probability per species)
  - Configurable iteration count and seed for reproducibility
- Sensitivity analysis — Run 3
  ([`ozkan_pto_sensitivity()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/ozkan_pto_sensitivity.md))
  - Species-specific inclusion probabilities derived from Run 2
  - Overall maximum across Run 1 + 2 + 3
- Jackknife leave-one-out analysis
  ([`ozkan_pto_jackknife()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/ozkan_pto_jackknife.md))
- Full pipeline in one call
  ([`ozkan_pto_full()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/ozkan_pto_full.md))
  - Combines Run 1 + Run 2 + Run 3 + jackknife
- Classical indices: Shannon H’
  ([`shannon()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/shannon.md))
  and Gini-Simpson
  ([`simpson()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/simpson.md))
  - Bias correction methods: Miller-Madow, Grassberger, Chao-Shen
- Clarke & Warwick taxonomic distinctness:
  - Taxonomic diversity Delta
    ([`delta()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/delta.md))
  - Taxonomic distinctness Delta\*
    ([`delta_star()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/delta_star.md))
  - Average taxonomic distinctness AvTD/Delta+
    ([`avtd()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/avtd.md))
  - Variation in taxonomic distinctness VarTD/Lambda+
    ([`vartd()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/vartd.md))
- Multi-community comparison
  ([`compare_indices()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/compare_indices.md))
  — 14 indices side by side
- Multi-site batch analysis
  ([`batch_analysis()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/batch_analysis.md))
  — automatic column detection
- Simulation-based significance testing
  ([`simulate_td()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/simulate_td.md))
  - Random subsampling from species pool for AvTD/VarTD confidence
    funnels
- Taxonomic rarefaction with bootstrap CI
  ([`rarefaction_taxonomic()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/rarefaction_taxonomic.md))
  - 8 index choices, configurable sample sizes and iterations
- Taxonomic distance matrix
  ([`tax_distance_matrix()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/tax_distance_matrix.md))
- Taxonomy tree builder
  ([`build_tax_tree()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/build_tax_tree.md))
- Parallel computing support for
  [`simulate_td()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/simulate_td.md),
  [`rarefaction_taxonomic()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/rarefaction_taxonomic.md),
  and
  [`batch_analysis()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/batch_analysis.md)
  via `parallel` and `n_cores` parameters

### S3 Class System

- 6 S3 classes: `compare_indices`, `batch_analysis`, `ozkan_pto`,
  `ozkan_pto_resample`, `ozkan_pto_sensitivity`,
  `rarefaction_taxonomic`, `td_simulation`
- 13 S3 methods: [`print()`](https://rdrr.io/r/base/print.html),
  [`summary()`](https://rdrr.io/r/base/summary.html), and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) for all main
  outputs

### Visualization

- Funnel plot for AvTD/VarTD significance testing
  ([`plot_funnel()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/plot_funnel.md))
- Rarefaction curves with bootstrap CI
  ([`plot_rarefaction()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/plot_rarefaction.md))
- Stochastic resampling iteration trajectories
  ([`plot_iteration()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/plot_iteration.md))
- Radar/spider chart for multi-community comparison
  ([`plot_radar()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/plot_radar.md))
- Taxonomic similarity heatmap
  ([`plot_heatmap()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/plot_heatmap.md))
- Bubble plot of community composition
  ([`plot_bubble()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/plot_bubble.md))
- Dendrogram of taxonomic hierarchy
  ([`plot_taxonomic_tree()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/plot_taxonomic_tree.md))

### Example Datasets

- `anatolian_trees` — Anatolian tree species with full taxonomic
  hierarchy
- `gazi_comm` — Community abundance data from Gazi University campus
- `gazi_gytk` — Taxonomic classification for Gazi campus species

### Documentation

- English vignette with 7 plots and worked examples
- Turkish (Turkce) vignette with 8 plots, formulas, and full walkthrough
- pkgdown website: <https://mgorgoz.github.io/taxonomic-diversity-r/>
- Pseudocode for pTO Run 1/2/3
  (`inst/pseudocode/taxonomic_diversity_pseudocode.md`)
- Pseudocode for Clarke & Warwick
  (`inst/pseudocode/taxonomic_distance_pseudocode.md`)

### Testing

- 610 assertions across 12 test files
- Integration tests covering full Run 1 -\> 2 -\> 3 pipeline
- Reproducibility tests with seed control
- Mathematical relationship checks (ordering constraints, boundary
  cases)
- R CMD check: 0 errors, 0 warnings, 0 notes

### Validation

- Verified against Ozkan (2018) hypothetical examples
- Verified against Ozkan & Mert (2022) reinforced estimator methodology
- Cross-validated with Kursad Ozkan’s Excel macro (TD_OMD.xlsm)
- Documented Excel macro bug: `Application.Calculate` missing inside
  tekerur2() and tekerur3() loops
