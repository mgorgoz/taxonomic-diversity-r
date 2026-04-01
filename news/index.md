# Changelog

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
- Cross-validated with Kursad Ozkan’s Excel macro (TD_OMD.xlsm)
- Documented Excel macro bug: `Application.Calculate` missing inside
  tekerur2() and tekerur3() loops
