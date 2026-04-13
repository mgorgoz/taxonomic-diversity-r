# taxdiv 0.1.0

Initial release of taxdiv — taxonomic diversity indices using Deng entropy.

## Core Functions

* Deng entropy at any taxonomic level (`deng_entropy_level()`)
* Ozkan (2018) pTO with slicing procedure (`ozkan_pto()`, `pto_components()`)
  - 8 indices: uTO, TO, uTO+, TO+, uTO_max, TO_max, uTO+_max, TO+_max
  - `max_level` parameter for controlling taxonomic depth (NULL, "auto", integer)
  - Presence-based entropy: equal weight (1/S) at each slice
* Ozkan & Mert (2022) reinforced estimators:
* Stochastic resampling — Run 2 (`ozkan_pto_resample()`)
  - Random species inclusion/exclusion (50% probability per species)
  - Configurable iteration count and seed for reproducibility
* Sensitivity analysis — Run 3 (`ozkan_pto_sensitivity()`)
  - Species-specific inclusion probabilities derived from Run 2
  - Overall maximum across Run 1 + 2 + 3
* Jackknife leave-one-out analysis (`ozkan_pto_jackknife()`)
* Full pipeline in one call (`ozkan_pto_full()`)
  - Combines Run 1 + Run 2 + Run 3 + jackknife
* Classical indices: Shannon H' (`shannon()`) and Gini-Simpson (`simpson()`)
  - Bias correction methods: Miller-Madow, Grassberger, Chao-Shen
* Clarke & Warwick taxonomic distinctness:
  - Taxonomic diversity Delta (`delta()`)
  - Taxonomic distinctness Delta* (`delta_star()`)
  - Average taxonomic distinctness AvTD/Delta+ (`avtd()`)
  - Variation in taxonomic distinctness VarTD/Lambda+ (`vartd()`)
* Multi-community comparison (`compare_indices()`) — 14 indices side by side
* Multi-site batch analysis (`batch_analysis()`) — automatic column detection
* Simulation-based significance testing (`simulate_td()`)
  - Random subsampling from species pool for AvTD/VarTD confidence funnels
* Taxonomic rarefaction with bootstrap CI (`rarefaction_taxonomic()`)
  - 8 index choices, configurable sample sizes and iterations
* Taxonomic distance matrix (`tax_distance_matrix()`)
* Taxonomy tree builder (`build_tax_tree()`)
* Parallel computing support for `simulate_td()`, `rarefaction_taxonomic()`,
  and `batch_analysis()` via `parallel` and `n_cores` parameters

## S3 Class System

* 6 S3 classes: `compare_indices`, `batch_analysis`, `ozkan_pto`,
  `ozkan_pto_resample`, `ozkan_pto_sensitivity`, `rarefaction_taxonomic`,
  `td_simulation`
* 13 S3 methods: `print()`, `summary()`, and `plot()` for all main outputs

## Visualization

* Funnel plot for AvTD/VarTD significance testing (`plot_funnel()`)
* Rarefaction curves with bootstrap CI (`plot_rarefaction()`)
* Stochastic resampling iteration trajectories (`plot_iteration()`)
* Radar/spider chart for multi-community comparison (`plot_radar()`)
* Taxonomic similarity heatmap (`plot_heatmap()`)
* Bubble plot of community composition (`plot_bubble()`)
* Dendrogram of taxonomic hierarchy (`plot_taxonomic_tree()`)

## Example Datasets

* `anatolian_trees` — Anatolian tree species with full taxonomic hierarchy
* `gazi_comm` — Community abundance data from Gazi University campus
* `gazi_gytk` — Taxonomic classification for Gazi campus species

## Documentation

* English vignette with 7 plots and worked examples
* Turkish (Turkce) vignette with 8 plots, formulas, and full walkthrough
* pkgdown website: https://mgorgoz.github.io/taxonomic-diversity-r/
* Pseudocode for pTO Run 1/2/3 (`inst/pseudocode/taxonomic_diversity_pseudocode.md`)
* Pseudocode for Clarke & Warwick (`inst/pseudocode/taxonomic_distance_pseudocode.md`)

## Testing

* 610 assertions across 12 test files
* Integration tests covering full Run 1 -> 2 -> 3 pipeline
* Reproducibility tests with seed control
* Mathematical relationship checks (ordering constraints, boundary cases)
* R CMD check: 0 errors, 0 warnings, 0 notes

## Validation

* Verified against Ozkan (2018) hypothetical examples
* Verified against Ozkan & Mert (2022) reinforced estimator methodology
* Cross-validated with Kursad Ozkan's Excel macro (TD_OMD.xlsm)
* Documented Excel macro bug: `Application.Calculate` missing inside
  tekerur2() and tekerur3() loops
