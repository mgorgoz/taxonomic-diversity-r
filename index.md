# taxdiv

## Overview

**taxdiv** is an R package for computing taxonomic diversity indices
from ecological community data. It provides a unified framework that
brings together classical diversity measures, Clarke & Warwick’s
taxonomic distinctness family, and the Deng entropy-based approach of
Ozkan (2018) — all in a single package.

Traditional diversity indices such as Shannon and Simpson treat all
species as equally distinct. However, a community of 10 species from 10
different families is taxonomically more diverse than 10 species from
the same genus. taxdiv accounts for this taxonomic structure using two
complementary frameworks:

- **Clarke & Warwick (1995, 1998, 2001)** — path-length-based taxonomic
  distinctness measures (Delta, Delta\*, AvTD, VarTD) with
  simulation-based significance testing
- **Ozkan (2018)** — Deng entropy-based taxonomic diversity (pTO) using
  Dempster-Shafer evidence theory and a slicing procedure, producing 8
  complementary indices

### Why taxdiv?

| Feature                                                       | vegan |   ape   | taxdiv  |
|---------------------------------------------------------------|:-----:|:-------:|:-------:|
| Shannon / Simpson                                             |  yes  |    –    |   yes   |
| Clarke & Warwick full suite (Delta, Delta\*, AvTD, VarTD)     |   –   | partial | **yes** |
| Ozkan pTO (8 indices, Run 1+2+3)                              |   –   |    –    | **yes** |
| Simulation-based significance testing (funnel plots)          |   –   |    –    | **yes** |
| Taxonomic rarefaction with bootstrap CI                       |   –   |    –    | **yes** |
| Stochastic resampling + sensitivity analysis                  |   –   |    –    | **yes** |
| Excel-to-results in one command                               |   –   |    –    | **yes** |
| Bias-corrected Shannon (Miller-Madow, Grassberger, Chao-Shen) |   –   |    –    | **yes** |

## Installation

``` r
# install.packages("devtools")
devtools::install_github("mgorgoz/taxonomic-diversity-r")
```

## Quick Start

### 1. From R vectors

``` r
library(taxdiv)

# Species abundances
community <- c(
  Quercus_robur       = 15,
  Pinus_nigra         = 8,
  Fagus_orientalis    = 12,
  Abies_nordmanniana  = 5,
  Juniperus_excelsa   = 3
)

# Taxonomic hierarchy
tax_tree <- build_tax_tree(
  species = names(community),
  Genus   = c("Quercus", "Pinus", "Fagus", "Abies", "Juniperus"),
  Family  = c("Fagaceae", "Pinaceae", "Fagaceae", "Pinaceae", "Cupressaceae"),
  Order   = c("Fagales", "Pinales", "Fagales", "Pinales", "Pinales")
)

# All 14 indices at once
compare_indices(community, tax_tree)

# Ozkan pTO — 8 values matching Excel macro output (Run 1+2+3)
ozkan_pto(community, tax_tree)
```

### 2. From Excel — one command

``` r
library(taxdiv)
library(readxl)

# Read your Excel file
data <- as.data.frame(read_excel("my_data.xlsx"))

# Compute all indices for all sites — automatic column detection
batch_analysis(data)
```

Output (16 columns):

      Site  N_Species Shannon Simpson  Delta Delta_star  AvTD  VarTD   uTO    TO uTO_plus TO_plus uTO_max TO_max uTO_plus_max TO_plus_max
      A1        6    1.494   0.757  1.622      2.138 2.333  0.667  2.14  3.49    3.891   5.244    2.14   3.49        3.891       5.244
      A2        5    1.577   0.784  1.719      2.243 2.500  0.500  1.98  3.21    3.456   4.872    1.98   3.21        3.456       4.872

A ready-to-use Excel template is included:
[`inst/templates/taxdiv_template.xlsx`](https://mgorgoz.github.io/taxonomic-diversity-r/inst/templates/taxdiv_template.xlsx)

## Features

### Diversity Indices (26 exported functions)

| Category             | Functions                                                                                                                                                                                                                                                                                                                                                                                                                                  | Description                                                                                                                                       |
|----------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------|
| **Classical**        | [`shannon()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/shannon.md), [`simpson()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/simpson.md)                                                                                                                                                                                                                                                                   | Shannon-Wiener H’ (with 3 bias corrections), Gini-Simpson                                                                                         |
| **Clarke & Warwick** | [`delta()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/delta.md), [`delta_star()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/delta_star.md), [`avtd()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/avtd.md), [`vartd()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/vartd.md)                                                                                                 | Taxonomic diversity (Delta), taxonomic distinctness (Delta\*), average taxonomic distinctness (AvTD), variation in taxonomic distinctness (VarTD) |
| **Ozkan pTO**        | [`ozkan_pto()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/ozkan_pto.md), [`pto_components()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/pto_components.md)                                                                                                                                                                                                                                                 | 8 Deng entropy-based indices: uTO, TO, uTO+, TO+ (all levels) and their max-informative-level variants                                            |
| **Ozkan Pipeline**   | [`ozkan_pto_resample()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/ozkan_pto_resample.md), [`ozkan_pto_sensitivity()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/ozkan_pto_sensitivity.md), [`ozkan_pto_jackknife()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/ozkan_pto_jackknife.md), [`ozkan_pto_full()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/ozkan_pto_full.md) | Stochastic resampling (Run 2), sensitivity analysis (Run 3), jackknife leave-one-out, full pipeline                                               |
| **Batch**            | [`batch_analysis()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/batch_analysis.md), [`compare_indices()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/compare_indices.md)                                                                                                                                                                                                                                     | Multi-site analysis from data frame, multi-community comparison                                                                                   |
| **Simulation**       | [`simulate_td()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/simulate_td.md)                                                                                                                                                                                                                                                                                                                                                | Random subsampling from species pool for significance testing                                                                                     |
| **Rarefaction**      | [`rarefaction_taxonomic()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/rarefaction_taxonomic.md)                                                                                                                                                                                                                                                                                                                            | Bootstrap rarefaction curves for 8 different indices                                                                                              |
| **Distance**         | [`tax_distance_matrix()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/tax_distance_matrix.md), [`build_tax_tree()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/build_tax_tree.md), [`deng_entropy_level()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/deng_entropy_level.md)                                                                                                                  | Taxonomic distance matrices, tree construction, per-level Deng entropy                                                                            |

### Visualization (7 plot types)

| Function                                                                                                    | Plot Type                                              |
|-------------------------------------------------------------------------------------------------------------|--------------------------------------------------------|
| [`plot_funnel()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/plot_funnel.md)                 | Funnel plot for AvTD/VarTD significance testing        |
| [`plot_rarefaction()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/plot_rarefaction.md)       | Rarefaction curves with bootstrap confidence intervals |
| [`plot_iteration()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/plot_iteration.md)           | Stochastic resampling iteration trajectories           |
| [`plot_radar()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/plot_radar.md)                   | Radar/spider chart for multi-community comparison      |
| [`plot_heatmap()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/plot_heatmap.md)               | Taxonomic similarity heatmap                           |
| [`plot_bubble()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/plot_bubble.md)                 | Bubble plot of community composition                   |
| [`plot_taxonomic_tree()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/plot_taxonomic_tree.md) | Dendrogram of taxonomic hierarchy                      |

### S3 Class System

All main output objects have dedicated
[`print()`](https://rdrr.io/r/base/print.html) and
[`summary()`](https://rdrr.io/r/base/summary.html) methods:

``` r
result <- batch_analysis(data)
print(result)    # Clean formatted output
summary(result)  # Min/max/mean/SD per index across sites

result <- ozkan_pto(community, tax_tree)
print(result)    # All 8 pTO values + Deng entropy by level
```

### Example Datasets

| Dataset           | Description                                          |
|-------------------|------------------------------------------------------|
| `anatolian_trees` | Anatolian tree species with full taxonomic hierarchy |
| `gazi_comm`       | Community abundance data from Gazi University campus |
| `gazi_gytk`       | Taxonomic classification for Gazi campus species     |

## Excel Macro Equivalence

taxdiv produces the same 8 Ozkan pTO values as the TD_OMD Excel macro:

    Excel Macro (TD_OMD)       R function output
    ────────────────────────────────────────────
    Run 1:  uT0+           ->  uTO_plus
            T0+            ->  TO_plus
    Run 2:  uT0            ->  uTO
            T0             ->  TO
    Run 3:  uT0+max        ->  uTO_plus_max
            T0+max         ->  TO_plus_max
            uT0max         ->  uTO_max
            T0max          ->  TO_max

The “max” variants use only informative taxonomic levels where Deng
entropy \> 0, matching the Excel macro’s Run 3 behavior.

## Theoretical Background

### The Problem: Why Species Counts Are Not Enough

Consider two forest plots, each containing 10 species. In the first
plot, all 10 species belong to the same genus. In the second, they span
5 families across 3 orders. Standard indices like Shannon and Simpson
would assign identical diversity scores to both — yet the second
community is clearly more diverse in an evolutionary and functional
sense. **Taxonomic diversity indices solve this problem by incorporating
the hierarchical relationships among species.**

### Ozkan pTO Method

Ozkan (2018) introduced a method that uses **Deng entropy** — a
generalization of Shannon entropy rooted in Dempster-Shafer evidence
theory — to measure how species are distributed across a taxonomic
hierarchy.

The method works in three stages:

**Run 1 — Deterministic calculation:** At each taxonomic level (genus,
family, order, etc.), Deng entropy measures how evenly species are
grouped. A level where all species fall into one group contributes zero
entropy (no diversity), while a level where species are spread evenly
across many groups contributes high entropy. The product of these
level-wise entropies gives **pTO** (taxonomic diversity). When
species-level Shannon entropy is also included in the product, the
result is **pTO+** (taxonomic distance), which captures both taxonomic
structure and within-community evenness.

**Run 2 — Stochastic resampling (slicing):** Species are removed one at
a time, starting with the least abundant. After each removal, all
indices are recalculated. This reveals each species’ contribution to
overall diversity: removing a “happy” species decreases diversity (it
was contributing positively), while removing an “unhappy” species
increases diversity (it was redundant in the taxonomic structure). The
maximum pTO value across all slicing steps represents the community’s
optimal taxonomic organization.

**Run 3 — Max-informative level variants:** Some taxonomic levels may
carry no information (e.g., when all species share the same order, Deng
entropy at that level is zero). Run 3 repeats the calculations using
only the levels where Deng entropy is greater than zero, producing the
`_max` variants of each index.

This three-stage pipeline yields **8 complementary indices**: `uTO`,
`TO`, `uTO_plus`, `TO_plus`, and their `_max` counterparts.

### Clarke & Warwick Taxonomic Distinctness

Clarke & Warwick (1995, 1998, 2001) proposed a family of indices based
on the **pairwise taxonomic path length** between species in a
classification tree.

- **Delta (Δ)** — average taxonomic distance between two randomly chosen
  individuals, weighted by abundance
- **Delta\* (Δ\*)** — same as Delta but excluding same-species pairs,
  isolating the pure taxonomic component
- **AvTD (Δ+)** — average taxonomic distinctness based on
  presence/absence only. Because it does not depend on abundance, AvTD
  is **independent of sample size**, making it ideal for comparing
  datasets collected with different sampling efforts
- **VarTD (Λ+)** — the variance in taxonomic path lengths. High VarTD
  indicates an uneven taxonomic tree (some branches are deep, others
  shallow)

To assess whether an observed AvTD or VarTD value is statistically
significant,
[`simulate_td()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/simulate_td.md)
draws random subsets of species from a master species pool and computes
expected distributions.
[`plot_funnel()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/plot_funnel.md)
then plots these as **95% confidence funnels** — if a community falls
below the funnel, its taxonomic diversity is significantly lower than
expected by chance.

## References

### Primary Methods

- Ozkan, K. (2018). A new proposed measure for estimating taxonomic
  diversity. *Turkish Journal of Forestry*, 19(4), 336-346. doi:
  [10.18182/tjf.441061](https://doi.org/10.18182/tjf.441061)

- Deng, Y. (2016). Deng entropy. *Chaos, Solitons & Fractals*, 91,
  549-553. doi:
  [10.1016/j.chaos.2016.08.011](https://doi.org/10.1016/j.chaos.2016.08.011)

### Taxonomic Distinctness

- Warwick, R.M. & Clarke, K.R. (1995). New ‘biodiversity’ measures
  reveal a decrease in taxonomic distinctness with increasing stress.
  *Marine Ecology Progress Series*, 129, 301-305. doi:
  [10.3354/meps129301](https://doi.org/10.3354/meps129301)

- Clarke, K.R. & Warwick, R.M. (1998). A taxonomic distinctness index
  and its statistical properties. *Journal of Applied Ecology*, 35(4),
  523-531. doi:
  [10.1046/j.1365-2664.1998.3540523.x](https://doi.org/10.1046/j.1365-2664.1998.3540523.x)

- Clarke, K.R. & Warwick, R.M. (1999). The taxonomic distinctness
  measure of biodiversity: weighting of step lengths between
  hierarchical levels. *Marine Ecology Progress Series*, 184, 21-29.
  doi: [10.3354/meps184021](https://doi.org/10.3354/meps184021)

- Clarke, K.R. & Warwick, R.M. (2001). A further biodiversity index
  applicable to species lists: variation in taxonomic distinctness.
  *Marine Ecology Progress Series*, 216, 265-278. doi:
  [10.3354/meps216265](https://doi.org/10.3354/meps216265)

### Classical Diversity

- Shannon, C.E. (1948). A mathematical theory of communication. *Bell
  System Technical Journal*, 27(3), 379-423. doi:
  [10.1002/j.1538-7305.1948.tb01338.x](https://doi.org/10.1002/j.1538-7305.1948.tb01338.x)

- Simpson, E.H. (1949). Measurement of diversity. *Nature*, 163, 688.
  doi: [10.1038/163688a0](https://doi.org/10.1038/163688a0)

### Evidence Theory

- Dempster, A.P. (1967). Upper and lower probabilities induced by a
  multivalued mapping. *The Annals of Mathematical Statistics*, 38(2),
  325-339. doi:
  [10.1214/aoms/1177698950](https://doi.org/10.1214/aoms/1177698950)

- Shafer, G. (1976). *A Mathematical Theory of Evidence*. Princeton
  University Press, Princeton, NJ.

### Bias Correction

- Chao, A. & Shen, T.-J. (2003). Nonparametric estimation of Shannon’s
  index of diversity when there are unseen species in sample.
  *Environmental and Ecological Statistics*, 10, 429-443. doi:
  [10.1023/A:1026096204727](https://doi.org/10.1023/A:1026096204727)

## Package Status

| Metric             | Value                         |
|--------------------|-------------------------------|
| R CMD check        | 0 errors, 0 warnings, 0 notes |
| Unit tests         | 610 passing                   |
| Exported functions | 26                            |
| S3 methods         | 13 (print, summary, plot)     |
| R source files     | 19                            |
| Test files         | 12                            |
| Example datasets   | 3                             |
| Vignettes          | 2 (English + Turkish)         |

## Roadmap

pkgdown documentation website

JOSS paper submission

CRAN submission

## Contributing

Contributions are welcome! If you encounter a bug or have an idea for a
new feature, please open an issue using our templates:

- [Bug
  Report](https://github.com/mgorgoz/taxonomic-diversity-r/issues/new?template=bug_report.md)
  — for errors, unexpected behavior, or incorrect results
- [Feature
  Request](https://github.com/mgorgoz/taxonomic-diversity-r/issues/new?template=feature_request.md)
  — for new indices, visualizations, or enhancements

For general questions about the package or taxonomic diversity methods,
feel free to start a [GitHub
Discussion](https://github.com/mgorgoz/taxonomic-diversity-r/discussions)
or open a blank issue.

## Citation

``` r
citation("taxdiv")
```

    Gorgoz MM, Ozkan K, Negiz MG (2026). taxdiv: Taxonomic Diversity Indices
    Using Deng Entropy. R package version 0.1.0.
    https://github.com/mgorgoz/taxonomic-diversity-r

    Ozkan K (2018). "A new proposed measure for estimating taxonomic diversity."
    Turkish Journal of Forestry, 19(4), 336-346. doi:10.18182/tjf.441061.

## License

MIT
