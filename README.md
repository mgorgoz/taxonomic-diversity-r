# taxdiv <img src="man/figures/logo.png" align="right" height="139" alt="" />

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check](https://github.com/mgorgoz/taxonomic-diversity-r/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mgorgoz/taxonomic-diversity-r/actions/workflows/R-CMD-check.yaml)
[![CRAN status](https://www.r-pkg.org/badges/version/taxdiv)](https://CRAN.R-project.org/package=taxdiv)
<!-- badges: end -->

## Overview

**taxdiv** provides tools for calculating taxonomic diversity indices with a focus on Deng entropy — a novel information-theoretic measure that incorporates taxonomic hierarchy into diversity assessment.

Unlike traditional diversity indices (Shannon, Simpson) that only consider species abundances, Deng entropy uses the Dempster-Shafer evidence theory framework to account for taxonomic relationships among species. This makes it particularly useful for ecological studies where taxonomic structure carries meaningful information about community composition.

### Key Features

- **Deng Entropy**: Calculate diversity using belief function framework that incorporates taxonomic hierarchy
- **Classical Indices**: Shannon (H'), Simpson (D, 1-D, 1/D) diversity indices
- **Taxonomic Distinctness**: Average taxonomic distinctness (AvTD/Δ+) and variation in taxonomic distinctness (VarTD/Λ+) sensu Clarke & Warwick
- **Taxonomic Distance**: Compute pairwise taxonomic distance matrices
- **Taxonomy Tools**: Build and manipulate taxonomic trees for analysis

## Installation

You can install the development version of taxdiv from GitHub:

```r
# install.packages("devtools")
devtools::install_github("mgorgoz/taxonomic-diversity-r")
```

## Quick Start

```r
library(taxdiv)

# Define species abundances
community <- c(
  Quercus_robur       = 15,
  Pinus_nigra         = 8,
  Fagus_orientalis    = 12,
  Abies_nordmanniana  = 5,
  Juniperus_excelsa   = 3
)

# Build taxonomic tree
tax_tree <- build_tax_tree(
  species = names(community),
  Genus   = c("Quercus", "Pinus", "Fagus", "Abies", "Juniperus"),
  Family  = c("Fagaceae", "Pinaceae", "Fagaceae", "Pinaceae", "Cupressaceae"),
  Order   = c("Fagales", "Pinales", "Fagales", "Pinales", "Pinales")
)

# Classical indices
shannon(community)          # Shannon H'
simpson(community)          # Gini-Simpson

# Taxonomic diversity
deng_entropy(community, tax_tree)  # Deng entropy
avtd(names(community), tax_tree)   # Average taxonomic distinctness
vartd(names(community), tax_tree)  # Variation in taxonomic distinctness

# Taxonomic distance matrix
tax_distance_matrix(tax_tree)
```

## Background

This package stems from research on taxonomic diversity of perennial plant species in relation to habitat characteristics, utilizing the Deng entropy framework as a novel approach to quantifying biodiversity. The underlying methodology is described in:

> Deng, Y. (2016). Deng entropy. *Chaos, Solitons & Fractals*, 91, 549-553.

> Clarke, K.R. & Warwick, R.M. (1998). A taxonomic distinctness index and its statistical properties. *Journal of Applied Ecology*, 35, 523-531.

## Roadmap

- [ ] Funnel plots for AvTD/VarTD significance testing
- [ ] Permutation-based null models
- [ ] Visualization functions (diversity profiles, rarefaction curves)
- [ ] Vignettes with ecological case studies
- [ ] CRAN submission
- [ ] JOSS software paper

## Contributing

Contributions are welcome! Please open an issue to discuss proposed changes or submit a pull request.

## License

MIT © [Murat Görgöz](https://github.com/mgorgoz)

## Citation

If you use this package in your research, please cite:

```
Görgöz, M. (2026). taxdiv: Taxonomic Diversity Indices Using Deng Entropy.
R package version 0.1.0. https://github.com/mgorgoz/taxonomic-diversity-r
```
