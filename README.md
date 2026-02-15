# taxdiv

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check](https://github.com/mgorgoz/taxonomic-diversity-r/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mgorgoz/taxonomic-diversity-r/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## Overview

**taxdiv** provides tools for calculating taxonomic diversity indices with a focus on Deng entropy-based measures from Ozkan (2018). Unlike traditional diversity indices (Shannon, Simpson) that only consider species abundances, Deng entropy uses the Dempster-Shafer evidence theory framework to account for taxonomic relationships among species.

### Key Features

- **Ozkan pTO**: Taxonomic diversity (pTO) and taxonomic distance (pTO+) using Deng entropy with slicing procedure
- **Deng Entropy**: Calculate diversity at each taxonomic level using belief function framework
- **Classical Indices**: Shannon (H'), Simpson (D, 1-D, 1/D) diversity indices
- **Taxonomic Distinctness**: Average taxonomic distinctness (AvTD/Delta+) and variation (VarTD/Lambda+) sensu Clarke & Warwick
- **Taxonomic Distance**: Pairwise taxonomic distance matrices

## Installation

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

# Ozkan (2018) Deng entropy-based taxonomic diversity
result <- ozkan_pto(community, tax_tree)
result$uTO_plus  # Unweighted taxonomic distance
result$TO_plus   # Weighted taxonomic distance
result$uTO       # Unweighted taxonomic diversity
result$TO        # Weighted taxonomic diversity

# Clarke & Warwick indices
avtd(names(community), tax_tree)   # Average taxonomic distinctness
vartd(names(community), tax_tree)  # Variation in taxonomic distinctness

# Taxonomic distance matrix
tax_distance_matrix(tax_tree)
```

## Background

This package implements the Deng entropy-based taxonomic diversity measure proposed by Ozkan (2018), which generalizes Shannon entropy through Dempster-Shafer evidence theory to incorporate taxonomic hierarchy information. The method uses a slicing procedure where species receive equal weight at each slice, and abundance information enters indirectly through which species survive progressive elimination steps.

### References

- Ozkan, K. (2018). A new proposed measure for estimating taxonomic diversity. *Turkish Journal of Forestry*, 19(4), 336-346.
- Deng, Y. (2016). Deng entropy. *Chaos, Solitons & Fractals*, 91, 549-553.
- Clarke, K.R. & Warwick, R.M. (1998). A taxonomic distinctness index and its statistical properties. *Journal of Applied Ecology*, 35, 523-531.

## Roadmap

- [x] Core Deng entropy calculation
- [x] Ozkan pTO formula with slicing procedure
- [x] Classical diversity indices (Shannon, Simpson)
- [x] Clarke & Warwick taxonomic distinctness (Delta, Delta*, AvTD, VarTD)
- [x] R CMD check CRAN compliance
- [x] GitHub Actions CI/CD
- [ ] Funnel plots for AvTD/VarTD significance testing
- [ ] Bias-corrected entropy estimators
- [ ] Vignettes with ecological case studies
- [ ] CRAN submission

## License

MIT
