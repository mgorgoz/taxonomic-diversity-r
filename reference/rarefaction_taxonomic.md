# Taxonomic Diversity Rarefaction

Computes rarefaction curves for taxonomic diversity indices by
subsampling individuals from the community at increasing sample sizes.
Uses bootstrap resampling to estimate expected diversity and confidence
intervals at each sample size.

## Usage

``` r
rarefaction_taxonomic(
  community,
  tax_tree,
  index = c("shannon", "simpson", "species", "uTO", "TO", "uTO_plus", "TO_plus", "avtd"),
  steps = 20,
  n_boot = 100,
  ci = 0.95,
  seed = NULL,
  correction = c("none", "miller_madow", "grassberger", "chao_shen"),
  parallel = FALSE,
  n_cores = NULL
)
```

## Arguments

- community:

  A named numeric vector of species abundances (integers). Names must
  match the first column of `tax_tree`.

- tax_tree:

  A data frame with taxonomic hierarchy. First column is species names,
  subsequent columns are taxonomic ranks.

- index:

  Which index to rarefy. One of `"shannon"`, `"simpson"`, `"uTO"`,
  `"TO"`, `"uTO_plus"`, `"TO_plus"`, `"avtd"`, `"species"` (default:
  `"shannon"`).

- steps:

  Number of sample-size steps along the curve (default: 20).

- n_boot:

  Number of bootstrap replicates per step (default: 100).

- ci:

  Confidence interval width (default: 0.95).

- seed:

  Optional random seed for reproducibility (default: NULL).

- correction:

  Bias correction for the Shannon index. One of `"none"` (default),
  `"miller_madow"`, `"grassberger"`, or `"chao_shen"`. Only used when
  `index = "shannon"`. Passed to
  [`shannon()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/shannon.md).
  See
  [`shannon()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/shannon.md)
  for details.

- parallel:

  Logical. If `TRUE`, use parallel processing to speed up bootstrap
  resampling across sample sizes. Default `FALSE`.

- n_cores:

  Number of CPU cores to use when `parallel = TRUE`. Default `NULL` uses
  up to 2 cores (CRAN policy limit).

## Value

A data frame with columns:

- sample_size:

  Number of individuals in the subsample

- mean:

  Mean index value across bootstrap replicates

- lower:

  Lower bound of the confidence interval

- upper:

  Upper bound of the confidence interval

- sd:

  Standard deviation across replicates

## Details

The algorithm works as follows:

1.  Expand the abundance vector into an individual-level vector (e.g.,
    c(sp1=3, sp2=2) becomes c("sp1","sp1","sp1","sp2","sp2")).

2.  For each sample size (from min to total N), draw `n_boot` random
    subsamples without replacement.

3.  For each subsample, count species abundances and compute the chosen
    diversity index.

4.  Return the mean and confidence interval at each step.

## References

Gotelli, N.J. & Colwell, R.K. (2001). Quantifying biodiversity:
procedures and pitfalls in the measurement and comparison of species
richness. Ecology Letters, 4, 379-391.

Ozkan, K. (2018). A new proposed measure for estimating taxonomic
diversity. Turkish Journal of Forestry, 19(4), 336-346.

Ozkan, K. & Mert, A. (2022). Comparisons of Deng entropy-based taxonomic
diversity measures with the other diversity measures and introduction to
the new proposed (reinforced) estimators. FORESTIST, 72(2). DOI:
10.5152/forestist.2021.21025

## See also

[`plot_rarefaction()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/plot_rarefaction.md)
for visualising the rarefaction curve,
[`ozkan_pto()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/ozkan_pto.md)
for full pTO calculation,
[`shannon()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/shannon.md)
and
[`simpson()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/simpson.md)
for classical indices,
[`avtd()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/avtd.md)
for average taxonomic distinctness.

## Examples

``` r
comm <- c(sp1 = 10, sp2 = 5, sp3 = 8, sp4 = 2, sp5 = 3)
tax <- data.frame(
  Species = paste0("sp", 1:5),
  Genus   = c("G1", "G1", "G2", "G2", "G3"),
  Family  = c("F1", "F1", "F1", "F2", "F2"),
  stringsAsFactors = FALSE
)
rarefaction_taxonomic(comm, tax, index = "shannon", n_boot = 50)
#> taxdiv -- Rarefaction Curve
#>   Index: shannon 
#>   Total N: 28 
#>   Bootstrap: 50 replicates
#>   CI: 95 %
#>   Steps: 20 
#> 
#>  sample_size      mean     lower     upper        sd
#>            2 0.4990660 0.0000000 0.6931472 0.3143820
#>            3 0.7519484 0.0000000 1.0986123 0.3129182
#>            5 1.0436130 0.5392395 1.5470547 0.2634342
#>            6 1.1578025 0.6365142 1.5607104 0.2258762
#>            7 1.1932895 0.8321740 1.5498260 0.1972654
#>            9 1.2095448 0.8486856 1.5229551 0.1958911
#>           10 1.2659255 0.8979457 1.5453400 0.1943499
#>           12 1.3698280 1.0905674 1.5607104 0.1331104
#>           13 1.3517378 1.1275231 1.5188263 0.1151802
#>           14 1.3991807 1.1635515 1.5457653 0.1019017
#> ... ( 20  rows total)
```
