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
comm <- c(sp1 = 9, sp2 = 5, sp3 = 8, sp4 = 2, sp5 = 3)
tax <- data.frame(
  Species = paste0("sp", 1:5),
  Genus   = c("G1", "G1", "G2", "G2", "G3"),
  Family  = c("F1", "F1", "F1", "F2", "F2"),
  stringsAsFactors = FALSE
)
rarefaction_taxonomic(comm, tax, index = "shannon", n_boot = 50)
#> taxdiv -- Rarefaction Curve
#>   Index: shannon 
#>   Total N: 27 
#>   Bootstrap: 50 replicates
#>   CI: 95 %
#>   Steps: 20 
#> 
#>  sample_size      mean     lower     upper        sd
#>            2 0.5129289 0.0000000 0.6931472 0.3071248
#>            3 0.7611903 0.0000000 1.0986123 0.3162456
#>            5 1.0921604 0.6730117 1.3321790 0.2429445
#>            6 1.1237115 0.8675632 1.3296613 0.1664130
#>            7 1.2150355 0.7084239 1.5498260 0.1957262
#>            9 1.2922652 0.9368883 1.5680125 0.1869525
#>           10 1.3065742 0.9081613 1.5571131 0.1991835
#>           11 1.3344961 1.0431492 1.5396483 0.1326990
#>           13 1.3883030 1.1417644 1.5559010 0.1226039
#>           14 1.4120225 1.1764582 1.5740973 0.1141675
#> ... ( 20  rows total)
```
