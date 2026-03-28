# Simulate Expected AvTD/VarTD Under Random Sampling

Generates the null distribution of Average Taxonomic Distinctness (AvTD)
and/or Variation in Taxonomic Distinctness (VarTD) by randomly drawing
species subsets from a regional species pool. Used to construct funnel
plots for statistical testing (Clarke & Warwick 1998, 2001).

## Usage

``` r
simulate_td(
  tax_tree,
  s_range = NULL,
  n_sim = 999L,
  index = c("both", "avtd", "vartd"),
  weights = NULL,
  ci = 0.95,
  seed = NULL,
  parallel = FALSE,
  n_cores = NULL
)
```

## Arguments

- tax_tree:

  A data frame representing the full regional species pool taxonomy.
  First column is species names, subsequent columns are taxonomic ranks
  from lowest to highest.

- s_range:

  Integer vector of species richness values to simulate. Default `NULL`
  uses `2:S` where S is the total number of species in `tax_tree`.

- n_sim:

  Number of random draws per species richness value (default 999).

- index:

  Which index to simulate: `"avtd"`, `"vartd"`, or `"both"` (default).

- weights:

  Optional numeric vector of weights for taxonomic levels. Passed to
  [`avtd()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/avtd.md)
  and
  [`vartd()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/vartd.md).

- ci:

  Confidence interval width (default 0.95).

- seed:

  Optional random seed for reproducibility.

- parallel:

  Logical. If `TRUE`, use parallel processing to speed up simulations.
  Default `FALSE`.

- n_cores:

  Number of CPU cores to use when `parallel = TRUE`. Default `NULL` uses
  up to 2 cores (CRAN policy limit).

## Value

A data frame with class `"td_simulation"` containing columns:

- s:

  Species richness (number of species drawn)

- mean_avtd:

  Mean simulated AvTD (if index includes avtd)

- lower_avtd:

  Lower CI bound for AvTD

- upper_avtd:

  Upper CI bound for AvTD

- mean_vartd:

  Mean simulated VarTD (if index includes vartd)

- lower_vartd:

  Lower CI bound for VarTD

- upper_vartd:

  Upper CI bound for VarTD

Attributes: `ci`, `index`, `n_sim`, `pool_size`.

## Details

For each value of S in `s_range`, `n_sim` random subsets of S species
are drawn (without replacement) from the full species pool in
`tax_tree`. AvTD and/or VarTD are computed for each random subset. The
mean and percentile-based confidence limits are recorded.

The resulting object can be passed to
[`plot_funnel()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/plot_funnel.md)
to produce the classic Clarke & Warwick funnel plot.

## References

Clarke, K.R. & Warwick, R.M. (1998). A taxonomic distinctness index and
its statistical properties. Journal of Applied Ecology, 35, 523-531.

Clarke, K.R. & Warwick, R.M. (2001). A further biodiversity index
applicable to species lists: variation in taxonomic distinctness. Marine
Ecology Progress Series, 216, 265-278.

## See also

[`plot_funnel()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/plot_funnel.md)
for visualisation,
[`avtd()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/avtd.md)
and
[`vartd()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/vartd.md)
for the underlying calculations.

## Examples

``` r
tax <- data.frame(
  Species = paste0("sp", 1:10),
  Genus   = rep(c("G1", "G2", "G3", "G4", "G5"), each = 2),
  Family  = rep(c("F1", "F1", "F2", "F2", "F3"), each = 2),
  Order   = rep(c("O1", "O1", "O2", "O2", "O2"), each = 2),
  stringsAsFactors = FALSE
)
sim <- simulate_td(tax, n_sim = 99, seed = 42)
sim
#> Taxonomic Distinctness Simulation
#>   Index: both 
#>   Species pool: 10 species
#>   Simulations per S: 99 
#>   Confidence interval: 95 %
#>   S range: 2 - 10 
#> 
#>   s mean_avtd lower_avtd upper_avtd mean_vartd lower_vartd upper_vartd
#>   2  2.595960   1.000000   3.000000  0.0000000   0.0000000   0.0000000
#>   3  2.649832   1.966667   3.000000  0.3120090   0.0000000   0.8888889
#>   4  2.597643   2.333333   2.833333  0.4189113   0.1388889   0.5833333
#>   5  2.616162   2.290000   2.800000  0.4321212   0.1600000   0.6400000
#>   6  2.601347   2.400000   2.733333  0.4465993   0.3288889   0.5155556
#>   7  2.593074   2.428571   2.666667  0.4594242   0.4126984   0.5351474
#>   8  2.596681   2.571429   2.642857  0.4620181   0.4438776   0.5306122
#>   9  2.599888   2.555556   2.611111  0.4617471   0.4598765   0.4691358
#>  10  2.600000   2.600000   2.600000  0.4622222   0.4622222   0.4622222
```
