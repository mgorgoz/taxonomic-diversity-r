# Stochastic Resampling of Ozkan's pTO Index (Islem 2 / Run 2)

Implements the stochastic resampling procedure from Ozkan's Excel macro
(Islem 2). First performs a jackknife (Islem 1) to identify "happy"
(contributing) and "unhappy" (non-contributing) species, then runs
stochastic resampling where unhappy species are always included and
happy species are randomly included with 50\\

## Usage

``` r
ozkan_pto_resample(community, tax_tree, n_iter = 101L, seed = NULL)
```

## Arguments

- community:

  A named numeric vector of species abundances. Names must match the
  first column of `tax_tree`.

- tax_tree:

  A data frame with taxonomic hierarchy. First column is species names,
  subsequent columns are taxonomic ranks.

- n_iter:

  Number of stochastic iterations to run (default: 101). Must be \>=
  101.

- seed:

  Optional random seed for reproducibility (default: NULL).

## Value

A named list with components:

- uTO_plus_max:

  Maximum unweighted taxonomic distance across iterations

- TO_plus_max:

  Maximum weighted taxonomic distance across iterations

- uTO_max:

  Maximum unweighted taxonomic diversity across iterations

- TO_max:

  Maximum weighted taxonomic diversity across iterations

- uTO_plus_det:

  Deterministic uTO+ (first iteration, same as
  [`ozkan_pto()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/ozkan_pto.md))

- TO_plus_det:

  Deterministic TO+ (first iteration)

- uTO_det:

  Deterministic uTO (first iteration)

- TO_det:

  Deterministic TO (first iteration)

- n_iter:

  Number of iterations performed

- species_status:

  Named logical vector from jackknife (`TRUE` = happy)

- jackknife:

  Full jackknife result from
  [`ozkan_pto_jackknife()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/ozkan_pto_jackknife.md)

- n_positive:

  Number of iterations with positive uTO+

- iteration_results:

  Data frame with all iteration results

## Details

The algorithm follows the Excel macro's Islem 1 + Islem 2 logic:

1.  Run jackknife
    ([`ozkan_pto_jackknife()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/ozkan_pto_jackknife.md))
    to classify each species as happy or unhappy.

2.  Iteration 1: Use the original community (deterministic).

3.  Iterations 2..n_iter: For each species:

    - Unhappy species (AA = 0): always included with original abundance.

    - Happy species (AA \> 0): randomly included (50\\ or excluded. Uses
      `RANDBETWEEN(0,1) * abundance`.

4.  Return the maximum of each component across all iterations.

## References

Ozkan, K. (2018). A new proposed measure for estimating taxonomic
diversity. Turkish Journal of Forestry, 19(4), 336-346. DOI:
10.18182/tjf.441061

Ozkan, K. & Mert, A. (2022). Comparisons of Deng entropy-based taxonomic
diversity measures with the other diversity measures and introduction to
the new proposed (reinforced) estimators. FORESTIST, 72(2). DOI:
10.5152/forestist.2021.21025

## See also

[`ozkan_pto_jackknife()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/ozkan_pto_jackknife.md)
for the jackknife step,
[`ozkan_pto_sensitivity()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/ozkan_pto_sensitivity.md)
for Run 3,
[`ozkan_pto_full()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/ozkan_pto_full.md)
for the full pipeline.

## Examples

``` r
comm <- c(sp1 = 4, sp2 = 2, sp3 = 3, sp4 = 1, sp5 = 2)
tax <- data.frame(
  Species = paste0("sp", 1:5),
  Genus   = c("G1", "G1", "G1", "G2", "G2"),
  Family  = c("F1", "F1", "F1", "F1", "F1"),
  stringsAsFactors = FALSE
)
set.seed(42)
result <- ozkan_pto_resample(comm, tax, n_iter = 101)
result$species_status  # Happy/unhappy classification
#>  sp1  sp2  sp3  sp4  sp5 
#> TRUE TRUE TRUE TRUE TRUE 
```
