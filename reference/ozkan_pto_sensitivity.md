# Sensitivity Analysis of Ozkan's pTO Index (Islem 3 / Run 3)

Implements the sensitivity analysis procedure from Ozkan's Excel macro
(Islem 3). Uses the jackknife results from Run 2 to apply
species-specific inclusion probabilities: unhappy species get \\P =
(S-1)/S\\, happy species get a data-driven probability derived from Run
2 iteration results.

## Usage

``` r
ozkan_pto_sensitivity(
  community,
  tax_tree,
  run2_result,
  n_iter = NULL,
  seed = NULL
)
```

## Arguments

- community:

  A named numeric vector of species abundances.

- tax_tree:

  A data frame with taxonomic hierarchy.

- run2_result:

  The result from
  [`ozkan_pto_resample()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/ozkan_pto_resample.md).

- n_iter:

  Number of iterations (default: same as Run 2).

- seed:

  Optional random seed for reproducibility.

## Value

A named list with components:

- uTO_plus_max:

  Maximum uTO+ across Run 1, 2, and 3

- TO_plus_max:

  Maximum TO+ across all runs

- uTO_max:

  Maximum uTO across all runs

- TO_max:

  Maximum TO across all runs

- run3_uTO_plus_max:

  Maximum uTO+ from Run 3 only

- run3_TO_plus_max:

  Maximum TO+ from Run 3 only

- run3_uTO_max:

  Maximum uTO from Run 3 only

- run3_TO_max:

  Maximum TO from Run 3 only

- n_iter:

  Number of iterations performed

- species_probs:

  Named numeric vector of inclusion probabilities

- prob_happy:

  Probability used for happy species

- prob_unhappy:

  Probability used for unhappy species

- iteration_results:

  Data frame with all Run 3 iteration results

## Details

The algorithm follows the Excel macro's Islem 3 logic:

For each species, the inclusion probability depends on its jackknife
classification from Islem 1:

- **Unhappy species** (AA = 0, non-contributing): Included with
  probability \\(S-1)/S\\, where S is total species count. In the Excel
  formula: `IF(RANDBETWEEN(1, S) > 1, H2, 0)`.

- **Happy species** (AA \> 0, contributing): Included with probability
  derived from Run 2 results. In the Excel formula:
  `IF(L25 >= RANDBETWEEN(0, K22), H2, 0)`, where L25 is a summary score
  from Run 2 and K22 is the iteration count.

The happy species probability is computed as: \$\$P\_{happy} =
\frac{\max(0, N\_{positive} - S) + 1}{N\_{iter} + 1}\$\$

where \\N\_{positive}\\ is the number of Run 2 iterations that produced
a positive uTO+ value and S is the species count.

The maximum pTO across all three runs (Run 1, 2, 3) is the final result.

## References

Ozkan, K. (2018). A new proposed measure for estimating taxonomic
diversity. Turkish Journal of Forestry, 19(4), 336-346. DOI:
10.18182/tjf.441061

Ozkan, K. & Mert, A. (2022). Comparisons of Deng entropy-based taxonomic
diversity measures with the other diversity measures and introduction to
the new proposed (reinforced) estimators. FORESTIST, 72(2). DOI:
10.5152/forestist.2021.21025

## See also

[`ozkan_pto_resample()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/ozkan_pto_resample.md)
for Run 2,
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
run2 <- ozkan_pto_resample(comm, tax, n_iter = 101)
ozkan_pto_sensitivity(comm, tax, run2, n_iter = 101)
#> taxdiv -- Sensitivity Analysis (Run 3)
#> 
#>   Iterations: 101 
#>   P(happy):   0.7745 
#>   P(unhappy): 0.8 
#> 
#>   Overall maximum (across Run 1-2-3):
#>     uTO+ : 3.060705 
#>     TO+  : 3.753852 
#>     uTO  : 2.632978 
#>     TO   : 3.305854 
```
