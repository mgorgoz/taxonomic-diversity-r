# Full Ozkan pTO Pipeline (Islem 1 + 2 + 3)

Runs the complete Ozkan taxonomic diversity analysis pipeline: jackknife
(Islem 1), stochastic resampling (Islem 2), and sensitivity analysis
(Islem 3), returning the maximum values across all three runs. This is
equivalent to running all three steps in the Excel macro sequentially.

## Usage

``` r
ozkan_pto_full(community, tax_tree, n_iter = 101L, seed = NULL)
```

## Arguments

- community:

  A named numeric vector of species abundances. Names must match the
  first column of `tax_tree`.

- tax_tree:

  A data frame with taxonomic hierarchy. First column is species names,
  subsequent columns are taxonomic ranks.

- n_iter:

  Number of stochastic iterations for Run 2 and Run 3 (default: 101,
  minimum: 101).

- seed:

  Optional random seed for reproducibility. If provided, Run 2 uses this
  seed and Run 3 uses `seed + 1` to ensure independent randomness.

## Value

A named list with components:

- uTO_plus:

  Final maximum uTO+ across all 3 runs

- TO_plus:

  Final maximum TO+ across all 3 runs

- uTO:

  Final maximum uTO across all 3 runs

- TO:

  Final maximum TO across all 3 runs

- run1:

  Deterministic pTO result (from
  [`ozkan_pto()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/ozkan_pto.md))

- run2:

  Full Run 2 result (from
  [`ozkan_pto_resample()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/ozkan_pto_resample.md))

- run3:

  Full Run 3 result (from
  [`ozkan_pto_sensitivity()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/ozkan_pto_sensitivity.md))

- jackknife:

  Jackknife result with species classifications

## Details

This function implements the full Excel macro pipeline in a single call:

1.  **Islem 1**: Leave-one-out jackknife to identify contributing
    (happy) vs non-contributing (unhappy) species, plus deterministic
    pTO calculation.

2.  **Islem 2**: Stochastic resampling – unhappy species are always
    included, happy species get 50\\

3.  **Islem 3**: Sensitivity analysis – unhappy species get \\P =
    (S-1)/S\\, happy species get a data-driven probability.

4.  **Final result**: Maximum values across all three runs.

## References

Ozkan, K. (2018). A new proposed measure for estimating taxonomic
diversity. Turkish Journal of Forestry, 19(4), 336-346. DOI:
10.18182/tjf.441061

Ozkan, K. & Mert, A. (2022). Comparisons of Deng entropy-based taxonomic
diversity measures with the other diversity measures and introduction to
the new proposed (reinforced) estimators. FORESTIST, 72(2). DOI:
10.5152/forestist.2021.21025

## See also

[`ozkan_pto()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/ozkan_pto.md)
for deterministic calculation only,
[`ozkan_pto_resample()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/ozkan_pto_resample.md)
for Run 2 only,
[`ozkan_pto_sensitivity()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/ozkan_pto_sensitivity.md)
for Run 3 only,
[`ozkan_pto_jackknife()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/ozkan_pto_jackknife.md)
for jackknife only.

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
result <- ozkan_pto_full(comm, tax, n_iter = 101)
result$uTO_plus  # Final maximum uTO+
#> [1] 3.060705
result$TO_plus   # Final maximum TO+
#> [1] 3.753852
```
