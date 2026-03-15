# Jackknife Analysis for Ozkan's pTO Index (Islem 1 / Run 1)

Implements the leave-one-out jackknife procedure from the Ozkan Excel
macro (Islem 1). Removes each species one at a time, recalculates pTO,
and identifies "happy" (contributing) and "unhappy" (non-contributing)
species. A species is "happy" if its removal decreases the pTO index,
indicating it positively contributes to the community's taxonomic
diversity.

## Usage

``` r
ozkan_pto_jackknife(community, tax_tree, component = "uTO_plus")
```

## Arguments

- community:

  A named numeric vector of species abundances. Names must match the
  first column of `tax_tree`.

- tax_tree:

  A data frame with taxonomic hierarchy. First column is species names,
  subsequent columns are taxonomic ranks.

- component:

  Character string specifying which pTO component to use for the
  happy/unhappy classification. One of `"uTO_plus"` (default),
  `"TO_plus"`, `"uTO"`, or `"TO"`.

## Value

A named list with components:

- full_result:

  The
  [`ozkan_pto()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/ozkan_pto.md)
  result for the full community

- jackknife_results:

  Data frame with leave-one-out results per species

- species_status:

  Named logical vector: `TRUE` = happy (contributing), `FALSE` = unhappy
  (non-contributing)

- n_happy:

  Number of happy species

- n_unhappy:

  Number of unhappy species

## Details

The jackknife procedure follows the Excel macro's Islem 1 logic:

1.  Compute pTO for the full community.

2.  For each species i, remove it and compute pTO for the remaining
    community (leave-one-out).

3.  Compare each leave-one-out result against the full-community value.

4.  If removing species i DECREASES the specified component (pTO becomes
    smaller), species i is classified as "happy" (contributing).

5.  If removing species i does NOT decrease the component, species i is
    classified as "unhappy" (non-contributing).

The happy/unhappy classification is used by
[`ozkan_pto_resample()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/ozkan_pto_resample.md)
(Islem 2) and
[`ozkan_pto_sensitivity()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/ozkan_pto_sensitivity.md)
(Islem 3) to apply different resampling probabilities to each species
category.

## References

Ozkan, K. (2018). A new proposed measure for estimating taxonomic
diversity. Turkish Journal of Forestry, 19(4), 336-346. DOI:
10.18182/tjf.441061

## See also

[`ozkan_pto()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/ozkan_pto.md)
for the core calculation,
[`ozkan_pto_resample()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/ozkan_pto_resample.md)
for Run 2,
[`ozkan_pto_full()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/ozkan_pto_full.md)
for the full 3-run pipeline.

## Examples

``` r
comm <- c(sp1 = 4, sp2 = 2, sp3 = 3, sp4 = 1, sp5 = 2)
tax <- data.frame(
  Species = paste0("sp", 1:5),
  Genus   = c("G1", "G1", "G1", "G2", "G2"),
  Family  = c("F1", "F1", "F1", "F1", "F1"),
  stringsAsFactors = FALSE
)
jk <- ozkan_pto_jackknife(comm, tax)
jk$species_status   # Which species are happy (contributing)?
#>  sp1  sp2  sp3  sp4  sp5 
#> TRUE TRUE TRUE TRUE TRUE 
jk$n_happy           # How many happy species?
#> [1] 5
```
