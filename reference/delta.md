# Taxonomic Diversity Index (Delta)

Calculates the taxonomic diversity index (Delta) from Warwick & Clarke
(1995). This is the average weighted path length between every pair of
individuals, including same-species pairs (weighted 0).

## Usage

``` r
delta(community, tax_tree, weights = NULL)
```

## Arguments

- community:

  A named numeric vector of species abundances.

- tax_tree:

  A data frame with taxonomic hierarchy.

- weights:

  Optional numeric vector of path weights for each taxonomic level. If
  NULL, a linear scale is used (1, 2, 3, ...).

## Value

A numeric value representing taxonomic diversity (Delta).

## Details

\$\$\Delta = \frac{\sum\sum\_{i\<j} w\_{ij} x_i x_j + \sum_i 0 \cdot
x_i(x_i-1)/2}{\sum\sum\_{i\<j} x_i x_j + \sum_i x_i(x_i-1)/2}\$\$

## References

Warwick, R.M. & Clarke, K.R. (1995). New 'biodiversity' measures reveal
a decrease in taxonomic distinctness with increasing stress. Marine
Ecology Progress Series, 129, 301-305.

## See also

[`delta_star()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/delta_star.md)
for taxonomic distinctness (excluding same-species),
[`avtd()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/avtd.md)
for presence/absence-based AvTD,
[`ozkan_pto()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/ozkan_pto.md)
for Deng entropy-based alternative.

## Examples

``` r
comm <- c(sp1 = 5, sp2 = 3, sp3 = 3, sp4 = 1, sp5 = 3)
tax <- data.frame(
  Species = paste0("sp", 1:5),
  Genus = c("G1", "G1", "G2", "G2", "G2"),
  Family = c("F1", "F1", "F1", "F2", "F2"),
  stringsAsFactors = FALSE
)
delta(comm, tax)
#> [1] 1.352381
```
