# Taxonomic Distinctness (Delta\*)

Calculates the taxonomic distinctness (Delta\*) from Warwick & Clarke
(1995). This is the average weighted path length between individuals of
different species only.

## Usage

``` r
delta_star(community, tax_tree, weights = NULL)
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

A numeric value representing taxonomic distinctness (Delta\*).

## Details

\$\$\Delta^\* = \frac{\sum\sum\_{i\<j} w\_{ij} x_i x_j}
{\sum\sum\_{i\<j} x_i x_j}\$\$

## References

Warwick, R.M. & Clarke, K.R. (1995). New 'biodiversity' measures reveal
a decrease in taxonomic distinctness with increasing stress. Marine
Ecology Progress Series, 129, 301-305.

## See also

[`delta()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/delta.md)
for taxonomic diversity (including same-species),
[`avtd()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/avtd.md)
and
[`vartd()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/vartd.md)
for presence/absence measures.

## Examples

``` r
comm <- c(sp1 = 5, sp2 = 3, sp3 = 3, sp4 = 1, sp5 = 3)
tax <- data.frame(
  Species = paste0("sp", 1:5),
  Genus = c("G1", "G1", "G2", "G2", "G2"),
  Family = c("F1", "F1", "F1", "F2", "F2"),
  stringsAsFactors = FALSE
)
delta_star(comm, tax)
#> [1] 1.651163
```
