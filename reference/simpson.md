# Simpson Diversity Index

Calculates the Simpson diversity index (1 - D) for a community.

## Usage

``` r
simpson(community, type = c("gini_simpson", "inverse", "dominance"))
```

## Arguments

- community:

  A numeric vector of species abundances.

- type:

  One of `"inverse"` (1/D), `"gini_simpson"` (1 - D), or `"dominance"`
  (D). Default is `"gini_simpson"`.

## Value

A numeric value representing the Simpson index.

## Details

Simpson's dominance index D is calculated as: \$\$D = \sum\_{i=1}^{S}
p_i^2\$\$ The Gini-Simpson index is \\1 - D\\ and the inverse Simpson is
\\1/D\\.

## See also

[`shannon()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/shannon.md)
for Shannon diversity.

## Examples

``` r
comm <- c(10, 5, 8, 3, 12)
simpson(comm)
#> [1] 0.7631579
simpson(comm, type = "inverse")
#> [1] 4.222222
```
