# Radar (Spider) Chart for Multi-Community Index Comparison

Creates a radar chart comparing diversity indices across multiple
communities. Each axis represents a different index, and each community
is drawn as a colored polygon. Values are normalized to 0-1 scale so
that indices with different magnitudes can be compared visually.

## Usage

``` r
plot_radar(communities, tax_tree, indices = NULL, title = NULL)
```

## Arguments

- communities:

  A named list of community vectors (named numeric).

- tax_tree:

  A data frame representing the taxonomic hierarchy, as produced by
  [`build_tax_tree`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/build_tax_tree.md).

- indices:

  Character vector specifying which indices to include. Default is all
  10 indices. Available: `"Shannon"`, `"Simpson"`, `"Delta"`,
  `"Delta_star"`, `"AvTD"`, `"VarTD"`, `"uTO"`, `"TO"`, `"uTO_plus"`,
  `"TO_plus"`.

- title:

  Optional character string for the plot title.

## Value

A `ggplot` object.

## Details

Each index value is normalized using min-max scaling across the
communities being compared:

\$\$x\_{norm} = \frac{x - x\_{min}}{x\_{max} - x\_{min}}\$\$

If all communities have the same value for an index (i.e., \\x\_{max} =
x\_{min}\\), the normalized value is set to 0.5.

The radar chart is built using polar coordinates in ggplot2. Each
community appears as a colored polygon overlay, making it easy to spot
which community scores higher on which indices.

## See also

[`compare_indices`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/compare_indices.md)
for tabular comparison

## Examples

``` r
# \donttest{
tax <- build_tax_tree(
  species = c("sp1", "sp2", "sp3", "sp4"),
  Genus   = c("G1", "G1", "G2", "G2"),
  Family  = c("F1", "F1", "F1", "F2")
)
comms <- list(
  Site_A = c(sp1 = 10, sp2 = 20, sp3 = 15, sp4 = 5),
  Site_B = c(sp1 = 5, sp2 = 5, sp3 = 5, sp4 = 5)
)
plot_radar(comms, tax)
#> Warning: Removed 2 rows containing missing values or values outside the scale range
#> (`geom_point()`).

# }
```
