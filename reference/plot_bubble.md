# Bubble Chart of Species Contributions to Diversity

Creates a bubble chart showing each species' abundance (x-axis), average
taxonomic distance to other species (y-axis), and relative contribution
to the community (bubble size). Species that are both abundant and
taxonomically distant from others contribute most to overall taxonomic
diversity.

## Usage

``` r
plot_bubble(community, tax_tree, color_by = NULL, title = NULL)
```

## Arguments

- community:

  Named numeric vector of species abundances.

- tax_tree:

  A data frame representing the taxonomic hierarchy, as produced by
  [`build_tax_tree`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/build_tax_tree.md).

- color_by:

  Character string specifying which taxonomic rank to use for coloring
  bubbles. Must match a column name in `tax_tree`. If `NULL` (default),
  the highest available rank is used.

- title:

  Optional character string for the plot title.

## Value

A `ggplot` object.

## Details

For each species \\i\\, the average taxonomic distance is calculated as:

\$\$\bar{\omega}\_i = \frac{1}{S-1} \sum\_{j \neq i} \omega\_{ij}\$\$

where \\\omega\_{ij}\\ is the pairwise taxonomic distance and \\S\\ is
the number of species. Bubble size represents the product of relative
abundance and average distance, indicating each species' contribution to
overall taxonomic diversity.

## See also

[`tax_distance_matrix`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/tax_distance_matrix.md),
[`compare_indices`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/compare_indices.md)

## Examples

``` r
# \donttest{
comm <- c(sp1 = 25, sp2 = 18, sp3 = 30, sp4 = 12, sp5 = 8)
tax <- build_tax_tree(
  species = paste0("sp", 1:5),
  Genus   = c("G1", "G1", "G2", "G2", "G3"),
  Family  = c("F1", "F1", "F1", "F2", "F2"),
  Order   = c("O1", "O1", "O1", "O1", "O1")
)
plot_bubble(comm, tax)

# }
```
