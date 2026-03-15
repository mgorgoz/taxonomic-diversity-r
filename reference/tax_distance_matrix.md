# Compute Taxonomic Distance Matrix

Calculates pairwise taxonomic distances between species based on their
positions in a taxonomic hierarchy. Distance is computed as the weighted
proportion of taxonomic levels at which two species diverge.

## Usage

``` r
tax_distance_matrix(tax_tree, species = NULL, weights = NULL)
```

## Arguments

- tax_tree:

  A data frame representing the taxonomic hierarchy. First column must
  be species names, subsequent columns are taxonomic ranks from lowest
  to highest.

- species:

  Optional character vector of species names to include. If NULL, all
  species in `tax_tree` are used.

- weights:

  Optional numeric vector of weights for each taxonomic level. If NULL,
  equal weights are assigned.

## Value

A symmetric matrix of taxonomic distances between species. With default
equal step weights (1, 2, 3, ...), values range from 0 (same species) to
the number of taxonomic levels (maximum distance when no common ancestor
is found at any level).

## Examples

``` r
tax <- data.frame(
  Species = c("Quercus_robur", "Pinus_nigra", "Fagus_orientalis"),
  Genus = c("Quercus", "Pinus", "Fagus"),
  Family = c("Fagaceae", "Pinaceae", "Fagaceae"),
  Order = c("Fagales", "Pinales", "Fagales"),
  stringsAsFactors = FALSE
)

tax_distance_matrix(tax)
#>                  Quercus_robur Pinus_nigra Fagus_orientalis
#> Quercus_robur                0           3                2
#> Pinus_nigra                  3           0                3
#> Fagus_orientalis             2           3                0
```
