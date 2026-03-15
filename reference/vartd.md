# Variation in Taxonomic Distinctness (Lambda+)

Calculates the variation in taxonomic distinctness (VarTD, Lambda+)
based on Clarke & Warwick (2001).

## Usage

``` r
vartd(species, tax_tree, weights = NULL)
```

## Arguments

- species:

  Character vector of species names present in the community.

- tax_tree:

  A data frame representing the taxonomic hierarchy.

- weights:

  Optional numeric vector of weights for taxonomic levels.

## Value

A numeric value representing the variation in taxonomic distinctness
(Lambda+).

## References

Clarke, K.R. & Warwick, R.M. (2001). A further biodiversity index
applicable to species lists: variation in taxonomic distinctness. Marine
Ecology Progress Series, 216, 265-278.

## Examples

``` r
tax <- data.frame(
  Species = c("Quercus_robur", "Pinus_nigra", "Fagus_orientalis",
              "Abies_nordmanniana"),
  Genus = c("Quercus", "Pinus", "Fagus", "Abies"),
  Family = c("Fagaceae", "Pinaceae", "Fagaceae", "Pinaceae"),
  Order = c("Fagales", "Pinales", "Fagales", "Pinales"),
  stringsAsFactors = FALSE
)

spp <- c("Quercus_robur", "Pinus_nigra", "Fagus_orientalis")
vartd(spp, tax)
#> [1] 0.2222222
```
