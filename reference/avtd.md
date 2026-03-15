# Average Taxonomic Distinctness (Delta+)

Calculates the average taxonomic distinctness (AvTD, Delta+) based on
Clarke & Warwick (1998). This is a presence/absence-based measure of the
average taxonomic distance between all pairs of species.

## Usage

``` r
avtd(species, tax_tree, weights = NULL)
```

## Arguments

- species:

  Character vector of species names present in the community
  (presence-only data).

- tax_tree:

  A data frame representing the taxonomic hierarchy.

- weights:

  Optional numeric vector of weights for taxonomic levels.

## Value

A numeric value representing the average taxonomic distinctness
(Delta+).

## References

Clarke, K.R. & Warwick, R.M. (1998). A taxonomic distinctness index and
its statistical properties. Journal of Applied Ecology, 35, 523-531.

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
avtd(spp, tax)
#> [1] 2.666667
```
