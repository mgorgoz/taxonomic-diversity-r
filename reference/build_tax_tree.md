# Build a Taxonomic Tree from Species Data

Creates a taxonomic hierarchy data frame from species classification
information. This is a convenience function for constructing the
`tax_tree` input required by other functions in the package.

## Usage

``` r
build_tax_tree(species, ...)
```

## Arguments

- species:

  Character vector of species names.

- ...:

  Named character vectors for each taxonomic rank, in order from lowest
  to highest (e.g., Genus, Family, Order).

## Value

A data frame with species as the first column and taxonomic ranks as
subsequent columns.

## Examples

``` r
tree <- build_tax_tree(
  species = c("Quercus_robur", "Pinus_nigra", "Fagus_orientalis"),
  Genus   = c("Quercus", "Pinus", "Fagus"),
  Family  = c("Fagaceae", "Pinaceae", "Fagaceae"),
  Order   = c("Fagales", "Pinales", "Fagales")
)
```
