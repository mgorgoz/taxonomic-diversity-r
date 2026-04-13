# Calculate All Eight pTO Components (Convenience Wrapper)

Returns a named numeric vector with all eight Ozkan (2018; 2022)
components: four using all taxonomic levels and four using only the
informative levels (max version), matching the Excel macro's Run 1+2+3
output.

## Usage

``` r
pto_components(community, tax_tree)
```

## Arguments

- community:

  A named numeric vector of species abundances. Names must match the
  first column of `tax_tree`.

- tax_tree:

  A data frame with taxonomic hierarchy. First column is species names,
  subsequent columns are taxonomic ranks from lowest to highest (e.g.,
  Species, Genus, Family, Order, Class, Phylum, Kingdom).

## Value

A named numeric vector with eight elements: `uTO`, `TO`, `uTO_plus`,
`TO_plus`, `uTO_max`, `TO_max`, `uTO_plus_max`, `TO_plus_max`.

## See also

[`ozkan_pto()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/ozkan_pto.md)
for the full result including per-level entropy.

## Examples

``` r
comm <- c(sp1 = 4, sp2 = 2, sp3 = 3, sp4 = 1, sp5 = 2)
tax <- data.frame(
  Species = paste0("sp", 1:5),
  Genus   = c("G1", "G1", "G1", "G2", "G2"),
  Family  = c("F1", "F1", "F1", "F1", "F1"),
  stringsAsFactors = FALSE
)
pto_components(comm, tax)
#>          uTO           TO     uTO_plus      TO_plus      uTO_max       TO_max 
#>     2.632978     3.300651     3.060705     3.753852     2.632978     3.300651 
#> uTO_plus_max  TO_plus_max 
#>     3.060705     3.753852 
```
