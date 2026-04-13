# Example Taxonomy: 8 Anatolian Tree Species

A data frame containing the taxonomic hierarchy for the 8 species in
[`gazi_comm`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/gazi_comm.md),
with 3 taxonomic ranks (Genus, Family, Order). This compact taxonomy
table is designed for quick demonstrations and unit testing.

## Usage

``` r
gazi_gytk
```

## Format

A data frame with 8 rows and 4 columns:

- Species:

  Binomial species name (character)

- Genus:

  Genus (character)

- Family:

  Family (character)

- Order:

  Order (character)

## Details

The taxonomy represents:

- 8 genera: Pinus, Cedrus, Quercus, Fagus, Juniperus, Carpinus

- 4 families: Pinaceae (3 spp), Fagaceae (3 spp), Cupressaceae (1),
  Betulaceae (1)

- 2 orders: Pinales (4 spp), Fagales (4 spp)

## References

Ozkan, K. (2018). A new proposed measure for estimating taxonomic
diversity. Turkish Journal of Forestry, 19(4), 336–346.

Ozkan, K. & Mert, A. (2022). Comparisons of Deng entropy-based taxonomic
diversity measures with the other diversity measures and introduction to
the new proposed (reinforced) estimators. FORESTIST, 72(2). DOI:
10.5152/forestist.2021.21025

## See also

[`gazi_comm`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/gazi_comm.md)
for the matching community vector,
[`build_tax_tree`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/build_tax_tree.md)
for building custom taxonomies.

## Examples

``` r
data(gazi_gytk)
gazi_gytk
#>             Species     Genus       Family   Order
#> 1       Pinus_nigra     Pinus     Pinaceae Pinales
#> 2      Pinus_brutia     Pinus     Pinaceae Pinales
#> 3     Cedrus_libani    Cedrus     Pinaceae Pinales
#> 4    Quercus_cerris   Quercus     Fagaceae Fagales
#> 5     Quercus_robur   Quercus     Fagaceae Fagales
#> 6  Fagus_orientalis     Fagus     Fagaceae Fagales
#> 7 Juniperus_excelsa Juniperus Cupressaceae Pinales
#> 8  Carpinus_betulus  Carpinus   Betulaceae Fagales

# Compute taxonomic distance matrix
tax_distance_matrix(gazi_gytk)
#>                   Pinus_nigra Pinus_brutia Cedrus_libani Quercus_cerris
#> Pinus_nigra                 0            1             2              3
#> Pinus_brutia                1            0             2              3
#> Cedrus_libani               2            2             0              3
#> Quercus_cerris              3            3             3              0
#> Quercus_robur               3            3             3              1
#> Fagus_orientalis            3            3             3              2
#> Juniperus_excelsa           3            3             3              3
#> Carpinus_betulus            3            3             3              3
#>                   Quercus_robur Fagus_orientalis Juniperus_excelsa
#> Pinus_nigra                   3                3                 3
#> Pinus_brutia                  3                3                 3
#> Cedrus_libani                 3                3                 3
#> Quercus_cerris                1                2                 3
#> Quercus_robur                 0                2                 3
#> Fagus_orientalis              2                0                 3
#> Juniperus_excelsa             3                3                 0
#> Carpinus_betulus              3                3                 3
#>                   Carpinus_betulus
#> Pinus_nigra                      3
#> Pinus_brutia                     3
#> Cedrus_libani                    3
#> Quercus_cerris                   3
#> Quercus_robur                    3
#> Fagus_orientalis                 3
#> Juniperus_excelsa                3
#> Carpinus_betulus                 0
```
