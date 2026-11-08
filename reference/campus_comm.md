# Example Community Vector: 8 Anatolian Tree Species

A named numeric vector of species abundances for a single forest
community with 8 Anatolian tree species. Abundance values follow the
Westhoff & van der Maarel (1973) scale (1–9). This vector mirrors the
hypothetical example in Ozkan (2018).

## Usage

``` r
campus_comm
```

## Format

A named numeric vector with 8 elements. Names are species binomials
(underscore-separated); values are integer abundances (1–4).

## Details

The species include 3 genera from Pinaceae, 2 from Fagaceae, 1 each from
Cupressaceae and Betulaceae, spanning 2 orders (Pinales, Fagales).

Pair with
[`campus_gytk`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/campus_gytk.md)
for analysis:


    ozkan_pto(campus_comm, campus_gytk)
    compare_indices(campus_comm, campus_gytk)

## References

Ozkan, K. (2018). A new proposed measure for estimating taxonomic
diversity. Turkish Journal of Forestry, 19(4), 336–346.

Ozkan, K. & Mert, A. (2022). Comparisons of Deng entropy-based taxonomic
diversity measures with the other diversity measures and introduction to
the new proposed (reinforced) estimators. FORESTIST, 72(2). DOI:
10.5152/forestist.2021.21025

## See also

[`campus_gytk`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/campus_gytk.md)
for the matching taxonomy,
[`anatolian_trees`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/anatolian_trees.md)
and
[`med_forest`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/med_forest.md)
for multi-site datasets.

## Examples

``` r
data(campus_comm)
data(campus_gytk)

# Ozkan pTO
ozkan_pto(campus_comm, campus_gytk)
#> taxdiv -- Ozkan pTO Result
#> 
#>   All levels:
#>     uTO  : 6.828164 
#>     TO   : 10.00622 
#>     uTO+ : 7.014319 
#>     TO+  : 10.19237 
#> 
#>   Max informative level (level 3):
#>     uTO_max  : 6.828164 
#>     TO_max   : 10.00622 
#>     uTO+_max : 7.014319 
#>     TO+_max  : 10.19237 
#> 
#>   Deng entropy by level:
#>     Species: 2.079442
#>     Genus: 2.282174
#>     Family: 2.714915
#>     Order: 3.401197

# All indices at once
compare_indices(campus_comm, campus_gytk)
#> taxdiv -- Index Comparison
#>   Communities: 1 
#>   Indices: Shannon, Simpson, Delta, Delta_star, AvTD, VarTD, uTO, TO, uTO_plus, TO_plus, uTO_max, TO_max, uTO_plus_max, TO_plus_max 
#> 
#>  Community N_Species  Shannon  Simpson    Delta Delta_star     AvTD    VarTD
#>  Community         8 2.013806 0.858726 2.444444   2.696774 2.714286 0.346939
#>       uTO       TO uTO_plus  TO_plus  uTO_max   TO_max uTO_plus_max TO_plus_max
#>  6.828164 10.00622 7.014319 10.19237 6.828164 10.00622     7.014319    10.19237
```
