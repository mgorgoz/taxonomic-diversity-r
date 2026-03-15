# Calculate Ozkan's Taxonomic Diversity Index (pTO)

Computes the four components of the Deng entropy-based taxonomic
diversity measure proposed by Ozkan (2018): weighted/unweighted
taxonomic diversity (TO, uTO) and weighted/unweighted taxonomic distance
(TO+, uTO+).

## Usage

``` r
ozkan_pto(community, tax_tree, max_level = NULL)
```

## Arguments

- community:

  A named numeric vector of species abundances. Names must match the
  first column of `tax_tree`.

- tax_tree:

  A data frame with taxonomic hierarchy. First column is species names,
  subsequent columns are taxonomic ranks from lowest to highest (e.g.,
  Species, Genus, Family, Order, Class, Phylum, Kingdom).

- max_level:

  Integer or `NULL`. Maximum number of taxonomic levels (above Species)
  to include in the product formula. When `NULL` (default), all
  available levels are used. When set to `"auto"`, the function
  automatically detects the highest informative level (where Deng
  entropy \> 0 at nk=0) and truncates the product there. A positive
  integer limits to that many levels (e.g., `max_level = 2` uses only
  Genus and Family).

## Value

A named list with components:

- uTO:

  Unweighted taxonomic diversity (all levels)

- TO:

  Weighted taxonomic diversity (all levels)

- uTO_plus:

  Unweighted taxonomic distance (all levels)

- TO_plus:

  Weighted taxonomic distance (all levels)

- uTO_max:

  Unweighted taxonomic diversity (max informative level)

- TO_max:

  Weighted taxonomic diversity (max informative level)

- uTO_plus_max:

  Unweighted taxonomic distance (max informative level)

- TO_plus_max:

  Weighted taxonomic distance (max informative level)

- Ed_levels:

  Deng entropy at each taxonomic level (nk=0 slice)

- max_informative_level:

  Integer: highest level with Ed \> 0

## Details

The method uses the slicing procedure from Ozkan (2018). At each slice
(nk = 0, 1, ..., n_s), species with abundance \<= nk are removed. The
surviving species receive EQUAL weight (1/count) — abundance information
enters indirectly through which species survive each slice.

Deng entropy at each taxonomic level is computed using these equal
proportions, where the mass function m(Fi) = count_in_group /
total_count and \|Fi\| = number of species in that taxonomic group.

The core product formula at each slice is:

\$\$\prod\_{i=1}^{L} \left( w_i \left( \frac{(e^{E_d^S})^2}
{e^{E_d^i}} + 1 \right) \right)\$\$

where \\E_d^S\\ is the Deng entropy at species level and \\E_d^i\\ is
the Deng entropy at level i, computed using presence/absence (equal
weight) proportions.

pTO+ (taxonomic distance) uses only the nk=0 slice: \$\$pT_O^+ = \ln
\prod\_{i=1}^{L} \left( w_i \left( \frac{(e^{E_d^S})^2}{e^{E_d^i}} + 1
\right) \right)\$\$

pTO (taxonomic diversity) aggregates across all slices: \$\$pT_O = \ln
\left( \frac{\sum\_{k=0}^{n_s} (n_s - n_k) \prod\_{i=1}^{L} \left( w_i
\left( \frac{(e^{E_d^S})^2} {e^{E_d^i}} + 1 \right) \right)}{n_s + \sum
n_k} \right)\$\$

## References

Ozkan, K. (2018). A new proposed measure for estimating taxonomic
diversity. Turkish Journal of Forestry, 19(4), 336-346. DOI:
10.18182/tjf.441061

Deng, Y. (2016). Deng entropy. Chaos, Solitons & Fractals, 91, 549-553.

## See also

[`deng_entropy_level()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/deng_entropy_level.md)
for the core Deng entropy calculation,
[`pto_components()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/pto_components.md)
for a convenience wrapper returning a named vector,
[`delta()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/delta.md)
and
[`avtd()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/avtd.md)
for Clarke & Warwick alternatives.

## Examples

``` r
# Simple example with 5 species
comm <- c(sp1 = 4, sp2 = 2, sp3 = 3, sp4 = 1, sp5 = 2)
tax <- data.frame(
  Species = paste0("sp", 1:5),
  Genus   = c("G1", "G1", "G1", "G2", "G2"),
  Family  = c("F1", "F1", "F1", "F1", "F1"),
  stringsAsFactors = FALSE
)
ozkan_pto(comm, tax)
#> taxdiv -- Ozkan pTO Result
#> 
#>   All levels:
#>     uTO  : 2.632978 
#>     TO   : 3.300651 
#>     uTO+ : 3.060705 
#>     TO+  : 3.753852 
#> 
#>   Max informative level (level 1):
#>     uTO_max  : 2.632978 
#>     TO_max   : 3.300651 
#>     uTO+_max : 3.060705 
#>     TO+_max  : 3.753852 
#> 
#>   Deng entropy by level:
#>     Species: 1.609438
#>     Genus: 2.280003
#>     Family: 0

# With auto max-level detection
ozkan_pto(comm, tax, max_level = "auto")
#> taxdiv -- Ozkan pTO Result
#> 
#>   All levels:
#>     uTO  : 2.632978 
#>     TO   : 3.300651 
#>     uTO+ : 3.060705 
#>     TO+  : 3.753852 
#> 
#>   Max informative level (level 1):
#>     uTO_max  : 2.632978 
#>     TO_max   : 3.300651 
#>     uTO+_max : 3.060705 
#>     TO+_max  : 3.753852 
#> 
#>   Deng entropy by level:
#>     Species: 1.609438
#>     Genus: 2.280003
#>     Family: 0
```
