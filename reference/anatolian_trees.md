# Anatolian Forest Trees: Multi-Site Species Data

A data frame containing 20 tree species from Anatolian forests,
distributed across three sample plots with varying community
compositions. Species abundances follow the Westhoff & van der Maarel
(1973) scale (1–9). Taxonomic classification includes seven ranks from
species to kingdom.

## Usage

``` r
anatolian_trees
```

## Format

A data frame with 33 rows and 9 columns:

- Site:

  Sample plot name (character)

- Species:

  Binomial species name with underscore separator (character)

- Genus:

  Genus (character)

- Family:

  Family (character)

- Order:

  Order (character)

- Class:

  Class (character)

- Phylum:

  Phylum / Division (character)

- Kingdom:

  Kingdom (character)

- Abundance:

  Westhoff abundance value, integer 1–9 (numeric)

## Details

The three sites represent different forest types:

- Karisik_Orman:

  Mixed forest – both conifers and broadleaves (12 species)

- Yaprakli_Orman:

  Broadleaf-dominated forest (13 species)

- Konifer_Orman:

  Conifer-dominated forest (8 species)

This dataset can be used directly with
[`batch_analysis`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/batch_analysis.md)
for multi-site analysis:

    batch_analysis(anatolian_trees)

To extract a single community for use with
[`ozkan_pto`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/ozkan_pto.md)
or
[`compare_indices`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/compare_indices.md):

    site1 <- anatolian_trees[anatolian_trees$Site == "Karisik_Orman", ]
    community <- setNames(site1$Abundance, site1$Species)
    tax_tree  <- site1[, c("Species", "Genus", "Family", "Order",
                            "Class", "Phylum", "Kingdom")]
    ozkan_pto(community, tax_tree)

## References

Westhoff, V. & van der Maarel, E. (1973). The Braun-Blanquet approach.
In: R.H. Whittaker (ed.), Ordination and classification of communities.
Handbook of Vegetation Science 5, 617–726.

## See also

[`batch_analysis`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/batch_analysis.md)
for multi-site analysis,
[`gazi_comm`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/gazi_comm.md)
and
[`gazi_gytk`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/gazi_gytk.md)
for a single-community example.

## Examples

``` r
data(anatolian_trees)
head(anatolian_trees)
#>            Site            Species   Genus   Family   Order         Class
#> 1 Karisik_Orman        Pinus_nigra   Pinus Pinaceae Pinales     Pinopsida
#> 2 Karisik_Orman      Cedrus_libani  Cedrus Pinaceae Pinales     Pinopsida
#> 3 Karisik_Orman Abies_nordmanniana   Abies Pinaceae Pinales     Pinopsida
#> 4 Karisik_Orman     Quercus_cerris Quercus Fagaceae Fagales Magnoliopsida
#> 5 Karisik_Orman    Quercus_petraea Quercus Fagaceae Fagales Magnoliopsida
#> 6 Karisik_Orman   Fagus_orientalis   Fagus Fagaceae Fagales Magnoliopsida
#>          Phylum Kingdom Abundance
#> 1     Pinophyta Plantae         5
#> 2     Pinophyta Plantae         4
#> 3     Pinophyta Plantae         3
#> 4 Magnoliophyta Plantae         7
#> 5 Magnoliophyta Plantae         6
#> 6 Magnoliophyta Plantae         8

# Multi-site analysis
batch_analysis(anatolian_trees)
#> taxdiv -- Batch Analysis
#>   Sites: 3 
#>   Indices: 14 
#> 
#>            Site N_Species  Shannon  Simpson    Delta Delta_star     AvTD
#>   Karisik_Orman        12 2.317570 0.889414 3.826087   4.208289 4.575758
#>   Konifer_Orman         8 1.868458 0.824142 2.526882   2.967172 4.071429
#>  Yaprakli_Orman        13 2.399399 0.897929 2.822775   3.083196 3.435897
#>     VarTD      uTO       TO uTO_plus  TO_plus  uTO_max   TO_max uTO_plus_max
#>  2.335170 8.984158 15.56321 9.717193 16.29644 9.551268 16.13051    10.393644
#>  3.066327 7.236543 13.80194 8.416052 14.99530 8.400411 14.97703     9.280162
#>  0.733070 6.887741 10.05336 7.872282 11.05034 7.468801 10.64251     8.506628
#>  TO_plus_max
#>     16.97290
#>     15.85941
#>     11.68468

# Single site extraction
site1 <- anatolian_trees[anatolian_trees$Site == "Karisik_Orman", ]
comm  <- setNames(site1$Abundance, site1$Species)
tax   <- site1[, c("Species", "Genus", "Family", "Order",
                    "Class", "Phylum", "Kingdom")]
ozkan_pto(comm, tax)
#> taxdiv -- Ozkan pTO Result
#> 
#>   All levels:
#>     uTO  : 8.984158 
#>     TO   : 15.56321 
#>     uTO+ : 9.717193 
#>     TO+  : 16.29644 
#> 
#>   Max informative level (level 5):
#>     uTO_max  : 8.984158 
#>     TO_max   : 15.56321 
#>     uTO+_max : 9.717193 
#>     TO+_max  : 16.29644 
#> 
#>   Deng entropy by level:
#>     Species: 2.484907
#>     Genus: 2.552484
#>     Family: 3.137316
#>     Order: 3.685721
#>     Class: 5.233373
#>     Phylum: 5.233373
#>     Kingdom: 0
```
