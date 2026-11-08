# Mediterranean Forest: Multi-Site Vegetation Data

A data frame of 72 woody and herbaceous plant species recorded across 92
sample plots in a Mediterranean forest / maquis ecosystem. Species
abundances follow a cover-abundance scale (1–9). Taxonomic
classification includes seven ranks from species to kingdom.

## Usage

``` r
med_forest
```

## Format

A data frame with 1545 rows and 9 columns:

- Site:

  Sample plot identifier, `plot_01`–`plot_92` (character)

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

  Cover-abundance value, integer 1–9 (numeric)

## Details

The data originate from a plot-based vegetation survey. Each of the 92
plots (`plot_01`–`plot_92`) records the cover-abundance of the species
present, yielding 1545 species-by-plot observations. Per-plot richness
ranges from 7 to 30 species (mean approximately 17).

Taxonomic coverage: 92 plots, 72 species, 66 genera, 32 families, 20
orders, 4 classes (Magnoliopsida, Liliopsida, Pinopsida,
Polypodiopsida), all within Tracheophyta (Plantae).

This dataset can be used directly with
[`batch_analysis`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/batch_analysis.md)
for multi-site analysis:

    batch_analysis(med_forest)

To extract a single plot for use with
[`ozkan_pto`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/ozkan_pto.md)
or
[`compare_indices`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/compare_indices.md):


    plot1     <- med_forest[med_forest$Site == "plot_01", ]
    community <- setNames(plot1$Abundance, plot1$Species)
    tax_tree  <- plot1[, c("Species", "Genus", "Family", "Order",
                           "Class", "Phylum", "Kingdom")]
    ozkan_pto(community, tax_tree)

## See also

[`batch_analysis`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/batch_analysis.md)
for multi-site analysis,
[`anatolian_trees`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/anatolian_trees.md)
for a smaller multi-site example,
[`campus_comm`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/campus_comm.md)
and
[`campus_gytk`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/campus_gytk.md)
for a single-community example.

## Examples

``` r
data(med_forest)
head(med_forest)
#>      Site               Species     Genus        Family       Order
#> 1 plot_01     Arbutus_andrachne   Arbutus     Ericaceae    Ericales
#> 2 plot_01 Asparagus_acutifolius Asparagus  Asparagaceae Asparagales
#> 3 plot_01     Ceratonia_siliqua Ceratonia      Fabaceae     Fabales
#> 4 plot_01       Cistus_creticus    Cistus     Cistaceae    Malvales
#> 5 plot_01     Cotinus_coggygria   Cotinus Anacardiaceae  Sapindales
#> 6 plot_01    Crataegus_monogyna Crataegus      Rosaceae     Rosales
#>           Class       Phylum Kingdom Abundance
#> 1 Magnoliopsida Tracheophyta Plantae         4
#> 2    Liliopsida Tracheophyta Plantae         4
#> 3 Magnoliopsida Tracheophyta Plantae         5
#> 4 Magnoliopsida Tracheophyta Plantae         5
#> 5 Magnoliopsida Tracheophyta Plantae         4
#> 6 Magnoliopsida Tracheophyta Plantae         4
length(unique(med_forest$Site))    # 92 plots
#> [1] 92
length(unique(med_forest$Species)) # 72 species
#> [1] 72

# Single plot extraction
plot1 <- med_forest[med_forest$Site == "plot_01", ]
comm  <- setNames(plot1$Abundance, plot1$Species)
tax   <- plot1[, c("Species", "Genus", "Family", "Order",
                    "Class", "Phylum", "Kingdom")]
ozkan_pto(comm, tax)
#> taxdiv -- Ozkan pTO Result
#> 
#>   All levels:
#>     uTO  : 11.42649 
#>     TO   : 16.21376 
#>     uTO+ : 11.62039 
#>     TO+  : 16.40788 
#> 
#>   Max informative level (level 4):
#>     uTO_max  : 11.42649 
#>     TO_max   : 16.21376 
#>     uTO+_max : 11.62039 
#>     TO+_max  : 16.40788 
#> 
#>   Deng entropy by level:
#>     Species: 3.091042
#>     Genus: 3.091042
#>     Family: 3.354025
#>     Order: 3.810515
#>     Class: 11.95889
#>     Phylum: 0
#>     Kingdom: 0

# \donttest{
# Multi-site analysis (reduce n_iter for a quick run)
batch_analysis(med_forest, n_iter = 100)
#> taxdiv -- Batch Analysis
#>   Sites: 92 
#>   Indices: 14 
#> 
#>     Site N_Species  Shannon  Simpson    Delta Delta_star     AvTD    VarTD
#>  plot_01        22 3.048622 0.950311 4.005360   4.172214 4.147186 0.402579
#>  plot_02        22 2.969729 0.943577 4.032364   4.230320 4.147186 0.402579
#>  plot_03        30 3.248865 0.956706 4.111555   4.262678 4.252874 0.450997
#>  plot_04        16 2.649498 0.923704 4.037203   4.322106 4.208333 0.548264
#>  plot_05        20 2.803143 0.931855 3.960876   4.201103 4.047368 0.413546
#>  plot_06        19 2.784268 0.932742 4.023740   4.268472 4.175439 0.682672
#>  plot_07        20 2.847504 0.935402 4.107951   4.345416 4.257895 0.454543
#>  plot_08        26 3.088315 0.948624 3.994431   4.172827 4.036923 0.509406
#>  plot_09        25 3.020219 0.944515 4.078507   4.279541 4.060000 0.489733
#>  plot_10        22 2.924331 0.940329 4.086189   4.305252 4.277056 0.486010
#>       uTO       TO uTO_plus  TO_plus  uTO_max   TO_max uTO_plus_max TO_plus_max
#>  11.42649 16.21376 11.62039 16.40788 11.42649 16.21376     11.62039    16.40788
#>  11.04168 15.82918 11.62039 16.40788 11.04168 15.82918     11.62039    16.40788
#>  11.63946 16.42695 12.45736 17.24486 11.63946 16.42695     12.45736    17.24486
#>  10.27902 15.06652 10.63379 15.42128 10.27902 15.06652     10.63379    15.42128
#>  10.40382 15.19131 11.06113 15.84862 10.42505 15.21254     11.13597    15.92346
#>  10.30577 15.09326 10.50537 15.29286 10.50143 15.28892     10.86950    15.65699
#>  10.77870 15.56619 11.42499 16.21248 10.82023 15.60773     11.47943    16.26693
#>  10.50674 15.29423 11.21123 15.99872 10.66086 15.44836     11.38966    16.17715
#>  10.72225 15.50975 11.29581 16.08330 10.77619 15.56369     11.34324    16.13073
#>  10.79208 15.57957 11.53628 16.32377 10.79208 15.57957     11.53628    16.32377
#> ... ( 92  sites total)
# }
```
