# Compare All Diversity Indices Side by Side

Computes all available diversity indices for one or more communities and
returns them in a single data frame. Optionally produces a grouped bar
plot for visual comparison.

## Usage

``` r
compare_indices(
  communities,
  tax_tree,
  correction = c("none", "miller_madow", "grassberger", "chao_shen"),
  plot = FALSE
)
```

## Arguments

- communities:

  A named list of community vectors (named numeric), or a single named
  numeric vector. When a single vector is provided, it is wrapped in a
  list with name "Community".

- tax_tree:

  A data frame representing the taxonomic hierarchy, as produced by
  [`build_tax_tree`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/build_tax_tree.md).

- correction:

  Bias correction for the Shannon index. One of `"none"` (default),
  `"miller_madow"`, `"grassberger"`, or `"chao_shen"`. Passed to
  [`shannon()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/shannon.md).
  See
  [`shannon()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/shannon.md)
  for details.

- plot:

  Logical. If `TRUE` and ggplot2 is available, returns a list with both
  the data frame and a ggplot object. Default is `FALSE`.

## Value

If `plot = FALSE`, a data frame with communities as rows and indices as
columns. If `plot = TRUE`, a list with two elements:

- table:

  The data frame of index values.

- plot:

  A `ggplot` object showing a grouped bar chart.

## Details

The function calculates the following indices:

- **Shannon**: Shannon-Wiener entropy
  ([`shannon`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/shannon.md))

- **Simpson**: Gini-Simpson index
  ([`simpson`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/simpson.md))

- **Delta**: Clarke & Warwick taxonomic diversity
  ([`delta`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/delta.md))

- **Delta_star**: Clarke & Warwick taxonomic distinctness
  ([`delta_star`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/delta_star.md))

- **AvTD**: Average taxonomic distinctness
  ([`avtd`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/avtd.md))

- **VarTD**: Variation in taxonomic distinctness
  ([`vartd`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/vartd.md))

- **uTO**: Unweighted taxonomic diversity (Ozkan pTO, all levels)

- **TO**: Weighted taxonomic diversity (Ozkan pTO, all levels)

- **uTO_plus**: Unweighted taxonomic distance (Ozkan pTO, all levels)

- **TO_plus**: Weighted taxonomic distance (Ozkan pTO, all levels)

- **uTO_max**: Unweighted taxonomic diversity (informative levels)

- **TO_max**: Weighted taxonomic diversity (informative levels)

- **uTO_plus_max**: Unweighted taxonomic distance (informative levels)

- **TO_plus_max**: Weighted taxonomic distance (informative levels)

## See also

[`shannon`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/shannon.md),
[`simpson`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/simpson.md),
[`delta`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/delta.md),
[`delta_star`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/delta_star.md),
[`avtd`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/avtd.md),
[`vartd`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/vartd.md),
[`ozkan_pto`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/ozkan_pto.md),
[`pto_components`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/pto_components.md)

## Examples

``` r
tax <- build_tax_tree(
  species = c("sp1", "sp2", "sp3", "sp4"),
  Genus   = c("G1", "G1", "G2", "G2"),
  Family  = c("F1", "F1", "F1", "F2"),
  Order   = c("O1", "O1", "O1", "O1")
)

# Single community
comm <- c(sp1 = 10, sp2 = 20, sp3 = 15, sp4 = 5)
compare_indices(comm, tax)
#> taxdiv -- Index Comparison
#>   Communities: 1 
#>   Indices: Shannon, Simpson, Delta, Delta_star, AvTD, VarTD, uTO, TO, uTO_plus, TO_plus, uTO_max, TO_max, uTO_plus_max, TO_plus_max 
#> 
#>  Community N_Species  Shannon Simpson    Delta Delta_star AvTD    VarTD
#>  Community         4 1.279854     0.7 1.326531   1.857143    2 0.666667
#>       uTO       TO uTO_plus  TO_plus  uTO_max   TO_max uTO_plus_max TO_plus_max
#>  3.413215 5.066836  4.04615 5.837909 3.413215 5.066836      4.04615    5.837909

# Multiple communities
comm_list <- list(
  Site_A = c(sp1 = 10, sp2 = 20, sp3 = 15, sp4 = 5),
  Site_B = c(sp1 = 5, sp2 = 5, sp3 = 5, sp4 = 5)
)
compare_indices(comm_list, tax)
#> taxdiv -- Index Comparison
#>   Communities: 2 
#>   Indices: Shannon, Simpson, Delta, Delta_star, AvTD, VarTD, uTO, TO, uTO_plus, TO_plus, uTO_max, TO_max, uTO_plus_max, TO_plus_max 
#> 
#>  Community N_Species  Shannon Simpson    Delta Delta_star AvTD    VarTD
#>     Site_A         4 1.279854    0.70 1.326531   1.857143    2 0.666667
#>     Site_B         4 1.386294    0.75 1.578947   2.000000    2 0.666667
#>       uTO       TO uTO_plus  TO_plus  uTO_max   TO_max uTO_plus_max TO_plus_max
#>  3.413215 5.066836  4.04615 5.837909 3.413215 5.066836      4.04615    5.837909
#>  4.046150 5.837909  4.04615 5.837909 4.046150 5.837909      4.04615    5.837909
```
