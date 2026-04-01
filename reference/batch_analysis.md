# Batch Analysis from a Single Data Frame

Computes all diversity indices for one or more sample sites from a
single data frame (e.g., imported from Excel). The function
automatically detects the site column, taxonomic columns, and abundance
column, splits the data by site, and returns a summary data frame with
species count and 14 diversity indices per site.

## Usage

``` r
batch_analysis(
  data,
  site_column = NULL,
  tax_columns = NULL,
  abundance_column = "Abundance",
  correction = c("none", "miller_madow", "grassberger", "chao_shen"),
  full = TRUE,
  n_iter = 101L,
  seed = 42L,
  parallel = FALSE,
  n_cores = NULL
)
```

## Arguments

- data:

  A data frame containing species data. Must include at minimum a
  species column, at least one taxonomic rank column, and an abundance
  column. Optionally includes a site/plot column for multi-site
  analysis.

- site_column:

  Character string specifying the name of the site column. If `NULL`
  (default), the function searches for columns named `"Site"`, `"site"`,
  `"Alan"`, `"alan"`, `"Plot"`, or `"plot"`. If no such column is found,
  all data is treated as a single site.

- tax_columns:

  Character vector specifying the names of the taxonomic columns (from
  Species to highest rank). If `NULL` (default), the function
  auto-detects columns named `"Species"`, `"Genus"`, `"Family"`,
  `"Order"`, `"Class"`, `"Phylum"`, and `"Kingdom"` (case-insensitive).

- abundance_column:

  Character string specifying the name of the abundance column. Default
  is `"Abundance"` (case-insensitive match).

- correction:

  Bias correction for the Shannon index. One of `"none"` (default),
  `"miller_madow"`, `"grassberger"`, or `"chao_shen"`. Passed to
  [`shannon()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/shannon.md).
  See
  [`shannon()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/shannon.md)
  for details.

- full:

  Logical. If `TRUE` (default), run the full Ozkan pipeline (Run 1+2+3)
  using
  [`ozkan_pto_full()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/ozkan_pto_full.md)
  instead of deterministic-only
  [`pto_components()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/pto_components.md).
  This produces max values across all three runs, matching the Excel
  macro output. Set to `FALSE` for deterministic Run 1 only (faster but
  incomplete).

- n_iter:

  Number of stochastic iterations for Run 2 and Run 3 when
  `full = TRUE`. Default `101`. Ignored when `full = FALSE`.

- seed:

  Random seed for reproducibility when `full = TRUE`. Default `42`. Set
  to `NULL` for non-deterministic results. Ignored when `full = FALSE`.

- parallel:

  Logical. If `TRUE`, use parallel processing to compute indices for
  multiple sites concurrently. Default `FALSE`.

- n_cores:

  Number of CPU cores to use when `parallel = TRUE`. Default `NULL` uses
  up to 2 cores (CRAN policy limit).

## Value

A data frame with one row per site and columns: `Site`, `N_Species`,
`Shannon`, `Simpson`, `Delta`, `Delta_star`, `AvTD`, `VarTD`, `uTO`,
`TO`, `uTO_plus`, `TO_plus`, `uTO_max`, `TO_max`, `uTO_plus_max`,
`TO_plus_max`. When `full = TRUE`, the max columns reflect the maximum
across Run 1, 2, and 3.

## Details

When no site column is present (or all values are identical), the entire
data set is treated as a single community.

The function calculates the following indices per site:

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

- **uTO_max**: Unweighted taxonomic diversity (informative levels only)

- **TO_max**: Weighted taxonomic diversity (informative levels only)

- **uTO_plus_max**: Unweighted taxonomic distance (informative levels
  only)

- **TO_plus_max**: Weighted taxonomic distance (informative levels only)

## See also

[`compare_indices`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/compare_indices.md)
for analysis with pre-built community vectors,
[`build_tax_tree`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/build_tax_tree.md)
for building taxonomic trees manually,
[`ozkan_pto_full`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/ozkan_pto_full.md)
for the full 3-run pipeline on a single community.

## Examples

``` r
# Single-site data (no Site column) — full pipeline by default
# \donttest{
df <- data.frame(
  Species   = c("sp1", "sp2", "sp3", "sp4"),
  Genus     = c("G1", "G1", "G2", "G2"),
  Family    = c("F1", "F1", "F1", "F2"),
  Order     = c("O1", "O1", "O1", "O1"),
  Abundance = c(10, 20, 15, 5),
  stringsAsFactors = FALSE
)
batch_analysis(df)
#> taxdiv -- Batch Analysis
#>   Sites: 1 
#>   Indices: 14 
#> 
#>  Site N_Species  Shannon Simpson    Delta Delta_star AvTD    VarTD      uTO
#>   All         4 1.279854     0.7 1.326531   1.857143    2 0.666667 3.413215
#>        TO uTO_plus  TO_plus  uTO_max   TO_max uTO_plus_max TO_plus_max
#>  5.066836  4.04615 5.837909 3.413215 5.066836      4.04615    5.837909

# Multi-site data (with Site column)
df2 <- data.frame(
  Site      = c("A", "A", "A", "B", "B", "B"),
  Species   = c("sp1", "sp2", "sp3", "sp1", "sp3", "sp4"),
  Genus     = c("G1", "G1", "G2", "G1", "G2", "G2"),
  Family    = c("F1", "F1", "F1", "F1", "F1", "F2"),
  Order     = c("O1", "O1", "O1", "O1", "O1", "O1"),
  Abundance = c(10, 20, 15, 5, 25, 10),
  stringsAsFactors = FALSE
)
batch_analysis(df2)
#> taxdiv -- Batch Analysis
#>   Sites: 2 
#>   Indices: 14 
#> 
#>  Site N_Species  Shannon  Simpson    Delta Delta_star     AvTD    VarTD
#>     A         3 1.060857 0.641975 1.111111   1.692308 1.666667 0.222222
#>     B         3 0.900256 0.531250 0.833333   1.529412 2.000000 0.666667
#>       uTO       TO uTO_plus  TO_plus  uTO_max   TO_max uTO_plus_max TO_plus_max
#>  2.442117 3.132154 2.577008 3.270155 2.442117 3.132154     2.577008    3.270155
#>  2.804266 4.533563 3.767722 5.559482 2.991176 4.771455     3.767722    5.559482

# Fast mode: deterministic Run 1 only (no stochastic runs)
batch_analysis(df, full = FALSE)
#> taxdiv -- Batch Analysis
#>   Sites: 1 
#>   Indices: 14 
#> 
#>  Site N_Species  Shannon Simpson    Delta Delta_star AvTD    VarTD      uTO
#>   All         4 1.279854     0.7 1.326531   1.857143    2 0.666667 3.413215
#>        TO uTO_plus  TO_plus  uTO_max   TO_max uTO_plus_max TO_plus_max
#>  5.066836  4.04615 5.837909 3.413215 5.066836      4.04615    5.837909
# }
```
