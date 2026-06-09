# Batch Analysis from a Single Data Frame

Computes selected diversity indices for one or more sample sites from a
single data frame (e.g., imported from Excel). The function
automatically detects the site column, taxonomic columns, and abundance
column, splits the data by site, and returns a summary data frame with
species count and diversity indices per site.

## Usage

``` r
batch_analysis(
  data,
  indices = c("classical", "clarke_warwick", "ozkan_pto"),
  site_column = NULL,
  tax_columns = NULL,
  abundance_column = "Abundance",
  correction = c("none", "miller_madow", "grassberger", "chao_shen"),
  full = TRUE,
  n_iter = 101L,
  seed = 42L,
  parallel = FALSE,
  n_cores = NULL,
  progress = TRUE,
  progress_fn = NULL
)
```

## Arguments

- data:

  A data frame containing species data. Must include at minimum a
  species column, at least one taxonomic rank column, and an abundance
  column. Optionally includes a site/plot column for multi-site
  analysis.

- indices:

  Character vector specifying which index groups to compute. One or more
  of `"classical"` (Shannon, Simpson), `"clarke_warwick"` (Delta,
  Delta\*, AvTD, VarTD), and `"ozkan_pto"` (uTO, TO, uTO+, TO+ with
  optional Run 2+3). Default is all three groups. Unambiguous
  abbreviations are allowed (e.g., `"clas"` for classical, `"clark"` for
  clarke_warwick, `"oz"` for ozkan_pto). Note that `"cl"` and `"cla"`
  are ambiguous and will produce an error.

- site_column:

  Character string specifying the name of the site column. If `NULL`
  (default), the function searches for columns named `"Site"`, `"site"`,
  `"Plot"`, or `"plot"`. If no such column is found, all data is treated
  as a single site.

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
  for details. Ignored when `"classical"` is not in `indices`.

- full:

  Logical. If `TRUE` (default), run the full Ozkan pipeline (Run 1+2+3)
  using
  [`ozkan_pto_full()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/ozkan_pto_full.md)
  instead of deterministic-only
  [`pto_components()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/pto_components.md).
  This produces max values across all three runs, matching the Excel
  macro output. Set to `FALSE` for deterministic Run 1 only (faster but
  incomplete). Ignored when `"ozkan_pto"` is not in `indices`.

- n_iter:

  Number of stochastic iterations for Run 2 and Run 3 when
  `full = TRUE`. Default `101`. Ignored when `full = FALSE` or when
  `"ozkan_pto"` is not in `indices`.

- seed:

  Random seed for reproducibility when `full = TRUE`. Default `42`. Set
  to `NULL` for non-deterministic results. Ignored when `full = FALSE`
  or when `"ozkan_pto"` is not in `indices`.

- parallel:

  Logical. If `TRUE`, use parallel processing to compute indices for
  multiple sites concurrently. Default `FALSE`.

- n_cores:

  Number of CPU cores to use when `parallel = TRUE`. Default `NULL` uses
  up to 2 cores (CRAN policy limit).

- progress:

  Logical. If `TRUE` (default), display a progress bar during sequential
  computation. Ignored when `parallel = TRUE`. Set to `FALSE` to
  suppress progress output.

- progress_fn:

  Optional callback function for custom progress reporting (e.g. Shiny).
  When provided, it is called after each site completes with named
  arguments `i` (current site index), `n` (total sites), and `site`
  (site name). Useful for integrating with
  [`shiny::withProgress()`](https://rdrr.io/pkg/shiny/man/withProgress.html).

## Value

A data frame with one row per site. Columns always include `Site` and
`N_Species`. Additional columns depend on the `indices` parameter:

- `"classical"`: `Shannon`, `Simpson`

- `"clarke_warwick"`: `Delta`, `Delta_star`, `AvTD`, `VarTD`

- `"ozkan_pto"`: `uTO`, `TO`, `uTO_plus`, `TO_plus`, `uTO_max`,
  `TO_max`, `uTO_plus_max`, `TO_plus_max`

## Details

When no site column is present (or all values are identical), the entire
data set is treated as a single community.

Three groups of indices are available via the `indices` parameter:

- `"classical"`:

  Shannon-Wiener entropy and Gini-Simpson index. These are species-level
  diversity measures that do not use taxonomic hierarchy.

- `"clarke_warwick"`:

  Delta, Delta\*, AvTD, and VarTD. Taxonomy-aware distinctness measures
  from Clarke & Warwick (1998).

- `"ozkan_pto"`:

  Deng entropy-based taxonomic diversity and distance (uTO, TO, uTO+,
  TO+) from Ozkan (2018). When `full = TRUE`, also computes max values
  across Run 1+2+3. The Westhoff-Maarel cover-abundance scale (integer
  values 1-9) is recommended for compatibility with the original paper,
  but any positive numeric abundance values are accepted.

## See also

[`compare_indices`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/compare_indices.md)
for analysis with pre-built community vectors,
[`build_tax_tree`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/build_tax_tree.md)
for building taxonomic trees manually,
[`ozkan_pto_full`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/ozkan_pto_full.md)
for the full 3-run pipeline on a single community.

## Examples

``` r
# All indices (default)
# \donttest{
df <- data.frame(
  Species   = c("sp1", "sp2", "sp3", "sp4"),
  Genus     = c("G1", "G1", "G2", "G2"),
  Family    = c("F1", "F1", "F1", "F2"),
  Order     = c("O1", "O1", "O1", "O1"),
  Abundance = c(4, 2, 3, 1),
  stringsAsFactors = FALSE
)
batch_analysis(df)
#> taxdiv -- Batch Analysis
#>   Sites: 1 
#>   Indices: 14 
#> 
#>  Site N_Species  Shannon Simpson    Delta Delta_star AvTD    VarTD      uTO
#>   All         4 1.279854     0.7 1.444444   1.857143    2 0.666667 3.357526
#>        TO uTO_plus  TO_plus  uTO_max   TO_max uTO_plus_max TO_plus_max
#>  5.002732  4.04615 5.837909 3.357526 5.002732      4.04615    5.837909

# Only classical indices (fast)
batch_analysis(df, indices = "classical")
#> taxdiv -- Batch Analysis
#>   Sites: 1 
#>   Indices: 2 
#> 
#>  Site N_Species  Shannon Simpson
#>   All         4 1.279854     0.7

# Classical + Clarke & Warwick (no pTO)
batch_analysis(df, indices = c("classical", "clarke_warwick"))
#> taxdiv -- Batch Analysis
#>   Sites: 1 
#>   Indices: 6 
#> 
#>  Site N_Species  Shannon Simpson    Delta Delta_star AvTD    VarTD
#>   All         4 1.279854     0.7 1.444444   1.857143    2 0.666667

# Only Ozkan pTO, deterministic Run 1
batch_analysis(df, indices = "ozkan_pto", full = FALSE)
#> taxdiv -- Batch Analysis
#>   Sites: 1 
#>   Indices: 8 
#> 
#>  Site N_Species      uTO       TO uTO_plus  TO_plus  uTO_max   TO_max
#>   All         4 3.357526 5.002732  4.04615 5.837909 3.357526 5.002732
#>  uTO_plus_max TO_plus_max
#>       4.04615    5.837909
# }
```
