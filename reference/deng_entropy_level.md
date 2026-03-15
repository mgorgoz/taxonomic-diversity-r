# Calculate Deng Entropy at a Single Taxonomic Level

Computes the Deng entropy (Ed) for a given set of group proportions at a
specific taxonomic level. This is the core entropy calculation from Deng
(2016), which generalizes Shannon entropy through the Dempster-Shafer
evidence theory framework.

## Usage

``` r
deng_entropy_level(
  abundances,
  group_sizes = NULL,
  correction = c("none", "miller_madow", "grassberger", "chao_shen")
)
```

## Arguments

- abundances:

  A numeric vector of abundances for each group (node) at the given
  taxonomic level.

- group_sizes:

  Optional integer vector of focal element sizes (`|Fi|`) for each
  group. At species level this is NULL (all sizes are 1, reducing to
  Shannon entropy). At higher taxonomic levels, each value represents
  the number of species within that group.

- correction:

  Bias correction method for Shannon entropy estimation. Only applied at
  species level (`group_sizes = NULL`). One of `"none"` (default),
  `"miller_madow"`, `"grassberger"`, or `"chao_shen"`. See
  [`shannon()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/shannon.md)
  for details. A warning is issued if correction is requested with
  non-NULL `group_sizes`.

## Value

A numeric value representing the Deng entropy at that level.

## Details

The Deng entropy is calculated as: \$\$E_d = -\sum\_{i} m(F_i) \ln
\frac{m(F_i)}{2^{\|F_i\|} - 1}\$\$

At species level, each focal element has cardinality 1, so Deng entropy
reduces to Shannon entropy: \$\$E_d^S = H = -\sum_i p_i \ln p_i\$\$

At higher levels (genus, family, etc.), \\\|F_i\|\\ equals the number of
species within each group, and the mass function is the normalized
proportion of total abundance in each group.

Bias correction is only meaningful at the species level where Deng
entropy equals Shannon entropy. At higher taxonomic levels the mass
function has a different structure and bias-correction formulas do not
apply.

## References

Deng, Y. (2016). Deng entropy. Chaos, Solitons & Fractals, 91, 549-553.

## See also

[`ozkan_pto()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/ozkan_pto.md)
which uses this function internally,
[`shannon()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/shannon.md)
for classical Shannon entropy and bias corrections.

## Examples

``` r
# Shannon entropy (species level, |Fi| = 1 for all)
deng_entropy_level(c(4, 2, 3, 1, 2, 3, 2, 2))
#> [1] 2.013806

# With bias correction at species level
deng_entropy_level(c(4, 2, 3, 1, 2), correction = "chao_shen")
#> [1] 1.704728

# Deng entropy at genus level with group sizes
deng_entropy_level(c(9, 3, 7), group_sizes = c(3, 2, 3))
#> [1] 2.825395
```
