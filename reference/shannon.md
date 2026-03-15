# Shannon Diversity Index

Calculates the Shannon-Wiener diversity index (H') for a community,
optionally applying a bias correction for small samples.

## Usage

``` r
shannon(
  community,
  base = exp(1),
  correction = c("none", "miller_madow", "grassberger", "chao_shen")
)
```

## Arguments

- community:

  A numeric vector of species abundances (counts).

- base:

  The logarithm base. Default is `exp(1)` (natural log). Use `2` for
  bits.

- correction:

  Bias correction method. One of `"none"` (default, naive MLE),
  `"miller_madow"`, `"grassberger"`, or `"chao_shen"`. See Details.

## Value

A numeric value representing the Shannon diversity index.

## Details

The naive (MLE) Shannon index is calculated as: \$\$H' =
-\sum\_{i=1}^{S} p_i \ln(p_i)\$\$ where \\p_i = n_i / N\\ is the
proportion of species \\i\\, \\N\\ is the total number of individuals,
and \\S\\ is the number of species observed.

The MLE estimator has a known negative bias that is significant for
small samples. Three bias-correction methods are available:

**Miller-Madow** (1955): Adds a first-order bias correction term:
\$\$H\_{MM} = H\_{MLE} + \frac{S\_{obs} - 1}{2N}\$\$

**Grassberger** (2003): Uses the digamma function instead of the
logarithm: \$\$H_G = \ln(N) - \frac{1}{N} \sum_i n_i \psi(n_i)\$\$ where
\\\psi\\ is the digamma function.

**Chao-Shen** (2003): Applies a Good-Turing coverage correction with
Horvitz-Thompson weighting: \$\$\hat{C} = 1 - f_1 / N\$\$ \$\$H\_{CS} =
-\sum_i \frac{\hat{p}\_i \ln \hat{p}\_i}{1 - (1 - \hat{p}\_i)^N}\$\$
where \\\hat{p}\_i = \hat{C} \cdot n_i / N\\ and \\f_1\\ is the number
of singletons.

Bias corrections require integer abundance counts. A warning is issued
if non-integer values are detected with `correction != "none"`.

## References

Miller, G.A. & Madow, W.G. (1954). On the maximum likelihood estimate of
the Shannon-Wiener index of diversity. AFCRC-TR-54-75.

Grassberger, P. (2003). Entropy estimates from insufficient samplings.
arXiv:physics/0307138.

Chao, A. & Shen, T.-J. (2003). Nonparametric estimation of Shannon's
index of diversity when there are unseen species in sample.
Environmental and Ecological Statistics, 10, 429-443.

## See also

[`simpson()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/simpson.md)
for Simpson diversity,
[`deng_entropy_level()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/deng_entropy_level.md)
for Deng entropy (a generalization of Shannon).

## Examples

``` r
comm <- c(10, 5, 8, 3, 12)
shannon(comm)
#> [1] 1.510657
shannon(comm, correction = "miller_madow")
#> [1] 1.563289
shannon(comm, correction = "grassberger")
#> [1] 1.578282
shannon(comm, correction = "chao_shen")
#> [1] 1.521172
```
