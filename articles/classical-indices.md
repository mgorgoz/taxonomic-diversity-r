# Classical Diversity Indices: Shannon and Simpson

## What Are Classical Diversity Indices?

Shannon and Simpson indices are the two most widely used diversity
measures in ecology. They quantify how species abundances are
distributed within a community — essentially answering: **“How diverse
is this community based on how individuals are distributed among
species?”**

Both indices consider only the abundance distribution. They do not
account for taxonomic, phylogenetic, or functional relationships between
species. A community of 10 species from the same genus receives the same
score as a community of 10 species spanning 10 different orders.

This is why taxdiv pairs them with taxonomic measures — classical
indices capture the **abundance structure**, while taxonomic indices
capture the **hierarchical structure**.

``` r
library(taxdiv)

# Example community
community <- c(
  Quercus_coccifera    = 25,
  Quercus_infectoria   = 18,
  Pinus_brutia         = 30,
  Pinus_nigra          = 12,
  Juniperus_excelsa    = 8,
  Juniperus_oxycedrus  = 6,
  Arbutus_andrachne    = 15,
  Styrax_officinalis   = 4,
  Cercis_siliquastrum  = 3,
  Olea_europaea        = 10
)
```

## Shannon-Wiener Index (H’)

### The idea

Shannon entropy, borrowed from information theory (Shannon, 1948),
measures the **uncertainty** in predicting the species identity of a
randomly chosen individual. High uncertainty means high diversity — if
species are evenly distributed, it is hard to guess which species the
next individual belongs to.

### The formula

$$H\prime = - \sum\limits_{i = 1}^{S}p_{i}\ln\left( p_{i} \right)$$

where $p_{i}$ is the proportion of species $i$ and $S$ is the total
number of species.

### Key properties

- **Minimum**: $H\prime = 0$ when there is only one species (no
  uncertainty)
- **Maximum**: $H\prime = \ln(S)$ when all species have equal abundance
  (maximum uncertainty)
- **Units**: Measured in “nats” when using natural logarithm, “bits”
  when using log base 2
- **Sensitivity**: Moderately sensitive to both rare and abundant
  species

### Usage in taxdiv

``` r
# Default: natural logarithm
H <- shannon(community)
cat("Shannon H':", round(H, 4), "\n")
#> Shannon H': 2.0948
cat("Maximum possible H' for", length(community), "species:",
    round(log(length(community)), 4), "\n")
#> Maximum possible H' for 10 species: 2.3026
cat("Evenness (H'/H'max):", round(H / log(length(community)), 4), "\n")
#> Evenness (H'/H'max): 0.9098
```

### Bias Correction

When sample sizes are small, the observed Shannon index underestimates
the true value because rare species are likely missing from the sample.
taxdiv provides three correction methods:

``` r
cat("Uncorrected:  ", round(shannon(community), 4), "\n")
#> Uncorrected:   2.0948
cat("Miller-Madow: ", round(shannon(community, correction = "miller_madow"), 4), "\n")
#> Miller-Madow:  2.1292
cat("Grassberger:  ", round(shannon(community, correction = "grassberger"), 4), "\n")
#> Grassberger:   2.1338
cat("Chao-Shen:    ", round(shannon(community, correction = "chao_shen"), 4), "\n")
#> Chao-Shen:     2.1014
```

**Which correction to use?**

- **No correction**: When sample size is large relative to species
  richness (N \>\> S). This is the standard approach used in most
  published studies.
- **Miller-Madow**: Simple first-order correction. Adds $(S - 1)/2N$ to
  the estimate. Appropriate when you want a lightweight adjustment.
- **Grassberger**: Uses the digamma function for a more accurate
  correction. Performs well across a range of sample sizes.
- **Chao-Shen**: Uses Horvitz-Thompson estimation to account for unseen
  species. Best when you suspect many rare species are missing from the
  sample.

## Simpson Index

### The idea

Simpson’s index (Simpson, 1949) measures the **probability that two
randomly chosen individuals belong to the same species**. A community
dominated by one species has a high probability (low diversity); an even
community has a low probability (high diversity).

### Three variants

taxdiv provides all three common Simpson variants:

``` r
# Dominance (D): probability of same-species pair
D <- simpson(community, type = "dominance")
cat("Simpson dominance (D):    ", round(D, 4), "\n")
#> Simpson dominance (D):     0.1424

# Gini-Simpson (1-D): probability of different-species pair
GS <- simpson(community, type = "gini_simpson")
cat("Gini-Simpson (1-D):       ", round(GS, 4), "\n")
#> Gini-Simpson (1-D):        0.8576

# Inverse Simpson (1/D): effective number of species
inv <- simpson(community, type = "inverse")
cat("Inverse Simpson (1/D):    ", round(inv, 4), "\n")
#> Inverse Simpson (1/D):     7.0246
```

### Understanding the variants

| Variant                   | Formula              | Range  | Interpretation                                |
|---------------------------|----------------------|--------|-----------------------------------------------|
| **Dominance (D)**         | $\sum p_{i}^{2}$     | 0 to 1 | Higher = less diverse (one species dominates) |
| **Gini-Simpson (1-D)**    | $1 - \sum p_{i}^{2}$ | 0 to 1 | Higher = more diverse (common choice)         |
| **Inverse Simpson (1/D)** | $1/\sum p_{i}^{2}$   | 1 to S | Effective number of equally abundant species  |

The **inverse Simpson** is often the most intuitive: a value of 6.5
means the community is as diverse as one with 6.5 perfectly even
species.

## Shannon vs Simpson: When to Use Which?

``` r
# Even community
even <- c(sp1 = 20, sp2 = 20, sp3 = 20, sp4 = 20, sp5 = 20)

# Uneven community (same species, different abundances)
uneven <- c(sp1 = 90, sp2 = 4, sp3 = 3, sp4 = 2, sp5 = 1)

cat("=== Even community ===\n")
#> === Even community ===
cat("Shannon:", round(shannon(even), 4), "\n")
#> Shannon: 1.6094
cat("Simpson (1-D):", round(simpson(even, type = "gini_simpson"), 4), "\n\n")
#> Simpson (1-D): 0.8

cat("=== Uneven community ===\n")
#> === Uneven community ===
cat("Shannon:", round(shannon(uneven), 4), "\n")
#> Shannon: 0.4531
cat("Simpson (1-D):", round(simpson(uneven, type = "gini_simpson"), 4), "\n")
#> Simpson (1-D): 0.187
```

**Key difference**: Shannon is more sensitive to **rare species**
(because of the logarithm), while Simpson is more sensitive to
**dominant species** (because of the squaring). When a community has
many rare species, Shannon will detect them; Simpson may not.

| Scenario                                    | Better index               |
|---------------------------------------------|----------------------------|
| Comparing sites with different rare species | Shannon                    |
| Detecting dominance shifts                  | Simpson                    |
| Need sample-size independence               | Neither (use AvTD)         |
| Need taxonomic information                  | Neither (use pTO or Delta) |

## The Limitation: Why You Need Taxonomic Indices Too

Classical indices treat all species as interchangeable. Consider:

``` r
# Community A: 5 species from 5 different orders
comm_A <- c(sp1 = 20, sp2 = 20, sp3 = 20, sp4 = 20, sp5 = 20)

# Community B: 5 species from the same genus
comm_B <- c(sp6 = 20, sp7 = 20, sp8 = 20, sp9 = 20, sp10 = 20)

cat("Community A (5 orders)  - Shannon:", round(shannon(comm_A), 4), "\n")
#> Community A (5 orders)  - Shannon: 1.6094
cat("Community B (1 genus)   - Shannon:", round(shannon(comm_B), 4), "\n")
#> Community B (1 genus)   - Shannon: 1.6094
cat("Identical scores, yet A is far more taxonomically diverse.\n")
#> Identical scores, yet A is far more taxonomically diverse.
```

This is exactly why taxdiv includes Clarke & Warwick and Ozkan pTO
indices — they incorporate the taxonomic hierarchy to distinguish
between these communities. See the [Clarke &
Warwick](https://mgorgoz.github.io/taxonomic-diversity-r/articles/clarke-warwick.md)
and [Ozkan
pTO](https://mgorgoz.github.io/taxonomic-diversity-r/articles/ozkan-pto.md)
articles for details.

## References

- Shannon, C.E. (1948). A mathematical theory of communication. *Bell
  System Technical Journal*, 27(3), 379-423.
- Simpson, E.H. (1949). Measurement of diversity. *Nature*, 163, 688.
- Chao, A. & Shen, T.-J. (2003). Nonparametric estimation of Shannon’s
  index of diversity when there are unseen species in sample.
  *Environmental and Ecological Statistics*, 10, 429-443.
