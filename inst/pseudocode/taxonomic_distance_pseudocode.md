# Clarke & Warwick Taxonomic Distance/Distinctness — Pseudocode and I/O Specification

> Based on:
> - Warwick, R.M. & Clarke, K.R. (1995). New 'biodiversity' measures reveal
>   a decrease in taxonomic distinctness with increasing stress.
>   *Marine Ecology Progress Series*, 129, 301-305.
> - Clarke, K.R. & Warwick, R.M. (1998). A taxonomic distinctness index and
>   its statistical properties. *Journal of Applied Ecology*, 35, 523-531.
> - Clarke, K.R. & Warwick, R.M. (2001). A further biodiversity index
>   applicable to species lists: variation in taxonomic distinctness.
>   *Marine Ecology Progress Series*, 216, 265-278.
>
> Implementation: `taxdiv::delta()`, `taxdiv::delta_star()`, `taxdiv::avtd()`,
> `taxdiv::vartd()`, `taxdiv::tax_distance_matrix()`

---

## 1. Overview

Clarke & Warwick proposed a family of four indices that measure how
taxonomically spread apart the species in a community are. Unlike Shannon
or Simpson, these indices use the **path length** between species through
a Linnean classification tree.

```
Δ   (Delta)        Taxonomic diversity      — abundance-weighted, includes same-species pairs
Δ*  (Delta star)   Taxonomic distinctness    — abundance-weighted, excludes same-species pairs
Δ+  (Delta plus)   Average taxonomic dist.   — presence/absence only (AvTD)
Λ+  (Lambda plus)  Variation in tax. dist.   — variance of Δ+ (VarTD)
```

The first two (Δ, Δ*) require **abundance data**. The latter two (Δ+, Λ+)
require only **species lists** (presence/absence).

---

## 2. Path Weight (ω_ij) — The Foundation

All four indices depend on a pairwise distance ω_ij between species i and j.
This distance is defined by the taxonomic hierarchy:

```
ω_ij = weight of the FIRST taxonomic level at which species i and j
       share a common group (searching from lowest to highest rank)
```

### 2.1 Weight Scale

With L taxonomic levels above species:

| Level (k)  | Rank      | Default weight (ω_k) |
|------------|-----------|----------------------|
| 1          | Genus     | 1                    |
| 2          | Family    | 2                    |
| 3          | Order     | 3                    |
| ...        | ...       | ...                  |
| L          | Kingdom   | L                    |

If species share a Genus → ω_ij = 1 (close relatives)
If species share a Family but not Genus → ω_ij = 2
If species share nothing up to the highest level → ω_ij = L (maximum distance)

### 2.2 Pseudocode: Path Weight Computation

```
FUNCTION path_weight(species_i, species_j, tax_tree, weights):

    L = number of taxonomic levels (columns - 1)

    FOR k = 1 TO L:
        IF tax_tree[species_i, level_k] == tax_tree[species_j, level_k]:
            RETURN weights[k]

    # No common group found at any level
    RETURN weights[L]
```

### 2.3 Worked Example

```
Taxonomy:
  Species         Genus      Family      Order
  Quercus_robur   Quercus    Fagaceae    Fagales
  Fagus_orient.   Fagus      Fagaceae    Fagales
  Pinus_nigra     Pinus      Pinaceae    Pinales
  Abies_nordm.    Abies      Pinaceae    Pinales

Weights: [1, 2, 3]

Path weights:
  Q.robur - F.orient.  : Genus≠, Family=Fagaceae  → ω = 2
  Q.robur - P.nigra    : Genus≠, Family≠, Order≠  → ω = 3
  Q.robur - A.nordm.   : Genus≠, Family≠, Order≠  → ω = 3
  F.orient. - P.nigra  : Genus≠, Family≠, Order≠  → ω = 3
  F.orient. - A.nordm. : Genus≠, Family≠, Order≠  → ω = 3
  P.nigra - A.nordm.   : Genus≠, Family=Pinaceae  → ω = 2
```

---

## 3. Taxonomic Distance Matrix (`tax_distance_matrix`)

### 3.1 Input / Output

**Input:**

| Parameter   | Type             | Description                                    |
|-------------|------------------|------------------------------------------------|
| `tax_tree`  | data.frame       | Col 1 = Species, Col 2..n = ranks (low→high)  |
| `species`   | character or NULL| Subset of species to include (NULL = all)      |
| `weights`   | numeric or NULL  | Weight per level (NULL = 1, 2, 3, ...)         |

**Output:**

| Return      | Type                | Description                            |
|-------------|---------------------|----------------------------------------|
| `dist_mat`  | numeric matrix (S×S)| Symmetric distance matrix, diag = 0   |

### 3.2 Pseudocode

```
FUNCTION tax_distance_matrix(tax_tree, species = NULL, weights = NULL):

    # --- Step 1: Setup ---
    IF species is NULL:
        species = all species in tax_tree column 1
    ELSE:
        VALIDATE all species exist in tax_tree

    S = number of species
    L = ncol(tax_tree) - 1     # number of taxonomic levels

    IF weights is NULL:
        weights = [1, 2, 3, ..., L]
    VALIDATE length(weights) == L

    # --- Step 2: Initialize S × S matrix with zeros ---
    dist_mat = MATRIX(S, S, fill = 0)
    rownames = colnames = species

    # --- Step 3: Compute pairwise distances ---
    FOR i = 1 TO (S - 1):
        FOR j = (i + 1) TO S:
            d = weights[L]           # default: maximum distance
            FOR k = 1 TO L:
                IF tax_tree[species_i, level_k] == tax_tree[species_j, level_k]:
                    d = weights[k]
                    BREAK            # found first matching level
            dist_mat[i, j] = d
            dist_mat[j, i] = d       # symmetric

    RETURN dist_mat
```

### 3.3 Complexity

- Time: O(S² × L)
- Space: O(S²)

For typical ecological datasets (S < 500, L < 7), this is negligible.

---

## 4. Taxonomic Diversity — Δ (Delta) (`delta`)

### 4.1 Mathematical Formula

```
        Σ Σ_{i<j}  ω_ij · x_i · x_j
Δ = ─────────────────────────────────────────
    Σ Σ_{i<j} x_i · x_j  +  Σ_i x_i·(x_i-1)/2
```

- **Numerator:** Sum of ω_ij × x_i × x_j for all species pairs i < j
- **Denominator:** Total number of individual pairs (between-species + within-species)
- Same-species pairs contribute 0 to numerator (distance to self = 0)
- This is the average path length across ALL pairs of individuals

### 4.2 Input / Output

**Input:**

| Parameter    | Type           | Description                    |
|--------------|----------------|--------------------------------|
| `community`  | named numeric  | Species abundances (x_i)       |
| `tax_tree`   | data.frame     | Taxonomic hierarchy            |
| `weights`    | numeric or NULL| Level weights                  |

**Output:**

| Return  | Type    | Description              |
|---------|---------|--------------------------|
| `Δ`     | numeric | Taxonomic diversity      |

### 4.3 Pseudocode

```
FUNCTION delta(community, tax_tree, weights = NULL):

    # --- Step 1: Validation and setup ---
    REMOVE species with abundance = 0
    IF fewer than 2 species: RETURN 0

    S = number of species
    x = abundance vector
    L = ncol(tax_tree) - 1
    IF weights is NULL: weights = [1, 2, 3, ..., L]

    # --- Step 2: Compute numerator and denominator ---
    numerator = 0
    denom_cross = 0

    FOR i = 1 TO (S - 1):
        FOR j = (i + 1) TO S:
            ω_ij = path_weight(species_i, species_j, tax_tree, weights)
            numerator   += ω_ij × x[i] × x[j]
            denom_cross += x[i] × x[j]

    # Within-species pairs: contribute 0 to numerator, but count in denominator
    denom_same = SUM( x[i] × (x[i] - 1) / 2 ) for all i

    # --- Step 3: Result ---
    Δ = numerator / (denom_cross + denom_same)
    RETURN Δ
```

### 4.4 Worked Example

```
Community: sp1=5, sp2=3, sp3=2
Taxonomy:
  sp1  G1  F1
  sp2  G1  F1
  sp3  G2  F2
Weights: [1, 2]

Path weights:
  ω(sp1,sp2) = 1  (same Genus G1)
  ω(sp1,sp3) = 2  (different at all levels)
  ω(sp2,sp3) = 2  (different at all levels)

Numerator = 1×5×3 + 2×5×2 + 2×3×2
          = 15 + 20 + 12
          = 47

Denominator:
  cross = 5×3 + 5×2 + 3×2 = 15 + 10 + 6 = 31
  same  = 5×4/2 + 3×2/2 + 2×1/2 = 10 + 3 + 1 = 14
  total = 31 + 14 = 45

Δ = 47 / 45 = 1.044
```

---

## 5. Taxonomic Distinctness — Δ* (Delta star) (`delta_star`)

### 5.1 Mathematical Formula

```
         Σ Σ_{i<j}  ω_ij · x_i · x_j
Δ* = ─────────────────────────────────
         Σ Σ_{i<j}  x_i · x_j
```

The only difference from Δ: the denominator **excludes** within-species
pairs. This removes the dependence on evenness of abundances and focuses
purely on the taxonomic spread between different species.

### 5.2 Pseudocode

```
FUNCTION delta_star(community, tax_tree, weights = NULL):

    # --- Step 1: Validation and setup ---
    REMOVE species with abundance = 0
    IF fewer than 2 species: RETURN 0

    S = number of species
    x = abundance vector
    L = ncol(tax_tree) - 1
    IF weights is NULL: weights = [1, 2, 3, ..., L]

    # --- Step 2: Compute numerator and denominator ---
    numerator = 0
    denominator = 0

    FOR i = 1 TO (S - 1):
        FOR j = (i + 1) TO S:
            ω_ij = path_weight(species_i, species_j, tax_tree, weights)
            numerator   += ω_ij × x[i] × x[j]
            denominator += x[i] × x[j]

    IF denominator == 0: RETURN 0

    # --- Step 3: Result ---
    Δ* = numerator / denominator
    RETURN Δ*
```

### 5.3 Worked Example (same data as Section 4.4)

```
Numerator = 47  (same as Δ)
Denominator = 31  (cross-species pairs only, no same-species)

Δ* = 47 / 31 = 1.516
```

### 5.4 Difference Between Δ and Δ*

```
Δ  = average path length across ALL individual pairs (including same-species)
Δ* = average path length across DIFFERENT-species individual pairs only

Δ ≤ Δ* always holds (since same-species pairs contribute 0 distance to Δ)

When all species have abundance = 1: Δ = Δ* (no within-species pairs)
```

---

## 6. Average Taxonomic Distinctness — Δ+ (AvTD) (`avtd`)

### 6.1 Mathematical Formula

```
         Σ Σ_{i<j}  ω_ij
Δ+ = ─────────────────────
           S(S - 1) / 2
```

This is the **presence/absence** version of Δ*. Each species has weight 1
(no abundance information), so x_i × x_j = 1 for all pairs.

### 6.2 Input / Output

**Input:**

| Parameter  | Type           | Description                          |
|------------|----------------|--------------------------------------|
| `species`  | character      | Species names present in community   |
| `tax_tree` | data.frame     | Taxonomic hierarchy                  |
| `weights`  | numeric or NULL| Level weights                        |

**Output:**

| Return | Type    | Description                             |
|--------|---------|-----------------------------------------|
| `Δ+`  | numeric | Average taxonomic distinctness          |

### 6.3 Pseudocode

```
FUNCTION avtd(species, tax_tree, weights = NULL):

    # --- Step 1: Validation ---
    IF count(species) < 2: RAISE ERROR "need at least 2 species"

    # --- Step 2: Compute distance matrix ---
    dist_mat = tax_distance_matrix(tax_tree, species, weights)

    # --- Step 3: Average of upper triangle ---
    S = count(species)
    n_pairs = S × (S - 1) / 2
    sum_dist = SUM( dist_mat[i, j] for all i < j )

    Δ+ = sum_dist / n_pairs
    RETURN Δ+
```

### 6.4 Worked Example

```
Species present: sp1, sp2, sp3 (same taxonomy as Section 4.4)

Distance matrix (upper triangle):
  ω(sp1,sp2) = 1
  ω(sp1,sp3) = 2
  ω(sp2,sp3) = 2

S = 3
n_pairs = 3 × 2 / 2 = 3
sum_dist = 1 + 2 + 2 = 5

Δ+ = 5 / 3 = 1.667
```

### 6.5 Key Property

Δ+ is **independent of sampling effort** — it depends only on which species
are present, not how many individuals were counted. This makes it suitable
for comparing communities with unequal survey effort.

---

## 7. Variation in Taxonomic Distinctness — Λ+ (VarTD) (`vartd`)

### 7.1 Mathematical Formula

```
         Σ Σ_{i<j}  (ω_ij - Δ+)²
Λ+ = ───────────────────────────────
              S(S - 1) / 2
```

This is the **variance** of pairwise distances around Δ+. High Λ+ means
some species pairs are very close and others very distant (uneven
taxonomic spread). Low Λ+ means species are evenly spread across the
taxonomy.

### 7.2 Pseudocode

```
FUNCTION vartd(species, tax_tree, weights = NULL):

    # --- Step 1: Validation ---
    IF count(species) < 2: RAISE ERROR "need at least 2 species"

    # --- Step 2: Compute distance matrix ---
    dist_mat = tax_distance_matrix(tax_tree, species, weights)

    # --- Step 3: Extract upper triangle values ---
    S = count(species)
    n_pairs = S × (S - 1) / 2
    upper_vals = all dist_mat[i, j] where i < j

    # --- Step 4: Compute Δ+ (mean) ---
    Δ+ = SUM(upper_vals) / n_pairs

    # --- Step 5: Compute Λ+ (variance around Δ+) ---
    Λ+ = SUM( (upper_vals - Δ+)² ) / n_pairs
    RETURN Λ+
```

### 7.3 Worked Example

```
upper_vals = [1, 2, 2]
Δ+ = 5/3 = 1.667

Λ+ = [(1 - 1.667)² + (2 - 1.667)² + (2 - 1.667)²] / 3
   = [0.444 + 0.111 + 0.111] / 3
   = 0.667 / 3
   = 0.222
```

---

## 8. Relationship Summary

```
                          Abundance data required?
                          ─────────────────────────
                          YES              NO
                     ┌──────────────┬──────────────┐
  Includes same-     │              │              │
  species pairs?     │   Δ (Delta)  │     N/A      │
  YES                │              │              │
                     ├──────────────┼──────────────┤
  Excludes same-     │              │              │
  species pairs?     │  Δ* (Delta*) │  Δ+ (AvTD)   │
  NO                 │              │  Λ+ (VarTD)  │
                     └──────────────┴──────────────┘
```

### Ordering Properties

For a given community and taxonomy:
```
Δ ≤ Δ*         (always, since same-species pairs add 0 distance)
Δ+ = Δ* when all abundances = 1
```

### Function Dependency Chain

```
tax_distance_matrix(tax_tree, species, weights)
    ^                     ^                ^
    |                     |                |
    |    ┌────────────────┘                |
    |    |                                 |
delta(community, tax_tree, weights)        |
delta_star(community, tax_tree, weights)   |
    |                                      |
    |    ┌─────────────────────────────────┘
    |    |
avtd(species, tax_tree, weights)
vartd(species, tax_tree, weights)
```

Note: `delta()` and `delta_star()` compute path weights inline (without
calling `tax_distance_matrix()`), while `avtd()` and `vartd()` call
`tax_distance_matrix()` internally.

---

## 9. Edge Cases

| Condition                          | Δ    | Δ*   | Δ+   | Λ+   |
|------------------------------------|------|------|------|------|
| 0-1 species                       | 0    | 0    | ERROR| ERROR|
| All species same genus             | ω=1  | ω=1  | ω=1  | 0    |
| All species different at all levels| ω=L  | ω=L  | ω=L  | 0    |
| Equal abundances (all x_i = c)     | = Δ* | same | N/A  | N/A  |
| Single dominant species            | → 0  | ≠ 0  | N/A  | N/A  |

---

## 10. Comparison with Ozkan pTO

| Feature            | Clarke & Warwick       | Ozkan pTO                 |
|--------------------|------------------------|---------------------------|
| Entropy basis      | None (path length avg) | Deng entropy              |
| Hierarchy use      | Pairwise distance      | Level-wise entropy + product|
| Abundance handling | Direct weighting       | Slicing procedure         |
| Presence/absence   | Δ+ and Λ+             | uTO+ and TO+              |
| Components         | 4 (Δ, Δ*, Δ+, Λ+)     | 4 (uTO, TO, uTO+, TO+)   |
| Stochastic runs    | No                     | Yes (Run 2, Run 3)        |
| R functions        | `delta()`, `delta_star()`, `avtd()`, `vartd()` | `ozkan_pto()`, `ozkan_pto_resample()`, `ozkan_pto_sensitivity()` |

---

## References

- Warwick, R.M. & Clarke, K.R. (1995). New 'biodiversity' measures reveal
  a decrease in taxonomic distinctness with increasing stress. *Marine
  Ecology Progress Series*, 129, 301-305.
- Clarke, K.R. & Warwick, R.M. (1998). A taxonomic distinctness index and
  its statistical properties. *Journal of Applied Ecology*, 35, 523-531.
- Clarke, K.R. & Warwick, R.M. (2001). A further biodiversity index
  applicable to species lists: variation in taxonomic distinctness. *Marine
  Ecology Progress Series*, 216, 265-278.
