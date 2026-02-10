# Taxonomic Diversity Formulas

> Mathematical foundations for the `taxdiv` R package.
> Extracted from the original publications.

---

## 1. Deng Entropy (Deng, 2016)

Deng Entropy is a generalized form of Shannon Entropy, developed within the
Dempster-Shafer Evidence Theory (DSET) framework. It measures uncertainty
of basic probability assignment (BPA).

$$
E_d = -\sum_{i} m(F_i) \ln \frac{m(F_i)}{2^{|F_i|} - 1}
$$

Where:
- $m(F_i)$ is the mass function (proportional value) assigned to focal element $F_i$
- $|F_i|$ is the cardinality (number of elements) of $F_i$
- $2^{|F_i|} - 1$ represents the potential number of states within $F_i$

**Key property:** When belief is assigned only to singletons (i.e., $|F_i| = 1$ for all $i$), Deng Entropy reduces to Shannon Entropy:

$$
E_d = -\sum_{i} m(\theta_i) \ln m(\theta_i) = H
$$

### Mass Function Calculation

For the taxonomic diversity index, the mass function is calculated through a slicing procedure:

$$
m(F_i) = m(F_i)^0 \big/ \sum m(F_i)^0
$$

Where $m(F_i)^0$ is the raw (unnormalized) mass value.

---

## 2. Özkan's Taxonomic Diversity Index — $pT_O$ (Özkan, 2018)

### 2.1 Main Formula (Equation 4 in Özkan, 2018)

The Deng Entropy-based taxonomic diversity measure ($pT_O$) is calculated as:

$$
pT_O = \ln \left( \frac{\sum_{k=0}^{n_s} (n_s - n_k) \prod_{i=1}^{7} \left( w_i \left( \frac{(e^{E_{d_S}})^2}{e^{E_{d_i}}} + 1 \right) \right)}{n_s + \sum n_k} \right), \quad n_{od_{i>1}} > 1
$$

Where:
- $e = 2.71828$ (Euler's number)
- $i$ = taxonomic level index (1=Species → 7=Kingdom)
- $E_{d_S}$ = Deng Entropy at species level (equals Shannon Entropy $H$)
- $E_{d_i}$ = Deng Entropy at taxonomic level $i$
- $n_s$ = total number of slicing steps
- $n_k$ = the $k$-th step value
- $n_{od_{i>1}}$ = number of nodes at genus level and above
- $w_i$ = weight for taxonomic level $i$

### 2.2 Weight System

Weights increase from species to kingdom level:

| Level | $i$ | $w_i$ (weighted) | $w_i$ (unweighted) |
|-------|-----|-------------------|---------------------|
| Species (Tür) | 1 | 1 | 1 |
| Genus (Cins) | 2 | 2 | 1 |
| Family (Familya) | 3 | 3 | 1 |
| Order (Takım) | 4 | 4 | 1 |
| Class (Sınıf) | 5 | 5 | 1 |
| Phylum (Şube) | 6 | 6 | 1 |
| Kingdom (Alem) | 7 | 7 | 1 |

When weights are used: $pT_O$ is denoted as $T_O$ (weighted taxonomic diversity)
When $w_i = 1$ for all levels: $pT_O$ is denoted as $uT_O$ (unweighted taxonomic diversity)

### 2.3 Slicing Procedure

The slicing procedure determines $n_s$ and the step values:

1. **Step 1** ($n_k = 0$): All species with abundance > 0 are marked as "1" (present)
2. **Step 2** ($n_k = 1$): Subtract 1 from each species' abundance. Species with remaining abundance > 0 → "1", otherwise → "0"
3. **Step 3** ($n_k = 2$): Subtract 2 from each species' abundance. Same rule applies.
4. Continue until no species has positive abundance remaining.
5. $n_s$ = total number of steps where at least one species is present.

**Example:** For abundances [4, 2, 3, 1, 2, 3, 2, 2]:
- Step 1: [1,1,1,1,1,1,1,1] → $(n_s - n_k) = 4$
- Step 2: [1,1,1,0,1,1,1,1] → $(n_s - n_k) = 3$
- Step 3: [1,0,1,0,0,1,0,0] → $(n_s - n_k) = 2$
- Step 4: [1,0,0,0,0,0,0,0] → $(n_s - n_k) = 1$
- Step 5: all zero → stop. $n_s = 4$

### 2.4 Deng Entropy at Each Taxonomic Level

$E_d$ is computed separately at each level of the Linnean hierarchy:
- $E_{d_S}$ = Deng Entropy at Species level (equals Shannon $H$, since $|F_i|=1$)
- $E_{d_G}$ = Deng Entropy at Genus level
- $E_{d_F}$ = Deng Entropy at Family level
- $E_{d_O}$ = Deng Entropy at Order level
- $E_{d_C}$ = Deng Entropy at Class level
- $E_{d_P}$ = Deng Entropy at Phylum level
- $E_{d_K}$ = Deng Entropy at Kingdom level (always 0 if only 1 kingdom)

**Termination rule:** If any level $i > 1$ has only one node ($n_{od_{i>1}} = 1$), then $E_{d_i} = 0$ and computation stops at that level.

### 2.5 Taxonomic Distance — $pT_O^+$ (Equation 5 in Özkan, 2018)

For presence/absence data (ignoring abundances), the taxonomic distance is:

$$
pT_O^+ = \prod_{i=1}^{7} \left( w_i \left( \frac{(e^{E_{d_S}})^2}{e^{E_{d_i}}} + 1 \right) \right), \quad n_{od_{i>1}} > 1
$$

This is essentially the core of the $pT_O$ formula when $n_s = 1$ and $n_k = 0$.

**Notation:**
- $T_O^+$ = weighted taxonomic distance
- $uT_O^+$ = unweighted taxonomic distance

### 2.6 Four Components Summary

| Component | Symbol | Description | Data Type |
|-----------|--------|-------------|-----------|
| Weighted taxonomic diversity | $T_O$ | Uses abundance + taxonomic weights | Abundance |
| Unweighted taxonomic diversity | $uT_O$ | Uses abundance, all $w_i = 1$ | Abundance |
| Weighted taxonomic distance | $T_O^+$ | Uses presence/absence + weights | Presence/Absence |
| Unweighted taxonomic distance | $uT_O^+$ | Uses presence/absence, all $w_i = 1$ | Presence/Absence |

---

## 3. Clarke & Warwick Indices (Warwick & Clarke, 1995)

### 3.1 Taxonomic Diversity — $\Delta$ (Equation 1)

$$
\Delta = \frac{\sum \sum_{i<j} w_{ij} x_i x_j + \sum_i 0 \cdot x_i(x_i - 1)/2}{\sum \sum_{i<j} x_i x_j + \sum_i x_i(x_i - 1)/2}
$$

Where:
- $x_i$ = abundance of species $i$
- $w_{ij}$ = distinctness weight (path length) between species $i$ and $j$
- The numerator zero term emphasizes: path length within same species = 0

### 3.2 Taxonomic Distinctness — $\Delta^*$ (Equation 2)

$$
\Delta^* = \frac{\sum \sum_{i<j} w_{ij} x_i x_j}{\sum \sum_{i<j} x_i x_j} = \frac{\sum w_k f_k}{\sum f_k}
$$

Where:
- $f_k$ = sum of cross-products of counts from all species pairs connected at level $k$
- $w_k$ = path weight for level $k$

**Weight scale (linear):**
- $w_1 = 1$ (species in same genus)
- $w_2 = 2$ (same order, different family)
- $w_3 = 3$ (same class, different order)
- ... up to $w_K$ (different phyla)

### 3.3 Average Taxonomic Distinctness — $\Delta^+$ (Clarke & Warwick, 1998)

For presence/absence data:

$$
\Delta^+ = \frac{\sum \sum_{i<j} w_{ij}}{S(S-1)/2}
$$

Where $S$ = number of species.

### 3.4 Variation in Taxonomic Distinctness — $\Lambda^+$ (Clarke & Warwick, 2001)

$$
\Lambda^+ = \frac{\sum \sum_{i<j} (w_{ij} - \Delta^+)^2}{S(S-1)/2}
$$

---

## 4. Relationship Between Indices

At species level, Deng Entropy equals Shannon Entropy:
$$E_{d_S} = H' = -\sum_i p_i \ln p_i$$

The $pT_O$ components show positive correlations with species richness ($S$) but are **not** indicators of species richness — they capture taxonomic hierarchy information that $S$ alone cannot represent.

From Özkan (2018), correlation values for 8 sample complexes:
- $r_{S-uT_O} = 0.965$
- $r_{S-T_O} = 0.827$
- $r_{S-uT_O^+} = 0.946$
- $r_{S-T_O^+} = 0.785$

---

## References

- Clarke, K.R. & Warwick, R.M. (1998). A taxonomic distinctness index and its statistical properties. *Journal of Applied Ecology*, 35, 523-531.
- Clarke, K.R. & Warwick, R.M. (2001). A further biodiversity index applicable to species lists: variation in taxonomic distinctness. *Marine Ecology Progress Series*, 216, 265-278.
- Deng, Y. (2016). Deng entropy. *Chaos, Solitons & Fractals*, 91, 549-553.
- Özkan, K. (2018). Taksonomik çeşitliliğin belirlenmesi için yeni önerilen bir eşitlik. *Turkish Journal of Forestry*, 19(4), 336-346. DOI: 10.18182/tjf.441061
- Warwick, R.M. & Clarke, K.R. (1995). New 'biodiversity' measures reveal a decrease in taxonomic distinctness with increasing stress. *Marine Ecology Progress Series*, 129, 301-305.
