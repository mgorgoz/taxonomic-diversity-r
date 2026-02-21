# Deng Entropy — Pseudocode and I/O Specification

> Based on: Deng, Y. (2016). Deng entropy. *Chaos, Solitons & Fractals*, 91, 549-553.
> Implementation: `taxdiv::deng_entropy_level()`

---

## 1. Mathematical Definition

Deng Entropy is a generalization of Shannon Entropy within the Dempster-Shafer
Evidence Theory (DSET) framework:

```
Ed = -SUM_i [ m(Fi) * ln( m(Fi) / (2^|Fi| - 1) ) ]
```

Where:
- `m(Fi)` = mass function (normalized proportion) for focal element `Fi`
- `|Fi|`  = cardinality of focal element `Fi` (number of elements it contains)
- `ln`    = natural logarithm

**Special case (Species level):** When `|Fi| = 1` for all elements:
```
Ed = -SUM_i [ m(Fi) * ln( m(Fi) / (2^1 - 1) ) ]
   = -SUM_i [ m(Fi) * ln( m(Fi) / 1 ) ]
   = -SUM_i [ m(Fi) * ln( m(Fi) ) ]
   = H   (Shannon Entropy)
```

---

## 2. Input / Output Specification

### Input

| Parameter      | Type              | Description                                      | Required |
|----------------|-------------------|--------------------------------------------------|----------|
| `abundances`   | numeric vector    | Abundance (or count) of each group at the level  | Yes      |
| `group_sizes`  | integer vector    | Number of species in each group (`|Fi|`)         | No*      |

*If `group_sizes` is NULL, species level is assumed (`|Fi| = 1` for all).

### Output

| Return         | Type    | Description                              |
|----------------|---------|------------------------------------------|
| `Ed`           | numeric | Deng entropy value at the given level    |

### Edge Cases

| Condition                    | Return Value | Reason                              |
|------------------------------|--------------|-------------------------------------|
| No positive abundances       | 0            | No information                      |
| Single group (length = 1)    | 0            | No uncertainty with one element     |
| All `|Fi| = 1`               | Shannon H    | Deng reduces to Shannon             |

---

## 3. Pseudocode

```
FUNCTION deng_entropy_level(abundances, group_sizes = NULL):

    # --- Step 1: Input validation ---
    IF abundances is not numeric OR any value < 0:
        RAISE ERROR "abundances must be numeric and non-negative"

    # --- Step 2: Remove zero-abundance groups ---
    IF group_sizes is provided:
        keep = indices where abundances > 0
        abundances = abundances[keep]
        group_sizes = group_sizes[keep]
    ELSE:
        abundances = abundances[abundances > 0]

    # --- Step 3: Handle trivial cases ---
    IF length(abundances) == 0:
        RETURN 0
    IF length(abundances) == 1:
        RETURN 0

    # --- Step 4: Compute mass function (proportions) ---
    total = SUM(abundances)
    m = abundances / total
    # Now: m[i] = proportion of total abundance in group i

    # --- Step 5: Compute entropy ---
    IF group_sizes is NULL:
        # Species level: |Fi| = 1 for all -> Shannon Entropy
        # Ed = -SUM( m[i] * ln(m[i]) )
        Ed = -SUM( m * ln(m) )
    ELSE:
        # Higher taxonomic level: use Deng formula
        # Ed = -SUM( m[i] * ln( m[i] / (2^|Fi[i]| - 1) ) )
        FOR each group i:
            denominator[i] = 2^group_sizes[i] - 1
        Ed = -SUM( m * ln( m / denominator ) )

    RETURN Ed
```

---

## 4. Worked Example

### Example 1: Species Level (reduces to Shannon)

8 species with abundances: `[4, 2, 3, 1, 2, 3, 2, 2]`

```
total = 4+2+3+1+2+3+2+2 = 19
m = [4/19, 2/19, 3/19, 1/19, 2/19, 3/19, 2/19, 2/19]

group_sizes = NULL  ->  all |Fi| = 1

Ed = -SUM( m[i] * ln(m[i]) )
   = -(4/19)*ln(4/19) - (2/19)*ln(2/19) - (3/19)*ln(3/19)
     - (1/19)*ln(1/19) - (2/19)*ln(2/19) - (3/19)*ln(3/19)
     - (2/19)*ln(2/19) - (2/19)*ln(2/19)
   = 2.0025 (approximately)
```

This equals Shannon H for the same distribution.

### Example 2: Genus Level with Group Sizes

3 genera containing [9, 3, 7] species respectively (total = 19 species):

```
total = 9+3+7 = 19
m = [9/19, 3/19, 7/19]
group_sizes = [3, 2, 3]

denominators:
  2^3 - 1 = 7
  2^2 - 1 = 3
  2^3 - 1 = 7

Ed = -[ (9/19)*ln((9/19)/7) + (3/19)*ln((3/19)/3) + (7/19)*ln((7/19)/7) ]
   = -[ (9/19)*ln(9/133) + (3/19)*ln(3/57) + (7/19)*ln(7/133) ]
   = -[ (9/19)*(-2.693) + (3/19)*(-2.944) + (7/19)*(-2.944) ]
   = -[ -1.275 - 0.465 - 1.084 ]
   = 2.824 (approximately)
```

Note: Ed_genus > Ed_species is possible because the `2^|Fi| - 1` denominator
inflates entropy for groups containing multiple species. This is by design —
it captures the "potential diversity" within each focal element.

---

## 5. Relationship to Other Functions

```
deng_entropy_level()
    ^
    |   used internally by
    |
ozkan_pto()  -- computes Ed at each taxonomic level (Species, Genus, ..., Kingdom)
    |
    v
pto_components()  -- convenience wrapper returning named numeric vector
```

At each slice `nk` of the pTO formula, `deng_entropy_level()` is called for
every active taxonomic level using the equal-weight (presence-based) proportions
of the surviving species.

---

## References

- Deng, Y. (2016). Deng entropy. *Chaos, Solitons & Fractals*, 91, 549-553.
- Ozkan, K. (2018). A new proposed measure for estimating taxonomic diversity.
  *Turkish Journal of Forestry*, 19(4), 336-346. DOI: 10.18182/tjf.441061
