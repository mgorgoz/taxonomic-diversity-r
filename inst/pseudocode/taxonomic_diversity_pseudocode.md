# Özkan pTO Taxonomic Diversity — Pseudocode and I/O Specification

> Based on: Özkan, K. (2018). A new proposed measure for estimating taxonomic
> diversity. *Turkish Journal of Forestry*, 19(4), 336-346.
> DOI: 10.18182/tjf.441061
>
> Implementation: `taxdiv::ozkan_pto()`, `taxdiv::ozkan_pto_resample()`,
> `taxdiv::ozkan_pto_sensitivity()`

---

## 1. Overview

The Özkan (2018) method produces four taxonomic diversity/distance components
through a three-stage pipeline:

```
Run 1 (İşlem 1):  Deterministic calculation     -> ozkan_pto()
Run 2 (İşlem 2):  Stochastic resampling (50%)   -> ozkan_pto_resample()
Run 3 (İşlem 3):  Sensitivity analysis           -> ozkan_pto_sensitivity()
```

Each run computes four components:

| Symbol  | Name                            | Data Type          | Weight   |
|---------|---------------------------------|--------------------|----------|
| uTO+    | Unweighted taxonomic distance   | Presence/Absence   | wi = 1   |
| TO+     | Weighted taxonomic distance     | Presence/Absence   | wi = i   |
| uTO     | Unweighted taxonomic diversity  | Abundance          | wi = 1   |
| TO      | Weighted taxonomic diversity    | Abundance          | wi = i   |

---

## 2. Core Mathematical Formulas

### 2.1 Taxonomic Distance (pTO+) — Equation 5

```
pTO+ = ln( PRODUCT_{i=1}^{L} [ wi * ( (e^{Ed_S})^2 / e^{Ed_i} + 1 ) ] )
```

where L = number of active taxonomic levels.

### 2.2 Taxonomic Diversity (pTO) — Equation 4

```
pTO = ln( SUM_{k=0}^{ns} [ (ns - nk) * PRODUCT_{i=1}^{L} [ wi * ( (e^{Ed_S})^2 / e^{Ed_i,k} + 1 ) ] ]
         / (ns + SUM(nk)) )
```

where:
- ns = max(abundance) = total number of slicing steps
- nk = threshold at step k (values: 0, 1, 2, ..., ns-1)
- Ed_S = Deng entropy at species level for slice k
- Ed_i,k = Deng entropy at level i for slice k
- wi = weight for level i (unweighted: wi=1; weighted: wi=i)

### 2.3 Weight System

| Level      | Index (i) | Weighted (wi) | Unweighted (wi) |
|------------|-----------|---------------|-----------------|
| Species    | 1         | 1             | 1               |
| Genus      | 2         | 2             | 1               |
| Family     | 3         | 3             | 1               |
| Order      | 4         | 4             | 1               |
| Class      | 5         | 5             | 1               |
| Phylum     | 6         | 6             | 1               |
| Kingdom    | 7         | 7             | 1               |

---

## 3. Run 1 — Deterministic Calculation (`ozkan_pto`)

### 3.1 Input / Output

**Input:**

| Parameter    | Type             | Description                                     |
|--------------|------------------|-------------------------------------------------|
| `community`  | named numeric    | Species abundances, names match taxonomy        |
| `tax_tree`   | data.frame       | Col 1 = Species, Col 2..n = Genus..Kingdom      |

**Output:**

| Field       | Type    | Description                                    |
|-------------|---------|------------------------------------------------|
| `uTO`       | numeric | Unweighted taxonomic diversity                 |
| `TO`        | numeric | Weighted taxonomic diversity                   |
| `uTO_plus`  | numeric | Unweighted taxonomic distance                  |
| `TO_plus`   | numeric | Weighted taxonomic distance                    |
| `Ed_levels` | named numeric | Deng entropy at each taxonomic level     |

### 3.2 Pseudocode

```
FUNCTION ozkan_pto(community, tax_tree):

    # === STEP 0: Input validation ===
    VALIDATE community is named numeric, all values >= 0
    VALIDATE tax_tree is data.frame with >= 2 columns
    REMOVE species with abundance = 0
    MATCH community species to tax_tree rows
    SET n_species = number of species with abundance > 0
    SET n_levels  = ncol(tax_tree) - 1   # taxonomic levels (excl. species col)

    IF n_species <= 1:
        RETURN all zeros (no diversity with 0-1 species)

    # === STEP 1: Compute Deng entropy at nk=0 slice (all species present) ===
    # At species level with EQUAL WEIGHTS (presence/absence):
    #   Each species gets mass = 1/S
    #   |Fi| = 1 for all species
    #   Ed_S = -SUM((1/S) * ln(1/S)) = ln(S)
    Ed_S = ln(n_species)

    # At higher levels (Genus, Family, ...):
    FOR each level L from Genus to Kingdom:
        Count unique groups at level L
        IF only 1 group:
            Ed[L] = 0
            STOP processing higher levels (termination rule)
        ELSE:
            FOR each group g at level L:
                group_count[g] = number of present species in group g
                group_size[g]  = group_count[g]  # presence-based
            Ed[L] = deng_entropy_level(group_count, group_sizes = group_size)
            MARK level L as "active"

    # === STEP 2: Compute pTO+ (taxonomic distance) from nk=0 only ===
    core_unweighted = 1.0
    core_weighted   = 1.0

    FOR each active level (index i):
        ratio = (e^Ed_S)^2 / e^Ed[i] + 1
        core_unweighted = core_unweighted * (1 * ratio)    # wi = 1
        core_weighted   = core_weighted   * (i * ratio)    # wi = i

    uTO_plus = ln(core_unweighted)
    TO_plus  = ln(core_weighted)

    # === STEP 3: Slicing procedure for pTO (taxonomic diversity) ===
    ns = max(abundance)     # total number of slicing steps
    sum_unweighted = 0.0
    sum_weighted   = 0.0

    FOR nk = 0 TO (ns - 1):
        factor = ns - nk

        # Species survive if abundance > nk
        surviving_species = species WHERE abundance > nk

        IF count(surviving_species) == 0:
            BREAK

        IF count(surviving_species) == 1:
            # Single species: Ed=0, ratio = 0+1 = 2, product = 2
            core_u_step = 2
            core_w_step = 2
        ELSE:
            # Compute Deng entropy at all levels for surviving species
            # Using EQUAL WEIGHTS (presence/absence)
            Ed_S_step = ln(count(surviving_species))
            FOR each level L (Genus to Kingdom):
                ... (same as Step 1 but for surviving species only)
            # Compute product formula (same as Step 2)
            ... compute core_u_step, core_w_step

        sum_unweighted += factor * core_u_step
        sum_weighted   += factor * core_w_step

    # Denominator: ns + SUM(0, 1, 2, ..., ns-1) = ns + ns*(ns-1)/2
    denom = ns + ns * (ns - 1) / 2

    uTO = ln(sum_unweighted / denom)
    TO  = ln(sum_weighted / denom)

    RETURN {uTO, TO, uTO_plus, TO_plus, Ed_levels}
```

### 3.3 Slicing Procedure — Worked Example

Community: 8 species with abundances `[4, 2, 3, 1, 2, 3, 2, 2]`

```
ns = max(4,2,3,1,2,3,2,2) = 4

Step 1 (nk=0): All 8 species survive (abundance > 0)
  factor = 4 - 0 = 4
  Presence: [1,1,1,1,1,1,1,1]  -> 8 species present

Step 2 (nk=1): 7 species survive (abundance > 1)
  factor = 4 - 1 = 3
  Presence: [1,1,1,0,1,1,1,1]  -> sp4 (abundance=1) excluded

Step 3 (nk=2): 3 species survive (abundance > 2)
  factor = 4 - 2 = 2
  Presence: [1,0,1,0,0,1,0,0]  -> sp1,sp3,sp6 survive

Step 4 (nk=3): 1 species survives (abundance > 3)
  factor = 4 - 3 = 1
  Presence: [1,0,0,0,0,0,0,0]  -> only sp1 survives

Denominator = 4 + (0+1+2+3) = 4 + 6 = 10

Numerator = SUM( factor * product_at_step )
          = 4 * P(8sp) + 3 * P(7sp) + 2 * P(3sp) + 1 * P(1sp)

uTO = ln(Numerator / 10)
```

---

## 4. Run 2 — Stochastic Resampling (`ozkan_pto_resample`)

### 4.1 Input / Output

**Input:**

| Parameter    | Type             | Description                               |
|--------------|------------------|-------------------------------------------|
| `community`  | named numeric    | Species abundances                        |
| `tax_tree`   | data.frame       | Taxonomic hierarchy                       |
| `n_iter`     | integer (>=101)  | Number of iterations (default: 101)       |
| `seed`       | integer or NULL  | Random seed for reproducibility           |

**Output:**

| Field              | Type       | Description                                |
|--------------------|------------|--------------------------------------------|
| `uTO_plus_max`     | numeric    | MAX(uTO+) across all iterations            |
| `TO_plus_max`      | numeric    | MAX(TO+) across all iterations             |
| `uTO_max`          | numeric    | MAX(uTO) across all iterations             |
| `TO_max`           | numeric    | MAX(TO) across all iterations              |
| `uTO_plus_det`     | numeric    | Deterministic uTO+ (iteration 1)          |
| `TO_plus_det`      | numeric    | Deterministic TO+ (iteration 1)           |
| `uTO_det`          | numeric    | Deterministic uTO (iteration 1)           |
| `TO_det`           | numeric    | Deterministic TO (iteration 1)            |
| `n_iter`           | integer    | Number of iterations performed             |
| `iteration_results`| data.frame | All 4 components for each iteration        |

### 4.2 Pseudocode

```
FUNCTION ozkan_pto_resample(community, tax_tree, n_iter = 101, seed = NULL):

    # === STEP 0: Validation ===
    VALIDATE n_iter >= 101
    VALIDATE at least 2 species with positive abundance
    IF seed is provided: SET random seed

    n_species = count(species with abundance > 0)

    # Pre-allocate results matrix: n_iter rows x 4 columns
    results[n_iter, 4]  # columns: uTO+, TO+, uTO, TO

    # === STEP 1: Iteration 1 — Deterministic (original community) ===
    results[1, ] = ozkan_pto(community, tax_tree)

    # === STEP 2: Iterations 2..n_iter — Stochastic resampling ===
    FOR iter = 2 TO n_iter:

        # For each species: 50% chance of inclusion
        # Excel formula: IF(H2 > 0, RANDBETWEEN(0,1) * H2, 0)
        FOR each species s:
            coin = RANDOM_INTEGER(0, 1)   # 0 or 1 with equal probability
            resampled[s] = community[s] * coin

        # Remove excluded species (abundance = 0)
        surviving = resampled WHERE resampled > 0

        IF count(surviving) < 2:
            results[iter, ] = [0, 0, 0, 0]
            CONTINUE

        # Compute pTO for resampled community
        results[iter, ] = ozkan_pto(surviving, tax_tree)

    # === STEP 3: Compute maximums ===
    uTO_plus_max = MAX(results[, "uTO_plus"])
    TO_plus_max  = MAX(results[, "TO_plus"])
    uTO_max      = MAX(results[, "uTO"])
    TO_max       = MAX(results[, "TO"])

    RETURN {
        *_max values,
        *_det values (from iteration 1),
        n_iter,
        iteration_results (full data.frame)
    }
```

### 4.3 Resampling Logic — Diagram

```
Original community:    [4, 2, 3, 1, 2, 3, 2, 2]
                        |  |  |  |  |  |  |  |
Random coin (50/50):   [1, 0, 1, 1, 0, 1, 0, 1]
                        |  |  |  |  |  |  |  |
Resampled:             [4, 0, 3, 1, 0, 3, 0, 2]
                        |     |  |     |     |
After removing zeros:  [4,    3, 1,    3,    2]  -> 5 species survive
                        |
                    ozkan_pto() -> {uTO, TO, uTO+, TO+}
```

This is repeated n_iter times. MAX of each component across all iterations
is the final Run 2 result.

---

## 5. Run 3 — Sensitivity Analysis (`ozkan_pto_sensitivity`)

### 5.1 Input / Output

**Input:**

| Parameter      | Type             | Description                               |
|----------------|------------------|-------------------------------------------|
| `community`    | named numeric    | Species abundances (original)             |
| `tax_tree`     | data.frame       | Taxonomic hierarchy                       |
| `run2_result`  | list             | Output of `ozkan_pto_resample()`          |
| `n_iter`       | integer or NULL  | Iterations (default: same as Run 2)       |
| `seed`         | integer or NULL  | Random seed                               |

**Output:**

| Field                | Type       | Description                              |
|----------------------|------------|------------------------------------------|
| `uTO_plus_max`       | numeric    | MAX(uTO+) across Run 1, 2, AND 3        |
| `TO_plus_max`        | numeric    | MAX(TO+) across Run 1, 2, AND 3         |
| `uTO_max`            | numeric    | MAX(uTO) across Run 1, 2, AND 3         |
| `TO_max`             | numeric    | MAX(TO) across Run 1, 2, AND 3          |
| `run3_uTO_plus_max`  | numeric    | MAX(uTO+) from Run 3 only               |
| `run3_TO_plus_max`   | numeric    | MAX(TO+) from Run 3 only                |
| `n_iter`             | integer    | Number of iterations                     |
| `species_probs`      | named numeric | Inclusion probability per species     |
| `iteration_results`  | data.frame | All Run 3 iteration results              |

### 5.2 Pseudocode

```
FUNCTION ozkan_pto_sensitivity(community, tax_tree, run2_result,
                                n_iter = NULL, seed = NULL):

    # === STEP 0: Validation ===
    VALIDATE run2_result contains iteration_results (from Run 2)
    IF n_iter is NULL: n_iter = run2_result$n_iter
    IF seed is provided: SET random seed

    n_species = count(species with abundance > 0)
    S = n_species

    # === STEP 1: Determine species-specific inclusion probabilities ===
    # Excel VBA formula (Module8):
    #   IF(AA2=0,
    #     IF(RANDBETWEEN(1, K13) > 1, H2, 0),
    #     IF(L25 >= RANDBETWEEN(0, K22), H2, 0))
    #
    # Where:
    #   AA2 = selection count for species (from Run 2)
    #   K13 = total number of species (S)
    #   K22 = total iterations
    #   L25 = selection score
    #
    # Translation:
    #   Unselected species (AA=0): P(include) = (S-1)/S
    #   Selected species (AA>0):   P(include) = L25/K22

    FOR each species s:
        IF selection_count[s] == 0:
            species_probs[s] = (S - 1) / S
        ELSE:
            species_probs[s] = selection_score / total_iterations

    # === STEP 2: Iteration 1 — Deterministic ===
    results[1, ] = ozkan_pto(community, tax_tree)

    # === STEP 3: Iterations 2..n_iter — Probability-weighted resampling ===
    FOR iter = 2 TO n_iter:

        FOR each species s:
            random_val = UNIFORM(0, 1)
            IF random_val < species_probs[s]:
                resampled[s] = community[s]     # include with original abundance
            ELSE:
                resampled[s] = 0                # exclude

        surviving = resampled WHERE resampled > 0

        IF count(surviving) < 2:
            results[iter, ] = [0, 0, 0, 0]
            CONTINUE

        results[iter, ] = ozkan_pto(surviving, tax_tree)

    # === STEP 4: Compute maximums ===
    # Run 3 maximums
    run3_max = column_max(results)

    # Overall maximums: MAX across Run 2 max AND Run 3 max
    overall_uTO_plus_max = MAX(run2_result$uTO_plus_max, run3_max["uTO_plus"])
    overall_TO_plus_max  = MAX(run2_result$TO_plus_max,  run3_max["TO_plus"])
    overall_uTO_max      = MAX(run2_result$uTO_max,      run3_max["uTO"])
    overall_TO_max       = MAX(run2_result$TO_max,        run3_max["TO"])

    RETURN {
        overall *_max values,
        run3-specific *_max values,
        n_iter,
        species_probs,
        iteration_results
    }
```

### 5.3 Difference Between Run 2 and Run 3

```
Run 2: All species have EQUAL inclusion probability (50%)
       coin = RANDBETWEEN(0, 1)  ->  P = 0.5

Run 3: Species have DIFFERENT inclusion probabilities based on Run 2 results
       Unselected species: P = (S-1)/S  (close to 1 for large S)
       Selected species:   P = score/iterations (varies by species)
```

The rationale is that Run 3 explores the parameter space more efficiently
by favoring the species composition that produced high diversity in Run 2.

---

## 6. Full Pipeline Flow

```
                    ┌─────────────────────────────────┐
                    │       Original Community        │
                    │   [4, 2, 3, 1, 2, 3, 2, 2]     │
                    └─────────┬───────────────────────┘
                              │
                              v
                    ┌─────────────────────────────────┐
                    │  Run 1: ozkan_pto()              │
                    │  Deterministic calculation       │
                    │  -> uTO, TO, uTO+, TO+           │
                    └─────────┬───────────────────────┘
                              │
                              v
                    ┌─────────────────────────────────┐
                    │  Run 2: ozkan_pto_resample()     │
                    │  n_iter iterations (>= 101)      │
                    │  50% random inclusion per species│
                    │  -> MAX of each component        │
                    └─────────┬───────────────────────┘
                              │
                              │ passes run2_result
                              v
                    ┌─────────────────────────────────┐
                    │  Run 3: ozkan_pto_sensitivity()  │
                    │  n_iter iterations               │
                    │  Species-specific probabilities  │
                    │  -> MAX across Run 1+2+3         │
                    └─────────────────────────────────┘
                              │
                              v
                    ┌─────────────────────────────────┐
                    │        FINAL RESULT              │
                    │  Overall MAX of uTO, TO, uTO+,  │
                    │  TO+ across all three runs       │
                    └─────────────────────────────────┘
```

### Ordering Property

For any valid community with >1 taxonomic level:

```
Run 3 overall MAX >= Run 2 MAX >= Run 1 deterministic
```

This holds because each subsequent run takes the pmax (parallel maximum)
with the previous run's results. The final result is guaranteed to be at
least as large as the deterministic calculation.

---

## 7. Edge Cases and Termination Rules

| Condition                                | Behavior                       |
|------------------------------------------|--------------------------------|
| Only 1 species with abundance > 0        | All components = 0             |
| All species have equal abundance         | Slicing has only 1 step        |
| Only 1 node at a taxonomic level         | Ed = 0, computation stops      |
| All species excluded in a stochastic run | That iteration returns 0       |
| Only 1 species survives stochastic run   | That iteration returns 0       |
| n_iter < 101                             | ERROR: minimum 101 required    |

---

## 8. Relationship to Other Functions

```
deng_entropy_level()          # Core: Deng entropy at a single level
    ^
    |  called at each level for each slice
    |
ozkan_pto()                   # Run 1: deterministic pTO
    ^
    |  called at each iteration
    |
ozkan_pto_resample()          # Run 2: stochastic resampling
    |
    |  output feeds into
    v
ozkan_pto_sensitivity()       # Run 3: sensitivity analysis
    |
    v
pto_components()              # Convenience: named vector wrapper for Run 1
```

---

## 9. Excel VBA Macro Correspondence

| R Function                 | VBA Module   | VBA Procedure  | Button    |
|----------------------------|-------------|----------------|-----------|
| `ozkan_pto()`              | Module1     | `tekerur()`    | İşlem 1   |
| `ozkan_pto_resample()`     | Module3     | `tekerur2()`   | İşlem 2   |
| `ozkan_pto_sensitivity()`  | Module8     | `tekerur3()`   | İşlem 3   |

**Known Excel bug:** Both `tekerur2()` and `tekerur3()` lack
`Application.Calculate` inside their iteration loops. This causes
`RANDBETWEEN()` formulas not to regenerate between iterations, resulting
in identical random values for all iterations within a single session.
The R implementation correctly generates new random values at each iteration.

---

## References

- Deng, Y. (2016). Deng entropy. *Chaos, Solitons & Fractals*, 91, 549-553.
- Özkan, K. (2018). Taksonomik çeşitliliğin belirlenmesi için yeni önerilen
  bir eşitlik. *Turkish Journal of Forestry*, 19(4), 336-346.
  DOI: 10.18182/tjf.441061
- Özkan, K., Mert, A., Şenol, A., Özdemir, S. (2018). Macrotakdivozkan
  Excel macro. http://www.kantitatifekoloji.net/takdivozkan
