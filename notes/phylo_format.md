# The `ape::phylo` Format and Its Relevance to taxdiv

> Research note for the taxdiv package project.
> Reviewed: February 2026

---

## 1. What Is the `phylo` Class?

The `phylo` class, defined in the `ape` package (Paradis & Schliep, 2019),
is the de facto standard for representing phylogenetic trees in R. Almost
every phylogenetic/evolutionary ecology package in the R ecosystem depends
on it (picante, hillR, entropart, phytools, ggtree, etc.).

A `phylo` object is a list with the following components:

| Component    | Type              | Description                              |
|-------------|-------------------|------------------------------------------|
| `edge`      | integer matrix    | 2-column matrix: [parent, child] pairs   |
| `edge.length`| numeric vector   | Branch lengths for each edge (optional)  |
| `tip.label` | character vector  | Names of terminal nodes (species)        |
| `Nnode`     | integer           | Number of internal nodes                 |
| `node.label`| character vector  | Names of internal nodes (optional)       |
| `root.edge` | numeric           | Length of root edge (optional)           |

### Minimal Example

A tree with 4 species and 3 internal nodes:

```
          root(5)
         /       \
      (6)         (7)
     /   \       /   \
   sp1   sp2   sp3   sp4
```

```r
tree <- list(
  edge = matrix(c(
    5, 6,   # root -> node6
    5, 7,   # root -> node7
    6, 1,   # node6 -> sp1
    6, 2,   # node6 -> sp2
    7, 3,   # node7 -> sp3
    7, 4    # node7 -> sp4
  ), ncol = 2, byrow = TRUE),
  tip.label = c("sp1", "sp2", "sp3", "sp4"),
  Nnode = 3L,
  edge.length = c(2, 2, 1, 1, 1, 1)  # optional
)
class(tree) <- "phylo"
```

---

## 2. Converting Between Classification Tables and `phylo`

### 2.1 Classification Table → `phylo` (via `ape`)

The `ape` package provides `as.phylo.formula()`:

```r
library(ape)

tax <- data.frame(
  Species = c("sp1", "sp2", "sp3", "sp4", "sp5"),
  Genus   = c("G1",  "G1",  "G2",  "G2",  "G3"),
  Family  = c("F1",  "F1",  "F1",  "F2",  "F2"),
  Order   = c("O1",  "O1",  "O1",  "O1",  "O1")
)

# Formula reads right-to-left: ~Order/Family/Genus/Species
tree <- as.phylo(~Order/Family/Genus/Species, data = tax)
plot(tree)
```

**Properties of the resulting tree:**
- Ultrametric (all tips at same depth)
- Equal branch lengths within each level
- Internal nodes labeled with taxonomic group names
- No polytomies — strictly bifurcating only if each group has exactly 2 children

**Limitations:**
- Branch lengths are arbitrary (based on rank distance, not evolutionary time)
- `varstep = TRUE` in `taxa2dist()` adjusts for unequal subtree sizes but
  this is Clarke & Warwick's convention, not biological data

### 2.2 Classification Table → Distance Matrix (via `vegan`)

The `vegan` package provides `taxa2dist()`:

```r
library(vegan)

# tax_table: rows = species, columns = ranks (Genus, Family, Order, ...)
# NOTE: species column should NOT be included
tax_matrix <- tax[, -1]  # remove Species column
rownames(tax_matrix) <- tax$Species

d <- taxa2dist(tax_matrix, varstep = FALSE)
# Returns a dist object with pairwise taxonomic distances
```

This is used internally by `vegan::taxondive()` for Clarke & Warwick indices.

### 2.3 `phylo` → Distance Matrix

```r
library(ape)
d <- cophenetic.phylo(tree)
# Returns a symmetric matrix of pairwise tip distances
```

### 2.4 `taxize` Package Bridge

The `taxize` package (Chamberlain & Szocs, 2013) provides `class2tree()`:

```r
library(taxize)
# Retrieve classification from online databases
cls <- classification(c("Quercus robur", "Pinus sylvestris"), db = "ncbi")

# Convert to phylo + distance matrix
tree_result <- class2tree(cls)
tree_result$phylo       # phylo object
tree_result$distmat     # distance matrix
```

---

## 3. How Other R Packages Use `phylo`

| Package   | Requires `phylo`? | How taxonomy is handled              |
|-----------|-------------------|--------------------------------------|
| vegan     | No                | `taxa2dist()` → dist, OR data.frame  |
| hillR     | Yes (for phylo)   | `hill_phylo(comm, tree)` — `phylo`   |
| entropart | Yes               | `phylo`, `hclust`, or `phylog`       |
| picante   | Yes               | `pd(samp, tree)` — always `phylo`    |
| ape       | Defines it        | Foundational `phylo` class           |
| **taxdiv**| **No**            | **Classification table (data.frame)**|

---

## 4. taxdiv's Design Choice: Classification Tables

The taxdiv package deliberately does NOT depend on `ape` or use `phylo`
objects. Instead, it works directly with classification tables (data frames).

### 4.1 Advantages of This Approach

1. **No external dependency:** taxdiv has zero Imports beyond base R stats.
   Adding `ape` would pull in a significant dependency tree.

2. **Familiar data format:** Ecologists typically record taxonomy as tables
   (Excel spreadsheets, CSV files), not as Newick strings or edge matrices.
   No conversion step needed.

3. **Exact level information preserved:** A data frame column directly
   represents a Linnean rank (Genus, Family, Order). With `phylo`, rank
   information is lost — only branching topology and distances remain.

4. **Deng entropy requires level-aware computation:** The pTO formula
   computes Ed separately at each named taxonomic level (Species, Genus,
   Family, ...) with level-specific weights (wi = i). This requires knowing
   WHICH level each node belongs to — information that `phylo` does not
   natively store.

5. **Faster for small-to-medium datasets:** Direct data frame iteration
   avoids tree traversal overhead. For the typical ecological survey
   (10–200 species, 4–7 taxonomic levels), this is negligible but keeps
   the code simpler.

### 4.2 Disadvantages

1. **No interoperability:** Users cannot easily pipe taxdiv results into
   `picante::ses.mpd()` or `hillR::hill_phylo()` without manual conversion.

2. **No continuous branch lengths:** Classification tables assume equal
   (or linearly increasing) distances between ranks. True phylogenetic
   distances from molecular data cannot be represented.

3. **Limited visualization:** `ape::plot.phylo()` provides sophisticated
   tree plotting. taxdiv would need its own visualization or conversion
   utilities.

---

## 5. Potential Future Integration

If taxdiv were to support `phylo` objects in the future, the most practical
approach would be:

### Option A: Accept Both Formats (Recommended)

```r
# User can provide either:
ozkan_pto(comm, tax_tree)            # data.frame (current)
ozkan_pto(comm, phylo_tree)          # phylo object (future)

# Internal dispatch based on class:
if (inherits(tax_tree, "phylo")) {
  # Convert phylo -> classification table using node.label
  # Then proceed as normal
} else if (is.data.frame(tax_tree)) {
  # Current implementation
}
```

**Requirements:**
- Add `ape` to `Suggests:` (not `Imports:`)
- Write a `phylo_to_tax_table()` conversion utility
- Node labels in the `phylo` object must encode taxonomic rank information

### Option B: Provide Conversion Utilities Only

```r
# User converts explicitly:
tax_df <- phylo_to_tax_table(tree, ranks = c("Genus", "Family", "Order"))
ozkan_pto(comm, tax_df)

# And the reverse:
tree <- tax_table_to_phylo(tax_df)
```

This keeps the core functions unchanged and avoids conditional logic.

### Option C: Bridge to vegan's `taxondive`

```r
# For Clarke & Warwick comparison:
d <- tax_distance_matrix(tax_tree)        # taxdiv's output
vegan_result <- taxondive(comm, d)        # vegan's implementation
taxdiv_result <- delta(comm, tax_tree)    # taxdiv's implementation
```

This already works since `tax_distance_matrix()` returns a `dist` object
compatible with vegan.

---

## 6. The `phylo` Format in Detail

### 6.1 Newick String Representation

Trees can be read from Newick (parenthetical) format:

```r
tree <- ape::read.tree(text = "((sp1:1,sp2:1):1,(sp3:1,sp4:1):1);")
```

Newick string anatomy:
```
((sp1:1,sp2:1)Genus1:1,(sp3:1,sp4:1)Genus2:1)Family1;
  ^tip  ^len  ^internal ^len
```

### 6.2 Edge Matrix Encoding

The `edge` matrix is the core of the `phylo` object:

```
     [,1] [,2]
[1,]    5    6    # root(5) -> internal(6)
[2,]    5    7    # root(5) -> internal(7)
[3,]    6    1    # internal(6) -> tip(1) = sp1
[4,]    6    2    # internal(6) -> tip(2) = sp2
[5,]    7    3    # internal(7) -> tip(3) = sp3
[6,]    7    4    # internal(7) -> tip(4) = sp4
```

**Numbering convention:**
- Tips are numbered 1..n (matching `tip.label` order)
- Internal nodes are numbered (n+1)..(n+Nnode)
- Root is always node (n+1)

### 6.3 Ultrametric Trees and Taxonomy

A tree is ultrametric if all root-to-tip path lengths are equal. Trees
derived from classification tables are always ultrametric because each
taxonomic rank represents a fixed distance from the root.

```r
is.ultrametric(tree)  # Should be TRUE for taxonomy-derived trees
```

Non-ultrametric trees (from molecular phylogenetics) cannot be directly
used with methods assuming equal-rank distances (like Clarke & Warwick).

---

## 7. Key `ape` Functions for Taxonomy

| Function                          | Purpose                              |
|-----------------------------------|--------------------------------------|
| `as.phylo.formula(~A/B/C, data)` | Classification table → phylo         |
| `cophenetic.phylo(tree)`          | phylo → pairwise distance matrix     |
| `plot.phylo(tree)`                | Visualize tree                       |
| `is.ultrametric(tree)`            | Check equal root-to-tip distances    |
| `read.tree(file)`                 | Read Newick format                   |
| `write.tree(tree, file)`          | Write Newick format                  |
| `Ntip(tree)` / `Nnode(tree)`      | Count tips / internal nodes          |
| `drop.tip(tree, tip)`             | Remove species from tree             |
| `extract.clade(tree, node)`       | Extract subtree                      |

---

## References

- Paradis, E. & Schliep, K. (2019). ape 5.0: an environment for modern
  phylogenetics and evolutionary analyses in R. *Bioinformatics*, 35, 526-528.
- Chamberlain, S. & Szocs, E. (2013). taxize: taxonomic search and retrieval
  in R. *F1000Research*, 2, 191.
- Felsenstein, J. (2004). *Inferring Phylogenies*. Sinauer Associates.
