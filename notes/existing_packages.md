# Existing R Packages for Biodiversity / Taxonomic Diversity

> Research note for the taxdiv package project.
> Reviewed: February 2026

---

## 1. vegan — Community Ecology Package

**Purpose:** Most widely used R package for community ecology. Ordination,
diversity analysis, ecological null models.

**Key diversity functions:**

| Function | Purpose |
|----------|---------|
| `diversity(x, index)` | Shannon, Simpson, inverse-Simpson |
| `renyi()` / `tsallis()` | Renyi/Tsallis generalized entropy, Hill numbers |
| `taxondive(comm, dis)` | Clarke & Warwick Delta, Delta*, Delta+, Lambda+ |
| `taxa2dist(x, varstep)` | Classification table to taxonomic distance matrix |
| `rarefy()`, `specnumber()` | Rarefaction, species richness |

**Taxonomic hierarchy:** `taxondive()` implements Clarke & Warwick (1998, 2001).
Requires a pre-computed `dist` object from `taxa2dist()`. Treats taxonomy
purely as a distance matrix — no notion of hierarchical levels or Deng entropy.

**API design:** S3. Functions accept matrices/data.frames, return named vectors
or lightweight S3 objects with `plot()` and `summary()`.

**Strengths:** Mature, widely cited, industry standard. Only mainstream
implementation of Clarke & Warwick indices.

**Limitations:** No Deng entropy, no pTO indices (Ozkan 2018). No bias
correction for taxonomic indices. API somewhat inconsistent across functions.

---

## 2. hillR — Diversity Through Hill Numbers

**Purpose:** Taxonomic, functional, and phylogenetic alpha/beta diversity
through the Hill numbers framework (Chao, Chiu & Jost 2014).

**Key functions:**

| Function | Purpose |
|----------|---------|
| `hill_taxa(comm, q)` | Taxonomic Hill number per site |
| `hill_taxa_parti()` | Gamma/alpha/beta partitioning |
| `hill_func(comm, traits, q)` | Functional diversity |
| `hill_phylo(comm, tree, q)` | Phylogenetic Hill diversity |

**Taxonomic hierarchy:** `hill_taxa()` does NOT use taxonomy at all — treats
species as equivalent units. Phylogenetic structure handled only via `phylo`
tree objects (from `ape`), not classification tables.

**API design:** Purely functional, flat. Consistent naming:
`hill_{facet}` / `hill_{facet}_parti` / `hill_{facet}_parti_pairwise`.

**Strengths:** Clean API, unified Hill numbers across three facets, lightweight.

**Limitations:** No taxonomic hierarchy for the taxonomic facet.
No Clarke & Warwick, no Deng entropy, no pTO.

---

## 3. entropart — Entropy Partitioning to Measure Diversity

**Purpose:** Rigorous entropy-based diversity with state-of-the-art bias
correction. Supports species-neutral, phylogenetic, similarity-based entropy.

**Key functions:**

| Function | Purpose |
|----------|---------|
| `Shannon()` / `bcShannon()` | Shannon entropy + bias-corrected |
| `Simpson()` / `bcSimpson()` | Simpson + bias-corrected |
| `Diversity()` / `bcDiversity()` | Hill numbers |
| `Hqz()` / `bcHqz()` | Similarity-based entropy (Z matrix) |
| `PhyloEntropy()` | Phylogenetic entropy via tree slicing |
| `DivPart()` / `DivEst()` | Partitioning + bootstrap CI |
| `CommunityProfile()` | Entropy/diversity profile across q |

**Taxonomic hierarchy:** Handled through `phylo`, `hclust`, or `phylog` tree
objects. No direct support for Linnean classification tables — must convert
to tree first (e.g., `ape::as.phylo.formula()`).

**API design:** S3 with key classes (`MetaCommunity`, `AbdVector`,
`ProbaVector`). Auto-selects bias correction based on input type.

**Strengths:** Most rigorous entropy framework in R. State-of-the-art bias
correction (Chao & Jost estimators). Bootstrap confidence intervals.

**Limitations:** No Clarke & Warwick. No Deng entropy/pTO. Taxonomy must be
encoded as ultrametric tree. Heavier dependency footprint. Steeper learning curve.

---

## 4. Other Relevant Packages

### ape — Analyses of Phylogenetics and Evolution

Foundational package defining the `phylo` class. Key bridge function:
`as.phylo.formula(~Order/Family/Genus/Species, data=df)` converts a Linnean
classification table into a `phylo` tree. Does not compute diversity indices.

### picante — Integrating Phylogenies and Ecology

Phylogenetic community ecology: `pd()` (Faith's PD), `mpd()`, `mntd()`,
`ses.mpd()` (null model significance). Requires `phylo` tree objects.
No Clarke & Warwick, no entropy-based taxonomic measures.

### taxize — Taxonomic Search and Retrieval (rOpenSci)

Programmatic access to 13+ taxonomic databases (NCBI, GBIF, ITIS).
`classification()` retrieves full Linnean hierarchy. `class2tree()` converts
to `phylo` + distance matrix. Data retrieval only — no diversity calculations.

### abdiv — Alpha and Beta Diversity Measures

Comprehensive collection of alpha/beta indices in flat functional API.
No taxonomic hierarchy support.

### metacoder — Hierarchical Taxonomic Visualization

"Heat tree" visualization of hierarchical taxonomic data from metabarcoding.
Not for computing diversity indices.

---

## 5. Comparative Table

| Feature | vegan | hillR | entropart | picante | **taxdiv** |
|---------|-------|-------|-----------|---------|------------|
| Shannon / Simpson / Hill | Yes | Yes | Yes | No | Yes |
| Clarke & Warwick (Delta, Delta+) | Yes | No | No | No | **Yes** |
| Deng entropy (Ed) | No | No | No | No | **Yes** |
| Ozkan pTO (TO, uTO, TO+, uTO+) | No | No | No | No | **Yes** |
| Classification table input | Yes | No | No | No | **Yes** |
| Stochastic resampling (Run 2/3) | No | No | No | No | **Yes** |
| Bias correction | No | No | Yes | No | No |
| Null model significance | Limited | No | Bootstrap | Yes | No |
| Phylogenetic diversity | No | Yes | Yes | Yes | No |
| Functional diversity | No | Yes | Yes | No | No |

---

## 6. Unique Position of taxdiv

1. **Only R package implementing Deng entropy-based taxonomic diversity**
   (Ozkan 2018 pTO indices). No existing CRAN package covers this.

2. **Direct classification table input** — No tree conversion needed.
   Other packages (hillR, entropart, picante) require `phylo` objects.

3. **Combines Clarke & Warwick + Deng entropy** in a single package,
   enabling direct comparison between traditional and novel approaches.

4. **Stochastic resampling** (Run 2/3) for sensitivity analysis —
   unique feature not available in any other package.

5. **Potential future additions:**
   - Bias correction (Chao-Shen, Jackknife) like entropart's `bc*` functions
   - Null model significance testing like picante's `ses.*` functions
   - Hill number unification across facets like hillR
   - taxize integration for automatic taxonomy retrieval

---

## References

- Oksanen, J. et al. (2022). vegan: Community Ecology Package. CRAN.
- Li, D. (2018). hillR: taxonomic, functional, and phylogenetic diversity
  and similarity through Hill Numbers. JOSS, 3(31), 1041.
- Marcon, E. & Herault, B. (2015). entropart: An R Package to Measure and
  Partition Diversity. JSS, 67(8).
- Paradis, E. & Schliep, K. (2019). ape 5.0: an environment for modern
  phylogenetics and evolutionary analyses in R. Bioinformatics, 35, 526-528.
- Kembel, S.W. et al. (2010). Picante: R tools for integrating phylogenies
  and ecology. Bioinformatics, 26, 1463-1464.
- Chamberlain, S. & Szocs, E. (2013). taxize: taxonomic search and retrieval
  in R. F1000Research, 2, 191.
