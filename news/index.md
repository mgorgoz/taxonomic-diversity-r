# Changelog

## taxdiv 0.1.0

### New Features

- Core Deng entropy calculation at any taxonomic level
  ([`deng_entropy_level()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/deng_entropy_level.md))
- Ozkan (2018) pTO formula with slicing procedure
  ([`ozkan_pto()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/ozkan_pto.md),
  [`pto_components()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/pto_components.md))
  - Unweighted/weighted taxonomic diversity (uTO, TO)
  - Unweighted/weighted taxonomic distance (uTO+, TO+)
  - Presence-based entropy: species receive equal weight (1/S) at each
    slice
- Stochastic resampling — Run 2
  ([`ozkan_pto_resample()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/ozkan_pto_resample.md))
  - Random species inclusion/exclusion (50% probability per species)
  - Configurable iteration count and seed for reproducibility
  - Returns deterministic (Run 1) + maximum values across iterations
- Sensitivity analysis — Run 3
  ([`ozkan_pto_sensitivity()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/ozkan_pto_sensitivity.md))
  - Species-specific inclusion probabilities derived from Run 2 results
  - Overall maximum across Run 1 + 2 + 3
  - Species probability table output
- Classical diversity indices: Shannon
  ([`shannon()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/shannon.md))
  and Simpson
  ([`simpson()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/simpson.md))
- Clarke & Warwick taxonomic distinctness measures:
  - Taxonomic diversity Delta
    ([`delta()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/delta.md))
  - Taxonomic distinctness Delta\*
    ([`delta_star()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/delta_star.md))
  - Average taxonomic distinctness AvTD/Delta+
    ([`avtd()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/avtd.md))
  - Variation in taxonomic distinctness VarTD/Lambda+
    ([`vartd()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/vartd.md))
- Taxonomic distance matrix computation
  ([`tax_distance_matrix()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/tax_distance_matrix.md))
- Taxonomy tree builder
  ([`build_tax_tree()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/build_tax_tree.md))
- Introductory vignette with Mediterranean forest example

### Documentation

- Pseudocode for pTO Run 1/2/3 pipeline
  (`inst/pseudocode/taxonomic_diversity_pseudocode.md`)
- Pseudocode for Clarke & Warwick measures
  (`inst/pseudocode/taxonomic_distance_pseudocode.md`)
- Research note on ape phylo format (`notes/phylo_format.md`)
- Comparative review of existing R packages
  (`notes/existing_packages.md`)
- Turkish (Turkce) inline comments added to all test files
- Vignette expanded with Run 2 (stochastic resampling), Run 3
  (sensitivity analysis), species inclusion probabilities, and full
  pipeline summary sections

### Testing

- 218 assertions across 65 test blocks in 6 test files
- Integration tests covering full Run 1 → 2 → 3 pipeline
  (`test-integration.R`)
- Reproducibility tests with seed control
- Mathematical relationship checks (ordering constraints, boundary
  cases)

### Validation

- Verified against Ozkan (2018) hypothetical examples (Community A =
  4.000034, Community B = 3.555348)
- Cross-validated with Kursad Ozkan’s Excel macro (TD_OMD.xlsm)
- Documented Excel macro bug: `Application.Calculate` missing inside
  tekerur2() and tekerur3() loops
