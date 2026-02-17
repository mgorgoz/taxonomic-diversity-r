# taxdiv 0.1.0

## New Features

* Core Deng entropy calculation at any taxonomic level (`deng_entropy_level()`)
* Ozkan (2018) pTO formula with slicing procedure (`ozkan_pto()`, `pto_components()`)
  - Unweighted/weighted taxonomic diversity (uTO, TO)
  - Unweighted/weighted taxonomic distance (uTO+, TO+)
  - Presence-based entropy: species receive equal weight (1/S) at each slice
* Classical diversity indices: Shannon (`shannon()`) and Simpson (`simpson()`)
* Clarke & Warwick taxonomic distinctness measures:
  - Taxonomic diversity Delta (`delta()`)
  - Taxonomic distinctness Delta* (`delta_star()`)
  - Average taxonomic distinctness AvTD/Delta+ (`avtd()`)
  - Variation in taxonomic distinctness VarTD/Lambda+ (`vartd()`)
* Taxonomic distance matrix computation (`tax_distance_matrix()`)
* Taxonomy tree builder (`build_tax_tree()`)
* Introductory vignette with Mediterranean forest example

## Validation

* Verified against Ozkan (2018) hypothetical examples (Community A = 4.000034, Community B = 3.555348)
* Cross-validated with Kursad Ozkan's Excel macro (TD_OMD.xlsm)
* 75 unit tests covering all exported functions
