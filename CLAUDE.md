# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build, Test & Check Commands

```bash
# Run all 610+ tests
cd /tmp/taxdiv-repo && Rscript -e 'devtools::test()'

# Run a single test file
Rscript -e 'testthat::test_file("tests/testthat/test-ozkan_pto.R")'

# R CMD check (standard)
R CMD build . && R CMD check taxdiv_0.1.0.tar.gz

# R CMD check with CRAN-level strictness
R CMD build . && R CMD check --as-cran taxdiv_0.1.0.tar.gz

# Regenerate documentation from roxygen2 comments
Rscript -e 'devtools::document()'

# Spell check (uses inst/WORDLIST for allowed words)
Rscript -e 'spelling::spell_check_package()'

# Update WORDLIST after adding new technical terms
Rscript -e 'spelling::update_wordlist()'

# Build pkgdown site locally
Rscript -e 'pkgdown::build_site()'
```

## Architecture

**taxdiv** calculates taxonomic diversity indices for ecological communities using three approaches:

1. **Classical indices** (`classical_indices.R`): Shannon, Simpson — species-level diversity
2. **Clarke & Warwick** (`classical_indices.R`): Delta, Delta*, AvTD, VarTD — taxonomy-aware distinctness
3. **Ozkan pTO** (`ozkan_pto.R`): Deng entropy-based proportional taxonomic originality — the package's core contribution

The Ozkan pTO pipeline has three sequential runs:
- **Run 1** (`ozkan_pto.R`): Deterministic calculation of 8 indices (TO, TO+, uTO, pTO, etc.)
- **Run 2** (`ozkan_pto_resample.R`): Stochastic bootstrap resampling with confidence intervals
- **Run 3** (`ozkan_pto_sensitivity.R`): Jackknife sensitivity analysis (species removal impact)

`ozkan_pto_full()` chains all three runs together.

### Key data flow

```
community vector + tax_tree → ozkan_pto() → S3 object
                             → ozkan_pto_resample() → S3 object
                             → ozkan_pto_sensitivity() → S3 object
                             → plot/summary/print methods
```

### Input conventions

- **community**: Named numeric vector — species names as names, abundances as values
- **tax_tree**: Data frame built with `build_tax_tree()` — columns are Species, Genus, Family, Order, ... (lowest to highest rank). Prefer `build_tax_tree()` over raw `data.frame()`.
- **weights**: Numeric vector for taxonomic levels — defaults to equal step weights (1, 2, 3, ...)

### S3 class system

6 S3 classes (`s3_methods.R`) with print/summary/plot methods: `ozkan_pto`, `ozkan_pto_resample`, `ozkan_pto_sensitivity`, `compare_indices`, `batch_analysis`, `rarefaction_taxonomic`, `td_simulation`.

## Key Conventions

- Functions use a `seed` parameter for reproducibility — never use `set.seed()` directly
- All function names are `snake_case`; plot functions prefixed with `plot_`
- Minimal dependencies: only `stats` and `rlang` (ggplot2 is Suggests-only)
- Roxygen2 with markdown enabled — LaTeX formulas use `\deqn{}` and `\eqn{}`
- Tests validate mathematical constraints (e.g., `TO+ >= TO`, Shannon equivalence at species level)
- Turkish vignette (`giris_rehberi.Rmd`) is intentionally kept — all other content is English only
- `inst/WORDLIST` contains ~723 words (mostly Turkish from the vignette) for `spelling::spell_check_package()`

## Datasets

- `gazi_comm` + `gazi_gytk`: Single-community quick demos (8 species, 4 taxonomic levels)
- `anatolian_trees`: Multi-site dataset (20 species, 3 sites, 7 taxonomic levels, Westhoff-Maarel scale)
