# Package index

## Classical Diversity Indices

Standard diversity indices widely used in ecology. Shannon-Wiener
(Shannon, 1948) and Simpson (Simpson, 1949) indices.

- [`shannon()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/shannon.md)
  : Shannon Diversity Index
- [`simpson()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/simpson.md)
  : Simpson Diversity Index

## Clarke & Warwick Taxonomic Distinctness

Path-length-based taxonomic distinctness measures following Warwick &
Clarke (1995) and Clarke & Warwick (1998, 1999, 2001).

- [`delta()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/delta.md)
  : Taxonomic Diversity Index (Delta)
- [`delta_star()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/delta_star.md)
  : Taxonomic Distinctness (Delta\*)
- [`avtd()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/avtd.md)
  : Average Taxonomic Distinctness (Delta+)
- [`vartd()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/vartd.md)
  : Variation in Taxonomic Distinctness (Lambda+)

## Ozkan pTO Method

Deng entropy-based taxonomic diversity indices (Ozkan, 2018; Deng, 2016;
Dempster, 1967; Shafer, 1976).

- [`ozkan_pto()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/ozkan_pto.md)
  : Calculate Ozkan's Taxonomic Diversity Index (pTO)
- [`pto_components()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/pto_components.md)
  : Calculate All Eight pTO Components (Convenience Wrapper)
- [`deng_entropy_level()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/deng_entropy_level.md)
  : Calculate Deng Entropy at a Single Taxonomic Level

## Ozkan pTO Pipeline

Stochastic resampling, jackknife estimation, sensitivity analysis, and
the full three-run pipeline for Ozkan pTO.

- [`ozkan_pto_jackknife()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/ozkan_pto_jackknife.md)
  : Jackknife Analysis for Ozkan's pTO Index (Islem 1 / Run 1)
- [`ozkan_pto_resample()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/ozkan_pto_resample.md)
  : Stochastic Resampling of Ozkan's pTO Index (Islem 2 / Run 2)
- [`ozkan_pto_sensitivity()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/ozkan_pto_sensitivity.md)
  : Sensitivity Analysis of Ozkan's pTO Index (Islem 3 / Run 3)
- [`ozkan_pto_full()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/ozkan_pto_full.md)
  : Full Ozkan pTO Pipeline (Islem 1 + 2 + 3)

## Batch Analysis & Comparison

Multi-site batch processing and multi-community index comparison.
Computes species count and 14 diversity indices per site/community.

- [`batch_analysis()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/batch_analysis.md)
  : Batch Analysis from a Single Data Frame
- [`compare_indices()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/compare_indices.md)
  : Compare All Diversity Indices Side by Side

## Simulation & Rarefaction

Funnel plot simulation for expected index distributions and sample-based
rarefaction curves with bootstrap confidence intervals.

- [`simulate_td()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/simulate_td.md)
  : Simulate Expected AvTD/VarTD Under Random Sampling
- [`rarefaction_taxonomic()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/rarefaction_taxonomic.md)
  : Taxonomic Diversity Rarefaction

## Visualization

Seven specialized plot types for taxonomic diversity analysis.

- [`plot_funnel()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/plot_funnel.md)
  : Funnel Plot for AvTD/VarTD
- [`plot_rarefaction()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/plot_rarefaction.md)
  : Plot Taxonomic Rarefaction Curve
- [`plot_iteration()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/plot_iteration.md)
  : Plot pTO Iteration Results from Run 2 or Run 3
- [`plot_radar()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/plot_radar.md)
  : Radar (Spider) Chart for Multi-Community Index Comparison
- [`plot_heatmap()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/plot_heatmap.md)
  : Plot Taxonomic Distance Heatmap
- [`plot_bubble()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/plot_bubble.md)
  : Bubble Chart of Species Contributions to Diversity
- [`plot_taxonomic_tree()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/plot_taxonomic_tree.md)
  : Plot Taxonomic Tree as a Dendrogram

## Utility Functions

Helper functions for taxonomic tree construction and pairwise taxonomic
distance computation.

- [`build_tax_tree()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/build_tax_tree.md)
  : Build a Taxonomic Tree from Species Data
- [`tax_distance_matrix()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/tax_distance_matrix.md)
  : Compute Taxonomic Distance Matrix

## Datasets

Example ecological community datasets for demonstrations.

- [`anatolian_trees`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/anatolian_trees.md)
  : Anatolian Forest Trees: Multi-Site Species Data
- [`gazi_comm`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/gazi_comm.md)
  : Example Community Vector: 8 Anatolian Tree Species
- [`gazi_gytk`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/gazi_gytk.md)
  : Example Taxonomy: 8 Anatolian Tree Species
