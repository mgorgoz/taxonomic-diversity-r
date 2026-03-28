# Package index

## Quick Start

Start here. These two functions handle multi-site batch processing and
multi-community comparison — the fastest way to get results from your
data.

- [`batch_analysis()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/batch_analysis.md)
  : Batch Analysis from a Single Data Frame
- [`compare_indices()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/compare_indices.md)
  : Compare All Diversity Indices Side by Side

## Classical Diversity Indices

Shannon-Wiener (H’) and Gini-Simpson indices — the most widely used
diversity measures in ecology. Shannon includes three bias correction
methods: Miller-Madow, Grassberger, and Chao-Shen.

- [`shannon()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/shannon.md)
  : Shannon Diversity Index
- [`simpson()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/simpson.md)
  : Simpson Diversity Index

## Clarke & Warwick Taxonomic Distinctness

Path-length-based taxonomic distinctness measures. Delta and Delta\* are
abundance-weighted; AvTD and VarTD are presence/absence-based and
sample-size independent.

- [`delta()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/delta.md)
  : Taxonomic Diversity Index (Delta)
- [`delta_star()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/delta_star.md)
  : Taxonomic Distinctness (Delta\*)
- [`avtd()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/avtd.md)
  : Average Taxonomic Distinctness (Delta+)
- [`vartd()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/vartd.md)
  : Variation in Taxonomic Distinctness (Lambda+)

## Ozkan pTO — Core Functions

Deng entropy-based taxonomic diversity following Ozkan (2018). Produces
8 complementary indices (uTO, TO, uTO+, TO+ and their
max-informative-level variants).

- [`ozkan_pto()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/ozkan_pto.md)
  : Calculate Ozkan's Taxonomic Diversity Index (pTO)
- [`pto_components()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/pto_components.md)
  : Calculate All Eight pTO Components (Convenience Wrapper)
- [`deng_entropy_level()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/deng_entropy_level.md)
  : Calculate Deng Entropy at a Single Taxonomic Level

## Ozkan pTO — Full Pipeline

The complete three-run pipeline: deterministic calculation (Run 1),
stochastic resampling with slicing (Run 2), max-informative levels (Run
3), plus jackknife and sensitivity analysis.

- [`ozkan_pto_full()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/ozkan_pto_full.md)
  : Full Ozkan pTO Pipeline (Islem 1 + 2 + 3)
- [`ozkan_pto_resample()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/ozkan_pto_resample.md)
  : Stochastic Resampling of Ozkan's pTO Index (Islem 2 / Run 2)
- [`ozkan_pto_sensitivity()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/ozkan_pto_sensitivity.md)
  : Sensitivity Analysis of Ozkan's pTO Index (Islem 3 / Run 3)
- [`ozkan_pto_jackknife()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/ozkan_pto_jackknife.md)
  : Jackknife Analysis for Ozkan's pTO Index (Islem 1 / Run 1)

## Simulation & Significance Testing

Random subsampling from a master species pool to generate expected
distributions of AvTD and VarTD. Use with plot_funnel() for 95%
confidence funnels and statistical significance testing.

- [`simulate_td()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/simulate_td.md)
  : Simulate Expected AvTD/VarTD Under Random Sampling
- [`rarefaction_taxonomic()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/rarefaction_taxonomic.md)
  : Taxonomic Diversity Rarefaction

## Visualization

Seven ggplot2-based plot types covering significance testing,
rarefaction, resampling trajectories, community comparison, similarity
patterns, composition, and taxonomic structure.

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

## Taxonomic Tree & Distance

Build taxonomic classification trees from species-level data and compute
pairwise taxonomic distance matrices used by Clarke & Warwick and Ozkan
pTO methods.

- [`build_tax_tree()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/build_tax_tree.md)
  : Build a Taxonomic Tree from Species Data
- [`tax_distance_matrix()`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/tax_distance_matrix.md)
  : Compute Taxonomic Distance Matrix

## Example Datasets

Ready-to-use ecological community datasets from Mediterranean forest
ecosystems in Anatolia, Turkey. Includes abundance data and full
taxonomic classifications.

- [`anatolian_trees`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/anatolian_trees.md)
  : Anatolian Forest Trees: Multi-Site Species Data
- [`gazi_comm`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/gazi_comm.md)
  : Example Community Vector: 8 Anatolian Tree Species
- [`gazi_gytk`](https://mgorgoz.github.io/taxonomic-diversity-r/reference/gazi_gytk.md)
  : Example Taxonomy: 8 Anatolian Tree Species
