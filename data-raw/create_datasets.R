## ---------------------------------------------------------------
## data-raw/create_datasets.R
## Build the example datasets and save them under data/
## ---------------------------------------------------------------

# --- 1. anatolian_trees: multi-site dataset ----------------------
#
# 3 sites, 20 unique species, Westhoff & van der Maarel (1973)
# cover-abundance scale (1-9).
# Real Anatolian forest species with real taxonomic classification.

species_pool <- data.frame(
  Species = c(
    "Pinus_nigra",
    "Pinus_brutia",
    "Cedrus_libani",
    "Abies_nordmanniana",
    "Juniperus_excelsa",
    "Quercus_cerris",
    "Quercus_robur",
    "Quercus_petraea",
    "Fagus_orientalis",
    "Carpinus_betulus",
    "Alnus_glutinosa",
    "Castanea_sativa",
    "Tilia_tomentosa",
    "Acer_platanoides",
    "Fraxinus_excelsior",
    "Platanus_orientalis",
    "Populus_tremula",
    "Sorbus_torminalis",
    "Ostrya_carpinifolia",
    "Taxus_baccata"
  ),
  Genus = c(
    "Pinus",
    "Pinus",
    "Cedrus",
    "Abies",
    "Juniperus",
    "Quercus",
    "Quercus",
    "Quercus",
    "Fagus",
    "Carpinus",
    "Alnus",
    "Castanea",
    "Tilia",
    "Acer",
    "Fraxinus",
    "Platanus",
    "Populus",
    "Sorbus",
    "Ostrya",
    "Taxus"
  ),
  Family = c(
    "Pinaceae",
    "Pinaceae",
    "Pinaceae",
    "Pinaceae",
    "Cupressaceae",
    "Fagaceae",
    "Fagaceae",
    "Fagaceae",
    "Fagaceae",
    "Betulaceae",
    "Betulaceae",
    "Fagaceae",
    "Malvaceae",
    "Sapindaceae",
    "Oleaceae",
    "Platanaceae",
    "Salicaceae",
    "Rosaceae",
    "Betulaceae",
    "Taxaceae"
  ),
  Order = c(
    "Pinales",
    "Pinales",
    "Pinales",
    "Pinales",
    "Pinales",
    "Fagales",
    "Fagales",
    "Fagales",
    "Fagales",
    "Fagales",
    "Fagales",
    "Fagales",
    "Malvales",
    "Sapindales",
    "Lamiales",
    "Proteales",
    "Malpighiales",
    "Rosales",
    "Fagales",
    "Pinales"
  ),
  Class = c(
    "Pinopsida",
    "Pinopsida",
    "Pinopsida",
    "Pinopsida",
    "Pinopsida",
    "Magnoliopsida",
    "Magnoliopsida",
    "Magnoliopsida",
    "Magnoliopsida",
    "Magnoliopsida",
    "Magnoliopsida",
    "Magnoliopsida",
    "Magnoliopsida",
    "Magnoliopsida",
    "Magnoliopsida",
    "Magnoliopsida",
    "Magnoliopsida",
    "Magnoliopsida",
    "Magnoliopsida",
    "Pinopsida"
  ),
  Phylum = c(
    "Pinophyta",
    "Pinophyta",
    "Pinophyta",
    "Pinophyta",
    "Pinophyta",
    "Magnoliophyta",
    "Magnoliophyta",
    "Magnoliophyta",
    "Magnoliophyta",
    "Magnoliophyta",
    "Magnoliophyta",
    "Magnoliophyta",
    "Magnoliophyta",
    "Magnoliophyta",
    "Magnoliophyta",
    "Magnoliophyta",
    "Magnoliophyta",
    "Magnoliophyta",
    "Magnoliophyta",
    "Pinophyta"
  ),
  Kingdom = c(
    rep("Plantae", 20)
  ),
  stringsAsFactors = FALSE
)


# Site 1: Mixed forest (conifer + broadleaf, 12 species)
site1_species <- c(
  "Pinus_nigra", "Cedrus_libani", "Abies_nordmanniana",
  "Quercus_cerris", "Quercus_petraea", "Fagus_orientalis",
  "Carpinus_betulus", "Castanea_sativa", "Tilia_tomentosa",
  "Acer_platanoides", "Sorbus_torminalis", "Taxus_baccata"
)
site1_abundance <- c(5, 4, 3, 7, 6, 8, 4, 3, 2, 2, 1, 1)

# Site 2: Broadleaf forest (13 species, diverse families)
site2_species <- c(
  "Quercus_robur", "Quercus_cerris", "Quercus_petraea",
  "Fagus_orientalis", "Carpinus_betulus", "Alnus_glutinosa",
  "Castanea_sativa", "Tilia_tomentosa", "Acer_platanoides",
  "Fraxinus_excelsior", "Platanus_orientalis", "Populus_tremula",
  "Ostrya_carpinifolia"
)
site2_abundance <- c(6, 5, 4, 9, 7, 3, 5, 3, 4, 2, 2, 1, 1)

# Site 3: Conifer forest (8 species, low taxonomic diversity)
site3_species <- c(
  "Pinus_nigra", "Pinus_brutia", "Cedrus_libani",
  "Abies_nordmanniana", "Juniperus_excelsa", "Taxus_baccata",
  "Quercus_cerris", "Sorbus_torminalis"
)
site3_abundance <- c(8, 7, 5, 4, 3, 1, 2, 1)


# Combine into a single data frame
build_site_df <- function(site_name, species_list, abundances, pool) {
  idx <- match(species_list, pool$Species)
  df <- pool[idx, , drop = FALSE]
  df$Site <- site_name
  df$Abundance <- abundances
  rownames(df) <- NULL
  df[, c("Site", "Species", "Genus", "Family", "Order",
         "Class", "Phylum", "Kingdom", "Abundance")]
}

anatolian_trees <- rbind(
  build_site_df("Mixed_Forest",     site1_species, site1_abundance, species_pool),
  build_site_df("Broadleaf_Forest", site2_species, site2_abundance, species_pool),
  build_site_df("Conifer_Forest",   site3_species, site3_abundance, species_pool)
)

rownames(anatolian_trees) <- NULL

# Save
save(anatolian_trees, file = "data/anatolian_trees.rda", compress = "xz")

cat("anatolian_trees created:\n")
cat("  Size:", nrow(anatolian_trees), "rows,", ncol(anatolian_trees), "columns\n")
cat("  Sites:", paste(unique(anatolian_trees$Site), collapse = ", "), "\n")
cat("  Total unique species:", length(unique(anatolian_trees$Species)), "\n\n")


# --- 2. campus_comm / campus_gytk: single-site quick demo ------------
#
# A single community vector + taxonomy table.
# 8 species, 3 taxonomic levels.
# Close in structure to the example in Ozkan (2018).

campus_comm <- c(
  Pinus_nigra       = 4,
  Pinus_brutia      = 2,
  Cedrus_libani     = 3,
  Quercus_cerris    = 1,
  Quercus_robur     = 2,
  Fagus_orientalis  = 3,
  Juniperus_excelsa = 2,
  Carpinus_betulus  = 2
)

campus_gytk <- data.frame(
  Species = c("Pinus_nigra", "Pinus_brutia", "Cedrus_libani",
              "Quercus_cerris", "Quercus_robur", "Fagus_orientalis",
              "Juniperus_excelsa", "Carpinus_betulus"),
  Genus   = c("Pinus", "Pinus", "Cedrus",
              "Quercus", "Quercus", "Fagus",
              "Juniperus", "Carpinus"),
  Family  = c("Pinaceae", "Pinaceae", "Pinaceae",
              "Fagaceae", "Fagaceae", "Fagaceae",
              "Cupressaceae", "Betulaceae"),
  Order   = c("Pinales", "Pinales", "Pinales",
              "Fagales", "Fagales", "Fagales",
              "Pinales", "Fagales"),
  stringsAsFactors = FALSE
)

save(campus_comm, file = "data/campus_comm.rda", compress = "xz")
save(campus_gytk, file = "data/campus_gytk.rda", compress = "xz")

cat("campus_comm created:\n")
cat("  Species count:", length(campus_comm), "\n")
cat("  Abundances:", paste(campus_comm, collapse = ", "), "\n\n")

cat("campus_gytk created:\n")
cat("  Size:", nrow(campus_gytk), "rows,", ncol(campus_gytk), "columns\n")
cat("  Levels:", paste(names(campus_gytk)[-1], collapse = ", "), "\n")

# --- 3. med_forest: multi-site Mediterranean vegetation ----------
#
# 92 sample plots, 72 woody/herbaceous species from a Mediterranean
# forest / maquis ecosystem. Cover-abundance scale (1-9), seven
# taxonomic ranks (species to kingdom).
#
# Derived from the plot-based survey in data-raw/Med_Dat_source.xlsx
# (wide species-by-plot matrix). The derivation applied to that source:
#   * trimmed stray non-breaking spaces from all taxonomy cells;
#   * built Species as the Genus_epithet binomial from the full name;
#   * corrected one source error: the genus of Thymbra spicata, which
#     was mislabelled "Thymus", was set to "Thymbra" (Lamiaceae);
#   * pivoted wide -> long, keeping only present species (Abundance > 0);
#   * relabelled the 92 plots plot_01 .. plot_92.
# The resulting long table is stored in data-raw/med_forest.csv, which
# this script reads so the .rda can be rebuilt reproducibly.

med_forest <- read.csv("data-raw/med_forest.csv", stringsAsFactors = FALSE)
med_forest$Abundance <- as.integer(med_forest$Abundance)
rownames(med_forest) <- NULL

save(med_forest, file = "data/med_forest.rda", compress = "xz")

cat("med_forest created:\n")
cat("  Size:", nrow(med_forest), "rows,", ncol(med_forest), "columns\n")
cat("  Sites:", length(unique(med_forest$Site)), "plots\n")
cat("  Unique species:", length(unique(med_forest$Species)), "\n")
cat("  Families:", length(unique(med_forest$Family)),
    "| Orders:", length(unique(med_forest$Order)), "\n\n")


cat("\nAll datasets saved under data/.\n")
