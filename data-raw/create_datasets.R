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


# --- 2. gazi_comm / gazi_gytk: single-site quick demo ------------
#
# A single community vector + taxonomy table.
# 8 species, 3 taxonomic levels.
# Close in structure to the example in Ozkan (2018).

gazi_comm <- c(
  Pinus_nigra       = 4,
  Pinus_brutia      = 2,
  Cedrus_libani     = 3,
  Quercus_cerris    = 1,
  Quercus_robur     = 2,
  Fagus_orientalis  = 3,
  Juniperus_excelsa = 2,
  Carpinus_betulus  = 2
)

gazi_gytk <- data.frame(
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

save(gazi_comm, file = "data/gazi_comm.rda", compress = "xz")
save(gazi_gytk, file = "data/gazi_gytk.rda", compress = "xz")

cat("gazi_comm created:\n")
cat("  Species count:", length(gazi_comm), "\n")
cat("  Abundances:", paste(gazi_comm, collapse = ", "), "\n\n")

cat("gazi_gytk created:\n")
cat("  Size:", nrow(gazi_gytk), "rows,", ncol(gazi_gytk), "columns\n")
cat("  Levels:", paste(names(gazi_gytk)[-1], collapse = ", "), "\n")

cat("\nAll datasets saved under data/.\n")
