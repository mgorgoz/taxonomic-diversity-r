## ---------------------------------------------------------------
## data-raw/create_datasets.R
## Ornek veri setlerini olustur ve data/ klasorune kaydet
## ---------------------------------------------------------------

# --- 1. anatolian_trees: Coklu alan veri seti ----------------------
#
# 3 ornek alan, 20 benzersiz tur, Westhoff & van der Maarel (1973)
# bolluk degerleri (1-9).
# Gercek Anadolu orman turleri ve gercek taksonomik siniflandirma.

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


# Alan 1: Karisik orman (konifer + yaprakli, 12 tur)
site1_species <- c(
  "Pinus_nigra", "Cedrus_libani", "Abies_nordmanniana",
  "Quercus_cerris", "Quercus_petraea", "Fagus_orientalis",
  "Carpinus_betulus", "Castanea_sativa", "Tilia_tomentosa",
  "Acer_platanoides", "Sorbus_torminalis", "Taxus_baccata"
)
site1_abundance <- c(5, 4, 3, 7, 6, 8, 4, 3, 2, 2, 1, 1)

# Alan 2: Yaprakli orman (13 tur, cesitli familyalar)
site2_species <- c(
  "Quercus_robur", "Quercus_cerris", "Quercus_petraea",
  "Fagus_orientalis", "Carpinus_betulus", "Alnus_glutinosa",
  "Castanea_sativa", "Tilia_tomentosa", "Acer_platanoides",
  "Fraxinus_excelsior", "Platanus_orientalis", "Populus_tremula",
  "Ostrya_carpinifolia"
)
site2_abundance <- c(6, 5, 4, 9, 7, 3, 5, 3, 4, 2, 2, 1, 1)

# Alan 3: Konifer orman (8 tur, az taksonomik cesitlilik)
site3_species <- c(
  "Pinus_nigra", "Pinus_brutia", "Cedrus_libani",
  "Abies_nordmanniana", "Juniperus_excelsa", "Taxus_baccata",
  "Quercus_cerris", "Sorbus_torminalis"
)
site3_abundance <- c(8, 7, 5, 4, 3, 1, 2, 1)


# Veri cercevesini birlestir
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
  build_site_df("Karisik_Orman", site1_species, site1_abundance, species_pool),
  build_site_df("Yaprakli_Orman", site2_species, site2_abundance, species_pool),
  build_site_df("Konifer_Orman", site3_species, site3_abundance, species_pool)
)

rownames(anatolian_trees) <- NULL

# Kaydet
save(anatolian_trees, file = "data/anatolian_trees.rda", compress = "xz")

cat("anatolian_trees olusturuldu:\n")
cat("  Boyut:", nrow(anatolian_trees), "satir,", ncol(anatolian_trees), "sutun\n")
cat("  Alanlar:", paste(unique(anatolian_trees$Site), collapse = ", "), "\n")
cat("  Toplam benzersiz tur:", length(unique(anatolian_trees$Species)), "\n\n")


# --- 2. gazi_inger: Tek alan ornegi (batch_analysis icin degil) -------
#
# Tek bir topluluk vektoru + taksonomi tablosu
# Tez verisine benzer yapida, 8 tur, 3 taksonomik seviye
# Ozkan (2018) makalesindeki ornege yakin

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

cat("gazi_comm olusturuldu:\n")
cat("  Tur sayisi:", length(gazi_comm), "\n")
cat("  Bolluklar:", paste(gazi_comm, collapse = ", "), "\n\n")

cat("gazi_gytk olusturuldu:\n")
cat("  Boyut:", nrow(gazi_gytk), "satir,", ncol(gazi_gytk), "sutun\n")
cat("  Seviyeler:", paste(names(gazi_gytk)[-1], collapse = ", "), "\n")

cat("\nTum veri setleri data/ klasorune kaydedildi.\n")
