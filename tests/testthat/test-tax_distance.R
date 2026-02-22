# =============================================================================
# Taksonomik Mesafe Matrisi ve Taksonomi Ağacı Oluşturucu Testleri
#
# Bu dosya şu fonksiyonları test eder:
#   tax_distance_matrix() — Türler arası çiftli taksonomik mesafe matrisi
#   build_tax_tree()      — Sınıflandırma tablosu (data.frame) oluşturucu
#
# Mesafe hesaplama mantığı:
#   İki tür arasındaki mesafe = ilk eşleşen taksonomik seviyenin ağırlığı
#   Aynı cins → ω=1, Aynı familya farklı cins → ω=2, vs.
# =============================================================================

test_that("tax_distance_matrix computes correct distances", {
  # 3 tür: sp_a ve sp_b aynı cinste, sp_c farklı cinste ama aynı familyada
  tax <- data.frame(
    Species = c("sp_a", "sp_b", "sp_c"),
    Genus = c("Gen1", "Gen1", "Gen2"),
    Family = c("Fam1", "Fam1", "Fam1"),
    stringsAsFactors = FALSE
  )

  # Varsayılan ağırlıklar: Genus=1, Family=2
  d <- tax_distance_matrix(tax)

  # sp_a ve sp_b aynı cins (Gen1) → ilk eşleşme Genus seviyesinde → ω = 1
  expect_equal(d["sp_a", "sp_b"], 1)

  # sp_a ve sp_c farklı cins (Gen1 ≠ Gen2), aynı familya → eşleşme Family'de → ω = 2
  expect_equal(d["sp_a", "sp_c"], 2)

  # Matris simetrik olmalı: ω(a,c) = ω(c,a)
  expect_equal(d["sp_a", "sp_c"], d["sp_c", "sp_a"])

  # Köşegen sıfır olmalı: bir türün kendisiyle mesafesi = 0
  expect_equal(d["sp_a", "sp_a"], 0)
})

test_that("tax_distance_matrix works with custom weights", {
  # 2 tür, hiçbir seviyede eşleşmeyen (farklı cins, farklı familya)
  tax <- data.frame(
    Species = c("sp_a", "sp_b"),
    Genus = c("Gen1", "Gen2"),
    Family = c("Fam1", "Fam2"),
    stringsAsFactors = FALSE
  )

  # Varsayılan ağırlıklar [1, 2]: hiçbir seviyede eşleşme yok → maksimum ağırlık = 2
  d1 <- tax_distance_matrix(tax)
  expect_equal(d1["sp_a", "sp_b"], 2)

  # Özel ağırlıklar [1, 3]: hiçbir eşleşme yok → maksimum ağırlık = 3
  d2 <- tax_distance_matrix(tax, weights = c(1, 3))
  expect_equal(d2["sp_a", "sp_b"], 3)
})


# =============================================================================
# build_tax_tree() Testleri — Sınıflandırma Tablosu Oluşturucu
# =============================================================================

test_that("build_tax_tree creates correct structure", {
  # 2 tür, 2 taksonomik seviye (Genus, Family)
  tree <- build_tax_tree(
    species = c("sp_a", "sp_b"),
    Genus = c("Gen1", "Gen2"),
    Family = c("Fam1", "Fam1")
  )

  # 3 sütun olmalı: Species + Genus + Family
  expect_equal(ncol(tree), 3)
  expect_equal(names(tree), c("Species", "Genus", "Family"))

  # 2 satır olmalı (her tür için 1)
  expect_equal(nrow(tree), 2)
})

test_that("build_tax_tree validates input", {
  # Tür sayısı ile taksonomik seviye uzunluğu eşleşmeli
  # 2 tür var ama sadece 1 cins adı verilmiş → hata
  expect_error(
    build_tax_tree(species = c("a", "b"), Genus = c("G1")),
    "does not match"
  )

  # En az 1 taksonomik seviye verilmeli (sadece tür adı yetmez)
  expect_error(build_tax_tree(species = c("a")), "At least one")
})

test_that("build_tax_tree rejects duplicate species", {
  # Aynı tür adı iki kez verilmiş → hata
  # Çünkü tür adları benzersiz olmalı (mesafe matrisinde satır/sütun isimleri)
  expect_error(
    build_tax_tree(species = c("sp1", "sp2", "sp1"), Genus = c("G1", "G2", "G1")),
    "Duplicate species names: sp1"
  )
})


# =============================================================================
# Hata Kontrolü — Yanlış Ağırlık Uzunluğu
# =============================================================================

test_that("tax_distance_matrix errors on wrong weights length", {
  # 2 taksonomik seviye var (Genus, Family) → 2 ağırlık gerekli
  tax <- data.frame(
    Species = c("sp_a", "sp_b"),
    Genus = c("Gen1", "Gen2"),
    Family = c("Fam1", "Fam1"),
    stringsAsFactors = FALSE
  )

  # 1 ağırlık verilmiş ama 2 gerekli → hata
  expect_error(tax_distance_matrix(tax, weights = c(1)), "must have length 2")

  # 4 ağırlık verilmiş ama 2 gerekli → hata
  expect_error(tax_distance_matrix(tax, weights = c(1, 2, 3, 4)), "must have length 2")
})
