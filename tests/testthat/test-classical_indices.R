# =============================================================================
# Klasik Çeşitlilik İndeksleri ve Clarke & Warwick İndeksleri Testleri
#
# Bu dosya şu fonksiyonları test eder:
#   shannon()     — Shannon çeşitlilik indeksi (H')
#   simpson()     — Simpson çeşitlilik indeksi (D, 1-D, 1/D)
#   delta()       — Taksonomik çeşitlilik (Δ) — bolluk ağırlıklı
#   delta_star()  — Taksonomik ayırt edicilik (Δ*) — aynı tür çiftleri hariç
#   avtd()        — Ortalama taksonomik ayırt edicilik (Δ+/AvTD) — bulunma/bulunmama
#   vartd()       — Taksonomik ayırt edicilik varyansı (Λ+/VarTD)
# =============================================================================


# =============================================================================
# Shannon Çeşitlilik İndeksi (H') Testleri
# Formül: H' = -Σ pᵢ × ln(pᵢ)
# =============================================================================

test_that("shannon calculates correctly", {
  # Eşit bolluklu 4 tür: H' = ln(4) = 1.386
  # Her türün oranı p = 10/40 = 0.25
  # H' = -4 × 0.25 × ln(0.25) = ln(4)
  comm <- c(10, 10, 10, 10)
  expect_equal(shannon(comm), log(4), tolerance = 1e-10)

  # Tek tür: çeşitlilik yok, H' = 0
  expect_equal(shannon(c(10)), 0)

  # 2 tabanlı logaritma ile (bit cinsinden): 2 eşit tür → H' = 1 bit
  comm <- c(10, 10)
  expect_equal(shannon(comm, base = 2), 1, tolerance = 1e-10)
})

test_that("shannon handles edge cases", {
  # Sıfır bolluklu türler otomatik çıkarılmalı
  # [10, 0, 5, 0, 8] ile [10, 5, 8] aynı sonucu vermeli
  comm <- c(10, 0, 5, 0, 8)
  expect_equal(shannon(comm), shannon(c(10, 5, 8)))

  # Hepsi sıfır: tür yok, çeşitlilik = 0
  expect_equal(shannon(c(0, 0, 0)), 0)
})

test_that("shannon validates input", {
  # Negatif bolluk: doğada bolluk negatif olamaz
  expect_error(shannon(c(-1, 5)), "non-negative")

  # Sayısal olmayan değer: bolluk sayı olmalı
  expect_error(shannon("abc"), "numeric")
})


# =============================================================================
# Simpson Çeşitlilik İndeksi Testleri
# Dominans: D = Σ pᵢ²
# Gini-Simpson: 1 - D
# Ters Simpson: 1/D
# =============================================================================

test_that("simpson calculates correctly", {
  # Eşit bolluklu 4 tür:
  # D = 4 × (0.25)² = 4 × 0.0625 = 0.25
  # Gini-Simpson = 1 - 0.25 = 0.75
  # Ters Simpson = 1/0.25 = 4
  comm <- c(10, 10, 10, 10)
  expect_equal(simpson(comm, "dominance"), 0.25, tolerance = 1e-10)
  expect_equal(simpson(comm, "gini_simpson"), 0.75, tolerance = 1e-10)
  expect_equal(simpson(comm, "inverse"), 4, tolerance = 1e-10)

  # Tek tür: D = 1 (tam dominans — sadece 1 tür var)
  expect_equal(simpson(c(10), "dominance"), 1)
})

test_that("simpson validates input", {
  # Negatif bolluk: hata vermeli
  expect_error(simpson(c(-1, 5)), "non-negative")
})


# =============================================================================
# Shannon Yanlılık Düzeltme (Bias Correction) Testleri
# =============================================================================

test_that("shannon miller_madow adds positive correction", {
  comm <- c(10, 5, 8, 3, 12)
  H_naive <- shannon(comm)
  H_mm <- shannon(comm, correction = "miller_madow")
  expect_true(H_mm > H_naive)
  N <- sum(comm)
  S_obs <- length(comm)
  expected <- H_naive + (S_obs - 1) / (2 * N)
  expect_equal(H_mm, expected, tolerance = 1e-10)
})

test_that("shannon miller_madow respects base parameter", {
  comm <- c(10, 5, 8, 3, 12)
  H_mm_nat <- shannon(comm, correction = "miller_madow")
  H_mm_2 <- shannon(comm, base = 2, correction = "miller_madow")
  expect_equal(H_mm_2, H_mm_nat / log(2), tolerance = 1e-10)
})

test_that("shannon grassberger uses digamma", {
  comm <- c(10, 5, 8, 3, 12)
  H_g <- shannon(comm, correction = "grassberger")
  N <- sum(comm)
  expected <- log(N) - (1 / N) * sum(comm * digamma(comm))
  expect_equal(H_g, expected, tolerance = 1e-10)
})

test_that("shannon grassberger respects base parameter", {
  comm <- c(10, 5, 8, 3, 12)
  H_g_nat <- shannon(comm, correction = "grassberger")
  H_g_2 <- shannon(comm, base = 2, correction = "grassberger")
  expect_equal(H_g_2, H_g_nat / log(2), tolerance = 1e-10)
})

test_that("shannon chao_shen applies coverage correction", {
  comm <- c(10, 5, 8, 3, 12)
  H_cs <- shannon(comm, correction = "chao_shen")
  expect_type(H_cs, "double")
  expect_true(is.finite(H_cs))
  expect_true(H_cs > 0)
})

test_that("shannon chao_shen handles many singletons", {
  comm <- c(1, 1, 1, 1, 5, 10)
  H_cs <- shannon(comm, correction = "chao_shen")
  expect_true(is.finite(H_cs))
  expect_true(H_cs > 0)
})

test_that("shannon chao_shen handles all singletons", {
  comm <- c(1, 1, 1, 1, 1)
  H_cs <- shannon(comm, correction = "chao_shen")
  expect_true(is.finite(H_cs))
  expect_true(H_cs > 0)
})

test_that("all corrections return 0 for single species", {
  for (corr in c("none", "miller_madow", "grassberger", "chao_shen")) {
    expect_equal(shannon(c(10), correction = corr), 0)
  }
})

test_that("all corrections converge for large equal samples", {
  comm <- rep(1000, 5)
  H_naive <- shannon(comm)
  H_mm <- shannon(comm, correction = "miller_madow")
  H_g <- shannon(comm, correction = "grassberger")
  H_cs <- shannon(comm, correction = "chao_shen")
  expect_equal(H_mm, H_naive, tolerance = 0.01)
  expect_equal(H_g, H_naive, tolerance = 0.01)
  expect_equal(H_cs, H_naive, tolerance = 0.01)
})

test_that("shannon warns for non-integer input with correction", {
  comm <- c(0.5, 0.3, 0.2)
  expect_warning(shannon(comm, correction = "miller_madow"),
                 "non-integer|Non-integer")
})

test_that("shannon default correction is none (backward compatible)", {
  comm <- c(10, 5, 8, 3, 12)
  expect_equal(shannon(comm), shannon(comm, correction = "none"))
})


# =============================================================================
# Clarke & Warwick İndeksleri
# =============================================================================

# --- Δ (Delta): Taksonomik Çeşitlilik ---
# Tüm birey çiftleri arasındaki ortalama taksonomik mesafe
# Aynı türden birey çiftleri de dahildir (mesafe = 0)

test_that("delta calculates taxonomic diversity correctly", {
  # 3 tür, 2 cins, 1 familya, hepsi eşit bollukta
  comm <- c(sp1 = 3, sp2 = 3, sp3 = 3)
  tax <- data.frame(
    Species = c("sp1", "sp2", "sp3"),
    Genus = c("G1", "G1", "G2"),       # sp1-sp2 aynı cins, sp3 farklı
    Family = c("F1", "F1", "F1"),       # hepsi aynı familya
    stringsAsFactors = FALSE
  )

  d <- delta(comm, tax)

  # Sayısal değer dönmeli ve pozitif olmalı
  expect_type(d, "double")
  expect_true(d > 0)

  # Sadece aynı cinsteki türlerle karşılaştır: mesafe daha düşük olmalı
  # Çünkü taksonomik yayılım daha az
  comm2 <- c(sp1 = 3, sp2 = 3)
  tax2 <- data.frame(
    Species = c("sp1", "sp2"),
    Genus = c("G1", "G1"),              # ikisi de aynı cins
    Family = c("F1", "F1"),
    stringsAsFactors = FALSE
  )
  d2 <- delta(comm2, tax2)

  # Aynı cinsteki 2 tür, farklı cinsleri de içeren 3 türden daha düşük Δ vermeli
  expect_true(d2 < d)
})

test_that("delta returns 0 for single species", {
  # Tek tür: karşılaştırma yapılacak başka tür yok → Δ = 0
  comm <- c(sp1 = 10)
  tax <- data.frame(
    Species = "sp1", Genus = "G1", Family = "F1",
    stringsAsFactors = FALSE
  )
  expect_equal(delta(comm, tax), 0)
})


# --- Δ* (Delta star): Taksonomik Ayırt Edicilik ---
# Δ ile aynı mantık ama aynı türden birey çiftleri HARİÇ tutulur
# Bu sayede Δ* her zaman Δ'dan büyük veya eşit olur

test_that("delta_star calculates taxonomic distinctness correctly", {
  # 3 tür, farklı bolluklar
  comm <- c(sp1 = 5, sp2 = 3, sp3 = 2)
  tax <- data.frame(
    Species = c("sp1", "sp2", "sp3"),
    Genus = c("G1", "G1", "G2"),
    Family = c("F1", "F1", "F1"),
    stringsAsFactors = FALSE
  )

  ds <- delta_star(comm, tax)
  d <- delta(comm, tax)

  # Sayısal ve pozitif olmalı
  expect_type(ds, "double")
  expect_true(ds > 0)

  # MATEMATİKSEL KURAL: Δ* her zaman Δ'dan büyük veya eşit
  # Çünkü Δ'nın paydasında aynı tür çiftleri de var (mesafe = 0 katkısı)
  # bu ortalamayı aşağı çeker
  expect_true(ds >= d)
})

test_that("delta_star returns 0 for single species", {
  # Tek tür: karşılaştırma yapılamaz → Δ* = 0
  comm <- c(sp1 = 10)
  tax <- data.frame(
    Species = "sp1", Genus = "G1",
    stringsAsFactors = FALSE
  )
  expect_equal(delta_star(comm, tax), 0)
})


# --- Δ+ (AvTD): Ortalama Taksonomik Ayırt Edicilik ---
# Sadece bulunma/bulunmama verisi kullanır (bolluk bilgisi gerekmez)
# Formül: Δ+ = Σ Σ_{i<j} ω_ij / [S(S-1)/2]

test_that("avtd calculates average taxonomic distinctness", {
  # 4 tür, 2 cins, 2 familya
  tax <- data.frame(
    Species = c("sp1", "sp2", "sp3", "sp4"),
    Genus = c("G1", "G1", "G2", "G2"),
    Family = c("F1", "F1", "F1", "F2"),   # sp4 farklı familyada
    stringsAsFactors = FALSE
  )

  spp <- c("sp1", "sp2", "sp3", "sp4")
  result <- avtd(spp, tax)

  # Sayısal ve pozitif olmalı
  expect_type(result, "double")
  expect_true(result > 0)
})

test_that("avtd with all same-genus species equals weights[1]", {
  # Tüm türler aynı cinste → tüm çiftlerin mesafesi = weights[1] = 1
  # Dolayısıyla ortalama da = 1
  tax <- data.frame(
    Species = c("sp1", "sp2", "sp3"),
    Genus = c("G1", "G1", "G1"),     # hepsi aynı cins
    Family = c("F1", "F1", "F1"),
    stringsAsFactors = FALSE
  )

  result <- avtd(c("sp1", "sp2", "sp3"), tax)

  # Tüm mesafeler 1 → ortalama = 1
  expect_equal(result, 1)
})

test_that("avtd requires at least 2 species", {
  # En az 2 tür gerekli (çift oluşturmak için)
  tax <- data.frame(
    Species = "sp1", Genus = "G1",
    stringsAsFactors = FALSE
  )
  expect_error(avtd("sp1", tax), "At least 2")
})


# --- Λ+ (VarTD): Taksonomik Ayırt Edicilik Varyansı ---
# Çiftler arası mesafelerin Δ+ etrafındaki varyansı
# Yüksek Λ+ = bazı çiftler çok yakın, bazıları çok uzak (dengesiz dağılım)
# Düşük Λ+ = tüm çiftlerin mesafesi birbirine yakın (homojen dağılım)

test_that("vartd calculates variation in taxonomic distinctness", {
  # 4 tür: bazıları yakın, bazıları uzak → varyans > 0
  tax <- data.frame(
    Species = c("sp1", "sp2", "sp3", "sp4"),
    Genus = c("G1", "G1", "G2", "G2"),
    Family = c("F1", "F1", "F1", "F2"),
    stringsAsFactors = FALSE
  )

  spp <- c("sp1", "sp2", "sp3", "sp4")
  result <- vartd(spp, tax)

  # Sayısal ve sıfırdan büyük veya eşit olmalı (varyans negatif olamaz)
  expect_type(result, "double")
  expect_true(result >= 0)
})

test_that("vartd is zero when all pairwise distances are equal", {
  # Tüm çiftlerin mesafesi eşit olduğunda varyans = 0
  # 3 tür, hepsi farklı cinste ama aynı familyada → tüm mesafeler = 2
  tax <- data.frame(
    Species = c("sp1", "sp2", "sp3"),
    Genus = c("G1", "G2", "G3"),       # hepsi farklı cins
    Family = c("F1", "F1", "F1"),       # hepsi aynı familya
    stringsAsFactors = FALSE
  )

  result <- vartd(c("sp1", "sp2", "sp3"), tax)

  # Tüm mesafeler eşit (= 2) → varyans = 0
  expect_equal(result, 0, tolerance = 1e-10)
})

test_that("vartd requires at least 2 species", {
  # En az 2 tür gerekli
  tax <- data.frame(
    Species = "sp1", Genus = "G1",
    stringsAsFactors = FALSE
  )
  expect_error(vartd("sp1", tax), "At least 2")
})


# =============================================================================
# Hata Kontrolü Testleri — Yanlış Girdi Durumları
# =============================================================================

test_that("delta errors on species not in tax_tree", {
  # Taksonomide bulunmayan tür adı verildiğinde hata vermeli
  comm <- c(sp1 = 5, sp_unknown = 3)
  tax <- data.frame(
    Species = c("sp1", "sp2"),
    Genus = c("G1", "G2"),
    Family = c("F1", "F1"),
    stringsAsFactors = FALSE
  )
  expect_error(delta(comm, tax), "Species not found in tax_tree: sp_unknown")
})

test_that("delta_star errors on species not in tax_tree", {
  # Birden fazla eksik tür varsa hepsini listele
  comm <- c(sp1 = 5, sp_missing = 3, sp_also = 2)
  tax <- data.frame(
    Species = c("sp1", "sp2"),
    Genus = c("G1", "G2"),
    Family = c("F1", "F1"),
    stringsAsFactors = FALSE
  )
  expect_error(delta_star(comm, tax), "Species not found in tax_tree: sp_missing, sp_also")
})

test_that("delta errors on wrong weights length", {
  # Ağırlık vektörünün uzunluğu taksonomik seviye sayısına eşit olmalı
  # Burada 2 seviye var (Genus, Family) ama 3 ağırlık verilmiş
  comm <- c(sp1 = 5, sp2 = 3)
  tax <- data.frame(
    Species = c("sp1", "sp2"),
    Genus = c("G1", "G2"),
    Family = c("F1", "F1"),
    stringsAsFactors = FALSE
  )
  expect_error(delta(comm, tax, weights = c(1, 2, 3)), "must have length 2")
})

test_that("delta_star errors on wrong weights length", {
  # Aynı ağırlık uzunluğu kontrolü delta_star için de geçerli
  # 2 seviye var ama 1 ağırlık verilmiş
  comm <- c(sp1 = 5, sp2 = 3)
  tax <- data.frame(
    Species = c("sp1", "sp2"),
    Genus = c("G1", "G2"),
    Family = c("F1", "F1"),
    stringsAsFactors = FALSE
  )
  expect_error(delta_star(comm, tax, weights = c(1)), "must have length 2")
})
