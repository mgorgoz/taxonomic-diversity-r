# =============================================================================
# Özkan pTO (İşlem 1) Fonksiyonu Testleri
# ozkan_pto() ve pto_components() fonksiyonlarının doğruluğunu kontrol eder.
#
# Bu testler Özkan (2018) makalesindeki deterministik (İşlem 1) hesaplamayı
# test eder. 4 bileşen üretilir: uTO, TO, uTO+, TO+
# =============================================================================

test_that("ozkan_pto returns all four components", {
  # 8 türlü, 3 cinsli, 1 familyalı örnek topluluk
  comm <- c(sp1 = 4, sp2 = 2, sp3 = 3, sp4 = 1, sp5 = 2,
            sp6 = 3, sp7 = 2, sp8 = 2)
  tax <- data.frame(
    Species = paste0("sp", 1:8),
    Genus   = c("G1", "G1", "G1", "G2", "G2", "G3", "G3", "G3"),
    Family  = c("F1", "F1", "F1", "F1", "F1", "F1", "F1", "F1"),
    stringsAsFactors = FALSE
  )

  result <- ozkan_pto(comm, tax)

  # Fonksiyon bir liste dönmeli
  expect_type(result, "list")

  # Listede 10 bileşen olmalı:
  # uTO, TO, uTO_plus, TO_plus = tüm seviyeler
  # uTO_max, TO_max, uTO_plus_max, TO_plus_max = bilgi veren seviyeler
  # Ed_levels, max_informative_level
  expect_true("uTO" %in% names(result))
  expect_true("TO" %in% names(result))
  expect_true("uTO_plus" %in% names(result))
  expect_true("TO_plus" %in% names(result))
  expect_true("uTO_max" %in% names(result))
  expect_true("TO_max" %in% names(result))
  expect_true("uTO_plus_max" %in% names(result))
  expect_true("TO_plus_max" %in% names(result))
  expect_true("Ed_levels" %in% names(result))
  expect_true("max_informative_level" %in% names(result))
})

test_that("ozkan_pto satisfies ordering constraints", {
  # Özkan (2018)'e göre şu sıralama her zaman geçerli olmalı:
  # TO+ >= TO    (mesafe >= çeşitlilik, ağırlıklı)
  # uTO+ >= uTO  (mesafe >= çeşitlilik, ağırlıksız)
  # TO >= uTO    (ağırlıklı >= ağırlıksız, çeşitlilik)
  # TO+ >= uTO+  (ağırlıklı >= ağırlıksız, mesafe)
  comm <- c(sp1 = 4, sp2 = 2, sp3 = 3, sp4 = 1, sp5 = 2,
            sp6 = 3, sp7 = 2, sp8 = 2)
  tax <- data.frame(
    Species = paste0("sp", 1:8),
    Genus   = c("G1", "G1", "G1", "G2", "G2", "G3", "G3", "G3"),
    Family  = c("F1", "F1", "F1", "F1", "F1", "F1", "F1", "F1"),
    stringsAsFactors = FALSE
  )

  r <- ozkan_pto(comm, tax)

  expect_true(r$TO_plus >= r$TO)
  expect_true(r$uTO_plus >= r$uTO)
  expect_true(r$TO >= r$uTO)
  expect_true(r$TO_plus >= r$uTO_plus)
})

test_that("ozkan_pto returns zero for single species", {
  # Tek türlü toplulukta çeşitlilik sıfır olmalı
  # Çünkü karşılaştırma yapılacak başka tür yok
  comm <- c(sp1 = 10)
  tax <- data.frame(
    Species = "sp1",
    Genus = "G1",
    stringsAsFactors = FALSE
  )

  r <- ozkan_pto(comm, tax)
  expect_equal(r$uTO, 0)
  expect_equal(r$TO, 0)
  expect_equal(r$uTO_max, 0)
  expect_equal(r$TO_max, 0)
  expect_equal(r$uTO_plus_max, 0)
  expect_equal(r$TO_plus_max, 0)
  expect_equal(r$max_informative_level, 0L)
})

test_that("ozkan_pto validates input", {
  # Hatalı girdilerde fonksiyon anlamlı hata mesajları vermeli
  tax <- data.frame(
    Species = c("sp1", "sp2"),
    Genus = c("G1", "G2"),
    stringsAsFactors = FALSE
  )

  # İsimsiz vektör: tür isimleri olmadan taksonomi eşleştirilemez
  expect_error(ozkan_pto(c(10, 5), tax), "named vector")

  # Negatif bolluk: doğada bolluk negatif olamaz
  expect_error(ozkan_pto(c(sp1 = -1, sp2 = 5), tax), "non-negative")

  # Taksonomide olmayan tür adı: eşleşme bulunamaz
  expect_error(ozkan_pto(c(sp_x = 10), tax), "not found")
})

test_that("pto_components returns named numeric vector with 8 elements", {
  # pto_components() fonksiyonu ozkan_pto()'nun kısa yol versiyonu
  # 8 elemanlı isimli sayısal vektör döner (Excel makrosundaki Run 1+2+3)
  comm <- c(sp1 = 4, sp2 = 2, sp3 = 3)
  tax <- data.frame(
    Species = paste0("sp", 1:3),
    Genus   = c("G1", "G1", "G2"),
    Family  = c("F1", "F1", "F1"),
    stringsAsFactors = FALSE
  )

  result <- pto_components(comm, tax)

  # Sayısal (double) tipinde olmalı
  expect_type(result, "double")

  # 8 eleman olmalı
  expect_length(result, 8)

  # İsimleri doğru sırada olmalı
  expect_equal(names(result), c("uTO", "TO", "uTO_plus", "TO_plus",
                                "uTO_max", "TO_max", "uTO_plus_max", "TO_plus_max"))
})

test_that("ozkan_pto uses presence-based (equal weight) entropy at species level", {
  # pTO formülünde tür seviyesinde her türe EŞİT AĞIRLIK verilir (1/S)
  # Bolluk bilgisi doğrudan kullanılmaz — sadece "dilimleme" aşamasında
  # hangi türlerin hayatta kaldığını belirler
  #
  # Bu yüzden tür seviyesi Deng entropisi = ln(S) olmalı
  # Bolluklar ne olursa olsun (100, 1, 1, 1, 1) -> Ed_S = ln(5)
  comm <- c(sp1 = 100, sp2 = 1, sp3 = 1, sp4 = 1, sp5 = 1)
  tax <- data.frame(
    Species = paste0("sp", 1:5),
    Genus   = c("G1", "G1", "G2", "G2", "G3"),
    Family  = rep("F1", 5),
    stringsAsFactors = FALSE
  )
  r <- ozkan_pto(comm, tax)

  # nk=0 diliminde 5 tür var → Ed_S = ln(5) = 1.6094
  expect_equal(r$Ed_levels[["Species"]], log(5), tolerance = 1e-10)
})

test_that("ozkan_pto genus-level Deng entropy with equal proportions", {
  # 6 tür, 3 cins (her cinste 2'şer tür): cins seviyesi Deng entropisi
  # m(Fi) = 2/6 = 1/3 (her cinsin oranı)
  # |Fi| = 2 (her cinste 2 tür)
  #
  # Ed_cins = -3 × (1/3) × ln((1/3) / (2²-1))
  #         = -3 × (1/3) × ln(1/9)
  #         = -ln(1/9) = ln(9) = 2.197
  comm <- setNames(rep(1, 6), paste0("sp", 1:6))
  tax <- data.frame(
    Species = paste0("sp", 1:6),
    Genus = rep(c("G1", "G2", "G3"), each = 2),
    Family = rep("F1", 6),
    stringsAsFactors = FALSE
  )
  r <- ozkan_pto(comm, tax)

  # Cins seviyesi Deng entropisi = ln(9) olmalı
  expect_equal(r$Ed_levels[["Genus"]], log(9), tolerance = 1e-10)
})

test_that("ozkan_pto matches Ozkan (2018) hypothetical examples", {
  # Özkan (2018) makalesinin Bölüm 4'teki hipotetik örnek topluluklar
  # Bu değerler makaledeki tablolarla doğrulanmıştır

  # Topluluk A: 12 tür, 3 cins (her cinste 4 tür), 1 familya
  # Beklenen uTO+ = ln(54.6) ≈ 4.00003
  comm_A <- setNames(rep(1, 12), paste0("sp", 1:12))
  tax_A <- data.frame(
    Species = paste0("sp", 1:12),
    Genus = rep(c("G1", "G2", "G3"), each = 4),
    Family = rep("F1", 12),
    stringsAsFactors = FALSE
  )
  r_A <- ozkan_pto(comm_A, tax_A)
  expect_equal(r_A$uTO_plus, log(54.6), tolerance = 1e-4)

  # Topluluk B: 6 tür, 3 cins (her cinste 2 tür), 1 familya
  # Beklenen uTO+ = ln(35) ≈ 3.5554
  comm_B <- setNames(rep(1, 6), paste0("sp", 1:6))
  tax_B <- data.frame(
    Species = paste0("sp", 1:6),
    Genus = rep(c("G1", "G2", "G3"), each = 2),
    Family = rep("F1", 6),
    stringsAsFactors = FALSE
  )
  r_B <- ozkan_pto(comm_B, tax_B)
  expect_equal(r_B$uTO_plus, log(35), tolerance = 1e-4)

  # A topluluğu B'den daha çeşitli olmalı (12 tür > 6 tür)
  expect_true(r_A$uTO_plus > r_B$uTO_plus)
})


# --- max_level parametresi testleri ---

test_that("ozkan_pto max_level='auto' detects informative levels", {
  # 5 tür, 2 cins, 1 familya — Family seviyesinde Ed=0 (tek grup)
  # Bu yüzden max_informative_level = 1 (sadece Genus bilgi veriyor)
  comm <- c(sp1 = 4, sp2 = 2, sp3 = 3, sp4 = 1, sp5 = 2)
  tax <- data.frame(
    Species = paste0("sp", 1:5),
    Genus   = c("G1", "G1", "G1", "G2", "G2"),
    Family  = c("F1", "F1", "F1", "F1", "F1"),
    stringsAsFactors = FALSE
  )

  r <- ozkan_pto(comm, tax, max_level = "auto")

  # Family seviyesinde tek grup var (F1) -> Ed = 0
  # max_informative_level = 1 olmali (Genus)
  expect_equal(r$max_informative_level, 1L)

  # max versiyonlari <= full versiyonlar olmali
  # (daha az seviye carpimda = daha kucuk veya esit deger)
  expect_true(r$uTO_plus_max <= r$uTO_plus + 1e-10)
  expect_true(r$TO_plus_max <= r$TO_plus + 1e-10)
})


test_that("ozkan_pto max_level integer restricts levels", {
  # 5 tür, 3 seviye (Genus, Family, Order)
  comm <- c(sp1 = 4, sp2 = 2, sp3 = 3, sp4 = 1, sp5 = 2)
  tax <- data.frame(
    Species = paste0("sp", 1:5),
    Genus   = c("G1", "G1", "G2", "G2", "G3"),
    Family  = c("F1", "F1", "F1", "F2", "F2"),
    Order   = c("O1", "O1", "O1", "O1", "O1"),
    stringsAsFactors = FALSE
  )

  r_full <- ozkan_pto(comm, tax)
  r_lim1 <- ozkan_pto(comm, tax, max_level = 1)
  r_lim2 <- ozkan_pto(comm, tax, max_level = 2)

  # max_level=1 -> sadece Genus seviyesi
  # max versiyonlari daha kucuk veya esit
  expect_true(r_lim1$uTO_plus_max <= r_full$uTO_plus + 1e-10)

  # max_level=2 >= max_level=1 (daha fazla seviye = daha buyuk)
  expect_true(r_lim2$uTO_plus_max >= r_lim1$uTO_plus_max - 1e-10)

  # Full versiyonlar max_level'den bagimsiz — ayni kalmali
  expect_equal(r_lim1$uTO, r_full$uTO, tolerance = 1e-10)
  expect_equal(r_lim2$uTO_plus, r_full$uTO_plus, tolerance = 1e-10)
})


test_that("ozkan_pto max_level validates input", {
  comm <- c(sp1 = 4, sp2 = 2, sp3 = 3)
  tax <- data.frame(
    Species = paste0("sp", 1:3),
    Genus   = c("G1", "G1", "G2"),
    Family  = c("F1", "F1", "F1"),
    stringsAsFactors = FALSE
  )

  # max_level = 0 -> hata (en az 1 olmali)
  expect_error(ozkan_pto(comm, tax, max_level = 0), "must be between")

  # max_level > mevcut seviye sayisi -> hata
  expect_error(ozkan_pto(comm, tax, max_level = 10), "must be between")

  # Gecersiz tip
  expect_error(ozkan_pto(comm, tax, max_level = TRUE), "must be NULL")
})


test_that("ozkan_pto max versions equal full when all levels informative", {
  # Tum seviyeler bilgi veriyorsa (Ed > 0), max = full olmali
  comm <- c(sp1 = 4, sp2 = 2, sp3 = 3, sp4 = 1, sp5 = 2, sp6 = 1)
  tax <- data.frame(
    Species = paste0("sp", 1:6),
    Genus   = c("G1", "G1", "G2", "G2", "G3", "G3"),
    Family  = c("F1", "F1", "F1", "F2", "F2", "F2"),
    stringsAsFactors = FALSE
  )

  r <- ozkan_pto(comm, tax, max_level = "auto")

  # 2 family (F1, F2) -> Family seviyesinde de Ed > 0
  # max = full olmali
  expect_equal(r$uTO_max, r$uTO, tolerance = 1e-10)
  expect_equal(r$TO_max, r$TO, tolerance = 1e-10)
  expect_equal(r$uTO_plus_max, r$uTO_plus, tolerance = 1e-10)
  expect_equal(r$TO_plus_max, r$TO_plus, tolerance = 1e-10)
})


test_that("pto_components returns max values matching ozkan_pto auto", {
  comm <- c(sp1 = 4, sp2 = 2, sp3 = 3, sp4 = 1, sp5 = 2)
  tax <- data.frame(
    Species = paste0("sp", 1:5),
    Genus   = c("G1", "G1", "G1", "G2", "G2"),
    Family  = c("F1", "F1", "F1", "F1", "F1"),
    stringsAsFactors = FALSE
  )

  pto <- pto_components(comm, tax)
  full <- ozkan_pto(comm, tax, max_level = "auto")

  expect_equal(unname(pto["uTO_max"]), full$uTO_max, tolerance = 1e-10)
  expect_equal(unname(pto["TO_max"]), full$TO_max, tolerance = 1e-10)
  expect_equal(unname(pto["uTO_plus_max"]), full$uTO_plus_max, tolerance = 1e-10)
  expect_equal(unname(pto["TO_plus_max"]), full$TO_plus_max, tolerance = 1e-10)
})
