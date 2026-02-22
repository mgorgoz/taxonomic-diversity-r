# =============================================================================
# Deng Entropisi Fonksiyonu Testleri
# deng_entropy_level() fonksiyonunun doğru çalıştığını kontrol eder.
#
# Deng entropisi, Shannon entropisinin genelleştirilmiş halidir.
# Tür seviyesinde (|Fi|=1) Shannon'a eşittir.
# Üst taksonomik seviyelerde (|Fi|>1) Shannon'dan büyüktür.
# =============================================================================

test_that("deng_entropy_level returns Shannon entropy at species level", {
  # Tür seviyesinde Deng = Shannon olmalı (çünkü |Fi|=1)

  # 4 eşit bolluklu tür: H = ln(4) = 1.386
  # Her türün oranı p = 10/40 = 0.25
  # H = -4 × 0.25 × ln(0.25) = ln(4)
  result <- deng_entropy_level(c(10, 10, 10, 10))
  expect_equal(result, log(4), tolerance = 1e-10)

  # Tek tür: çeşitlilik yok, H = 0
  expect_equal(deng_entropy_level(c(10)), 0)

  # 2 eşit tür: H = ln(2) = 0.693
  expect_equal(deng_entropy_level(c(5, 5)), log(2), tolerance = 1e-10)
})

test_that("deng_entropy_level with group_sizes differs from Shannon", {
  # Üst taksonomik seviyelerde (cins, familya vb.) Deng ≠ Shannon
  # Çünkü her grubun içinde birden fazla tür var (|Fi| > 1)
  # Bu durumda formüldeki 2^|Fi|-1 payda devreye girer

  # 3 grup: bollukları [9, 3, 7], grup büyüklükleri [3, 2, 3]
  # Yani 1. grupta 3 tür, 2. grupta 2 tür, 3. grupta 3 tür var
  abund <- c(9, 3, 7)
  gs <- c(3, 2, 3)

  ed_deng <- deng_entropy_level(abund, group_sizes = gs)
  ed_shannon <- deng_entropy_level(abund)  # group_sizes yok = Shannon

  # Deng, Shannon'dan büyük olmalı (|Fi| > 1 olduğunda)
  # Çünkü 2^|Fi|-1 terimi entropi değerini artırır
  expect_true(ed_deng > ed_shannon)
})

test_that("deng_entropy_level handles zeros", {
  # Sıfır bolluklu türler hesaplamadan otomatik çıkarılmalı
  # [10, 0, 5, 0, 8] ile [10, 5, 8] aynı sonucu vermeli
  result <- deng_entropy_level(c(10, 0, 5, 0, 8))
  expected <- deng_entropy_level(c(10, 5, 8))
  expect_equal(result, expected)
})

test_that("deng_entropy_level validates input", {
  # Negatif değer verilince hata vermeli
  expect_error(deng_entropy_level(c(-1, 5)), "non-negative")

  # Sayısal olmayan değer verilince hata vermeli
  expect_error(deng_entropy_level("abc"), "numeric")
})
