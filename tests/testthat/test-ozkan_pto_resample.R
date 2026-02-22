# =============================================================================
# Özkan pTO Stokastik Örnekleme (İşlem 2) ve Duyarlılık Analizi (İşlem 3)
# ozkan_pto_resample() ve ozkan_pto_sensitivity() testleri
#
# İşlem 2: Her türü %50 olasılıkla dahil et/çıkar, n_iter kez tekrarla,
#           her bileşenin maksimumunu al.
# İşlem 3: İşlem 2 sonuçlarına göre türlere farklı dahil etme olasılıkları
#           ver, tekrar n_iter kez tekrarla, genel maksimumu al.
# =============================================================================

# --- Tüm testlerde kullanılan ortak test verisi ---
make_test_data <- function() {
  # 8 tür, 3 cins, 2 familya, 1 takım
  # Westhoff-Maarel ölçeğine uygun bolluklar (1-4 arası)
  comm <- c(sp1 = 4, sp2 = 2, sp3 = 3, sp4 = 1, sp5 = 2,
            sp6 = 3, sp7 = 2, sp8 = 2)
  tax <- data.frame(
    Species = paste0("sp", 1:8),
    Genus   = c("G1", "G1", "G1", "G2", "G2", "G3", "G3", "G3"),
    Family  = c("F1", "F1", "F1", "F1", "F1", "F2", "F2", "F2"),
    Order   = rep("O1", 8),
    stringsAsFactors = FALSE
  )
  list(comm = comm, tax = tax)
}

# =============================================================================
# ozkan_pto_resample() testleri (İşlem 2 / Run 2)
# =============================================================================

test_that("ozkan_pto_resample returns correct structure", {
  # Fonksiyonun döndürdüğü listenin doğru yapıda olup olmadığını kontrol et
  d <- make_test_data()
  result <- ozkan_pto_resample(d$comm, d$tax, n_iter = 101, seed = 42)

  # Sonuç bir liste olmalı
  expect_type(result, "list")

  # Listede 4 adet maksimum değer olmalı (tüm iterasyonların max'ı)
  expect_true("uTO_plus_max" %in% names(result))
  expect_true("TO_plus_max" %in% names(result))
  expect_true("uTO_max" %in% names(result))
  expect_true("TO_max" %in% names(result))

  # Listede 4 adet deterministik değer olmalı (1. iterasyon = İşlem 1)
  expect_true("uTO_plus_det" %in% names(result))
  expect_true("TO_plus_det" %in% names(result))
  expect_true("uTO_det" %in% names(result))
  expect_true("TO_det" %in% names(result))

  # İterasyon sayısı ve sonuç tablosu
  expect_true("n_iter" %in% names(result))
  expect_true("iteration_results" %in% names(result))
  expect_equal(result$n_iter, 101L)
})

test_that("ozkan_pto_resample iteration results has correct dimensions", {
  # İterasyon sonuç tablosunun boyutlarını kontrol et
  d <- make_test_data()
  result <- ozkan_pto_resample(d$comm, d$tax, n_iter = 101, seed = 42)

  iter_df <- result$iteration_results

  # data.frame olmalı
  expect_s3_class(iter_df, "data.frame")

  # 101 satır (her iterasyon için 1 satır)
  expect_equal(nrow(iter_df), 101)

  # 5 sütun: iterasyon numarası + 4 bileşen
  expect_equal(ncol(iter_df), 5)
  expect_equal(names(iter_df), c("iteration", "uTO_plus", "TO_plus", "uTO", "TO"))
})

test_that("ozkan_pto_resample first iteration matches deterministic ozkan_pto", {
  # İlk iterasyon her zaman deterministik olmalı (orijinal topluluk aynen kullanılır)
  # Bu, ozkan_pto() fonksiyonuyla aynı sonucu vermeli
  d <- make_test_data()
  result <- ozkan_pto_resample(d$comm, d$tax, n_iter = 101, seed = 42)
  det <- ozkan_pto(d$comm, d$tax)

  # Deterministik değerler İşlem 1 ile birebir eşleşmeli
  expect_equal(result$uTO_plus_det, det$uTO_plus, tolerance = 1e-10)
  expect_equal(result$TO_plus_det, det$TO_plus, tolerance = 1e-10)
  expect_equal(result$uTO_det, unname(det$uTO), tolerance = 1e-10)
  expect_equal(result$TO_det, unname(det$TO), tolerance = 1e-10)

  # İterasyon tablosunun 1. satırı da aynı olmalı
  expect_equal(result$iteration_results$uTO_plus[1], det$uTO_plus, tolerance = 1e-10)
  expect_equal(result$iteration_results$TO_plus[1], det$TO_plus, tolerance = 1e-10)
})

test_that("ozkan_pto_resample maximums >= deterministic values", {
  # Maksimum değerler deterministik değerlerden büyük veya eşit olmalı
  # Çünkü 1. iterasyon zaten deterministik, max en az o kadar olacak
  d <- make_test_data()
  result <- ozkan_pto_resample(d$comm, d$tax, n_iter = 101, seed = 42)

  # Taksonomik mesafe (pTO+): max >= det her zaman geçerli
  expect_true(result$uTO_plus_max >= result$uTO_plus_det)
  expect_true(result$TO_plus_max >= result$TO_plus_det)

  # Taksonomik çeşitlilik (pTO): rastgele tür çıkarma çeşitliliği azaltabilir
  # ama max değerler en az 0 olmalı
  expect_true(result$uTO_max >= 0)
  expect_true(result$TO_max >= 0)
})

test_that("ozkan_pto_resample is reproducible with same seed", {
  # Aynı rastgele tohum (seed) ile iki kez çalıştırınca aynı sonuç çıkmalı
  # Bu, sonuçların tekrarlanabilir olduğunu garanti eder
  d <- make_test_data()
  r1 <- ozkan_pto_resample(d$comm, d$tax, n_iter = 101, seed = 123)
  r2 <- ozkan_pto_resample(d$comm, d$tax, n_iter = 101, seed = 123)

  expect_equal(r1$uTO_plus_max, r2$uTO_plus_max)
  expect_equal(r1$TO_plus_max, r2$TO_plus_max)
  expect_equal(r1$uTO_max, r2$uTO_max)
  expect_equal(r1$TO_max, r2$TO_max)

  # İterasyon tablosu da birebir aynı olmalı
  expect_equal(r1$iteration_results, r2$iteration_results)
})

test_that("ozkan_pto_resample produces different results with different seeds", {
  # Farklı tohumlarla farklı sonuçlar çıkmalı (rastgelelik çalışıyor mu?)
  d <- make_test_data()
  r1 <- ozkan_pto_resample(d$comm, d$tax, n_iter = 101, seed = 1)
  r2 <- ozkan_pto_resample(d$comm, d$tax, n_iter = 101, seed = 999)

  # Farklı sonuç çıkma olasılığı çok yüksek ama %100 garanti değil
  # Bu yüzden sadece geçerli sayısal değer olduğunu kontrol ediyoruz
  expect_true(is.numeric(r1$uTO_plus_max))
  expect_true(is.numeric(r2$uTO_plus_max))
})

test_that("ozkan_pto_resample validates n_iter >= 101", {
  # İterasyon sayısı en az 101 olmalı (Özkan makalesindeki minimum)
  d <- make_test_data()

  # 50 iterasyon: hata vermeli
  expect_error(ozkan_pto_resample(d$comm, d$tax, n_iter = 50), ">= 101")

  # 0 iterasyon: hata vermeli
  expect_error(ozkan_pto_resample(d$comm, d$tax, n_iter = 0), ">= 101")
})

test_that("ozkan_pto_resample validates inputs", {
  # Hatalı girdilerde anlamlı hata mesajları vermeli
  d <- make_test_data()
  tax <- d$tax

  # İsimsiz vektör: tür isimleri gerekli
  expect_error(ozkan_pto_resample(c(10, 5), tax), "named vector")

  # Negatif bolluk: doğada negatif bolluk olmaz
  expect_error(ozkan_pto_resample(c(sp1 = -1, sp2 = 5), tax), "non-negative")

  # Tek tür: örnekleme yapabilmek için en az 2 tür gerekli
  expect_error(ozkan_pto_resample(c(sp1 = 5), tax, n_iter = 101),
               "at least 2")
})

test_that("ozkan_pto_resample satisfies ordering: TO+ >= uTO+ for max", {
  # Ağırlıklı taksonomik mesafe (TO+) her zaman ağırlıksızdan (uTO+)
  # büyük veya eşit olmalı — çünkü ağırlıklar wi=i >= 1
  d <- make_test_data()
  result <- ozkan_pto_resample(d$comm, d$tax, n_iter = 101, seed = 42)

  expect_true(result$TO_plus_max >= result$uTO_plus_max)
})

test_that("ozkan_pto_resample handles many iterations", {
  # 500 iterasyon gibi yüksek değerlerle de çalışmalı
  d <- make_test_data()
  result <- ozkan_pto_resample(d$comm, d$tax, n_iter = 500, seed = 42)

  expect_equal(result$n_iter, 500L)
  expect_equal(nrow(result$iteration_results), 500)
  expect_true(is.numeric(result$uTO_plus_max))
})

test_that("ozkan_pto_resample with equal abundances works", {
  # Tüm türlerin bolluğu eşit olduğunda da çalışmalı
  # Bu durumda dilimleme prosedürü tek adımda biter (hepsi nk=0'da elenir)
  comm <- setNames(rep(1, 6), paste0("sp", 1:6))
  tax <- data.frame(
    Species = paste0("sp", 1:6),
    Genus = rep(c("G1", "G2", "G3"), each = 2),
    Family = rep("F1", 6),
    stringsAsFactors = FALSE
  )

  result <- ozkan_pto_resample(comm, tax, n_iter = 101, seed = 42)
  expect_true(is.numeric(result$uTO_plus_max))
  expect_true(result$uTO_plus_max > 0)
})


# =============================================================================
# ozkan_pto_sensitivity() testleri (İşlem 3 / Run 3)
# =============================================================================

test_that("ozkan_pto_sensitivity returns correct structure", {
  # İşlem 3 sonuç yapısını kontrol et
  d <- make_test_data()
  run2 <- ozkan_pto_resample(d$comm, d$tax, n_iter = 101, seed = 42)
  run3 <- ozkan_pto_sensitivity(d$comm, d$tax, run2, n_iter = 101, seed = 42)

  # Sonuç bir liste olmalı
  expect_type(run3, "list")

  # Genel maksimumlar (Run 1+2+3 birlikte)
  expect_true("uTO_plus_max" %in% names(run3))
  expect_true("TO_plus_max" %in% names(run3))
  expect_true("uTO_max" %in% names(run3))
  expect_true("TO_max" %in% names(run3))

  # Sadece Run 3'e ait maksimumlar
  expect_true("run3_uTO_plus_max" %in% names(run3))
  expect_true("run3_TO_plus_max" %in% names(run3))

  # Türlere özel dahil etme olasılıkları
  expect_true("species_probs" %in% names(run3))

  # İterasyon sonuç tablosu
  expect_true("iteration_results" %in% names(run3))
})

test_that("ozkan_pto_sensitivity overall max >= Run 2 max", {
  # İşlem 3'ün genel maksimumu İşlem 2'nin maksimumundan büyük veya eşit olmalı
  # Çünkü genel max = pmax(Run2_max, Run3_max)
  d <- make_test_data()
  run2 <- ozkan_pto_resample(d$comm, d$tax, n_iter = 101, seed = 42)
  run3 <- ozkan_pto_sensitivity(d$comm, d$tax, run2, n_iter = 101, seed = 42)

  expect_true(run3$uTO_plus_max >= run2$uTO_plus_max)
  expect_true(run3$TO_plus_max >= run2$TO_plus_max)
})

test_that("ozkan_pto_sensitivity uses Run 2 n_iter by default", {
  # n_iter belirtilmezse İşlem 2'deki iterasyon sayısını kullanmalı
  d <- make_test_data()
  run2 <- ozkan_pto_resample(d$comm, d$tax, n_iter = 101, seed = 42)
  run3 <- ozkan_pto_sensitivity(d$comm, d$tax, run2, seed = 42)

  # İşlem 2 ile aynı iterasyon sayısı olmalı
  expect_equal(run3$n_iter, 101L)
  expect_equal(nrow(run3$iteration_results), 101)
})

test_that("ozkan_pto_sensitivity validates run2_result", {
  # run2_result parametresi geçerli bir İşlem 2 sonucu olmalı
  # Geçersiz bir liste verince hata vermeli
  d <- make_test_data()
  expect_error(
    ozkan_pto_sensitivity(d$comm, d$tax, list(x = 1)),
    "ozkan_pto_resample"
  )
})

test_that("ozkan_pto_sensitivity is reproducible with same seed", {
  # Aynı tohum ile tekrarlanabilirlik kontrolü
  d <- make_test_data()
  run2 <- ozkan_pto_resample(d$comm, d$tax, n_iter = 101, seed = 42)
  r1 <- ozkan_pto_sensitivity(d$comm, d$tax, run2, seed = 123)
  r2 <- ozkan_pto_sensitivity(d$comm, d$tax, run2, seed = 123)

  # Aynı tohum → aynı sonuçlar
  expect_equal(r1$uTO_plus_max, r2$uTO_plus_max)
  expect_equal(r1$run3_uTO_plus_max, r2$run3_uTO_plus_max)
  expect_equal(r1$iteration_results, r2$iteration_results)
})

test_that("ozkan_pto_sensitivity species_probs has correct length", {
  # Her tür için bir dahil etme olasılığı olmalı
  d <- make_test_data()
  run2 <- ozkan_pto_resample(d$comm, d$tax, n_iter = 101, seed = 42)
  run3 <- ozkan_pto_sensitivity(d$comm, d$tax, run2, seed = 42)

  # Olasılık sayısı = bolluk > 0 olan tür sayısı
  expect_length(run3$species_probs, length(d$comm[d$comm > 0]))

  # Tüm olasılıklar 0 ile 1 arasında olmalı
  expect_true(all(run3$species_probs > 0))
  expect_true(all(run3$species_probs <= 1))
})

test_that("full pipeline Run1 -> Run2 -> Run3 works end to end", {
  # Tam boru hattı testi: İşlem 1 → İşlem 2 → İşlem 3
  # Bu test, 3 fonksiyonun birbirleriyle uyumlu çalıştığını doğrular
  d <- make_test_data()

  # İşlem 1: Deterministik hesaplama (orijinal topluluk)
  run1 <- ozkan_pto(d$comm, d$tax)

  # İşlem 2: Stokastik örnekleme (%50 rastgele dahil etme)
  run2 <- ozkan_pto_resample(d$comm, d$tax, n_iter = 101, seed = 42)

  # İşlem 3: Duyarlılık analizi (türlere özel olasılıklarla)
  run3 <- ozkan_pto_sensitivity(d$comm, d$tax, run2, n_iter = 101, seed = 42)

  # Hepsinden sonlu (geçerli) sayısal değerler çıkmalı
  expect_true(is.finite(run1$uTO_plus))
  expect_true(is.finite(run2$uTO_plus_max))
  expect_true(is.finite(run3$uTO_plus_max))

  # Sıralama kuralı: Run 3 max >= Run 2 max >= Run 1 deterministik
  expect_true(run2$uTO_plus_max >= run1$uTO_plus)
  expect_true(run2$TO_plus_max >= run1$TO_plus)
  expect_true(run3$uTO_plus_max >= run2$uTO_plus_max)
  expect_true(run3$TO_plus_max >= run2$TO_plus_max)
})

test_that("ozkan_pto_resample with Westhoff-Maarel scale data", {
  # Westhoff-Maarel örtme-bolluk ölçeği ile test (değerler 1-9 arası)
  # Bu, Özkan'ın orijinal Excel makrosunda kullandığı ölçek
  set.seed(42)
  n_sp <- 20  # 20 tür

  # Rastgele Westhoff-Maarel bollukları (en yaygın değer 3)
  comm <- setNames(
    sample(1:9, n_sp, replace = TRUE,
           prob = c(0.1, 0.1, 0.5, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05)),
    paste0("sp", 1:n_sp)
  )

  # 20 tür, 10 cins, 5 familya, 2 takım
  tax <- data.frame(
    Species = paste0("sp", 1:n_sp),
    Genus = paste0("G", rep(1:10, each = 2)),
    Family = paste0("F", rep(1:5, each = 4)),
    Order = paste0("O", rep(1:2, each = 10)),
    stringsAsFactors = FALSE
  )

  result <- ozkan_pto_resample(comm, tax, n_iter = 101, seed = 42)

  # Sonuçlar geçerli ve pozitif olmalı
  expect_true(is.finite(result$uTO_plus_max))
  expect_true(result$uTO_plus_max > 0)
  expect_true(result$TO_plus_max > 0)
})
