# =============================================================================
# test-compare_indices.R
# compare_indices() fonksiyonunun testleri
# =============================================================================

# --- Ortak test verisi ---
# 4 türlü basit taksonomik ağaç
tax <- data.frame(
  Species = c("sp1", "sp2", "sp3", "sp4"),
  Genus   = c("G1", "G1", "G2", "G2"),
  Family  = c("F1", "F1", "F1", "F2"),
  Order   = c("O1", "O1", "O1", "O1"),
  stringsAsFactors = FALSE
)

# İki farklı topluluk: biri eşit dağılımlı, biri baskın türlü
comm_equal   <- c(sp1 = 5, sp2 = 5, sp3 = 5, sp4 = 5)
comm_unequal <- c(sp1 = 50, sp2 = 2, sp3 = 1, sp4 = 1)


# ---- Test 1: Tek topluluk girişi doğru data.frame döndürmeli ----
# Tek bir vektör verildiğinde "Community" ismiyle sarmalanır
test_that("tek topluluk vektoru dogru tablo dondurur", {
  result <- compare_indices(comm_equal, tax)

  # Data frame mi?
  expect_true(is.data.frame(result))

  # 1 satır olmalı (tek topluluk)
  expect_equal(nrow(result), 1)

  # 16 sütun: Community + N_Species + 14 indeks (6 klasik + 4 PTO + 4 PTO max)
  expect_equal(ncol(result), 16)

  # N_Species kontrol
  expect_equal(result$N_Species, 4L)

  # Community sütunu "Community" olmalı
  expect_equal(result$Community, "Community")

  # Tüm indeks sütunları sayısal olmalı
  index_cols <- setdiff(names(result), "Community")
  for (col in index_cols) {
    expect_true(is.numeric(result[[col]]),
                info = paste(col, "sayisal olmali"))
  }
})


# ---- Test 2: Çoklu topluluk girişi ----
# Liste olarak verilen birden fazla topluluk
test_that("coklu topluluk listesi dogru tablo dondurur", {
  comm_list <- list(
    Esit    = comm_equal,
    Baskin  = comm_unequal
  )
  result <- compare_indices(comm_list, tax)

  # 2 satır olmalı
  expect_equal(nrow(result), 2)

  # Topluluk isimleri korunmalı
  expect_equal(result$Community, c("Esit", "Baskin"))
})


# ---- Test 3: Bilinen değerlerle doğrulama ----
# Tek tek hesaplanan değerlerle compare_indices sonuçları eşleşmeli
test_that("degerler tek tek hesaplanan indekslerle eslesiyor", {
  result <- compare_indices(comm_equal, tax)

  # Shannon
  expect_equal(result$Shannon, round(shannon(comm_equal), 6))

  # Simpson
  expect_equal(result$Simpson, round(simpson(comm_equal), 6))

  # Delta
  expect_equal(result$Delta, round(delta(comm_equal, tax), 6))

  # Delta*
  expect_equal(result$Delta_star, round(delta_star(comm_equal, tax), 6))

  # AvTD
  sp <- names(comm_equal[comm_equal > 0])
  expect_equal(result$AvTD, round(avtd(sp, tax), 6))

  # VarTD
  expect_equal(result$VarTD, round(vartd(sp, tax), 6))

  # pTO bileşenleri (unname ile isim farkını kaldır)
  pto <- pto_components(comm_equal, tax)
  expect_equal(result$uTO, round(unname(pto["uTO"]), 6))
  expect_equal(result$TO, round(unname(pto["TO"]), 6))
  expect_equal(result$uTO_plus, round(unname(pto["uTO_plus"]), 6))
  expect_equal(result$TO_plus, round(unname(pto["TO_plus"]), 6))

  # pTO max bileşenleri
  expect_equal(result$uTO_max, round(unname(pto["uTO_max"]), 6))
  expect_equal(result$TO_max, round(unname(pto["TO_max"]), 6))
  expect_equal(result$uTO_plus_max, round(unname(pto["uTO_plus_max"]), 6))
  expect_equal(result$TO_plus_max, round(unname(pto["TO_plus_max"]), 6))
})


# ---- Test 4: Eşit dağılım vs baskın tür karşılaştırması ----
# Eşit dağılımlı toplulukta Shannon daha yüksek olmalı
test_that("esit dagilimda Shannon ve Simpson daha yuksek", {
  comm_list <- list(Esit = comm_equal, Baskin = comm_unequal)
  result <- compare_indices(comm_list, tax)

  esit  <- result[result$Community == "Esit", ]
  baskin <- result[result$Community == "Baskin", ]

  # Eşit dağılımda Shannon daha yüksek
  expect_gt(esit$Shannon, baskin$Shannon)

  # Eşit dağılımda Simpson daha yüksek
  expect_gt(esit$Simpson, baskin$Simpson)
})


# ---- Test 5: plot = TRUE ile liste döndürmeli ----
# ggplot2 yüklüyse hem tablo hem grafik döner
test_that("plot = TRUE ile liste dondurur", {
  skip_if_not_installed("ggplot2")

  result <- compare_indices(comm_equal, tax, plot = TRUE)

  # Liste olmalı
  expect_true(is.list(result))
  expect_true("table" %in% names(result))
  expect_true("plot" %in% names(result))

  # table bir data.frame olmalı
  expect_true(is.data.frame(result$table))

  # plot bir ggplot nesnesi olmalı
  expect_s3_class(result$plot, "ggplot")
})


# ---- Test 6: Hatalı girdi kontrolü ----
# Yanlış veri türü verildiğinde hata fırlatmalı
test_that("hatali girdi durumunda hata firlatir", {
  # tax_tree data.frame değilse
  expect_error(compare_indices(comm_equal, "not_a_df"),
               "tax_tree")

  # community sayısal değilse
  expect_error(compare_indices("not_numeric", tax),
               "must be")

  # community isimsizse
  expect_error(compare_indices(c(1, 2, 3), tax),
               "must be")
})
