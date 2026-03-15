# =============================================================================
# test-batch_analysis.R
# batch_analysis() fonksiyonunun testleri
# =============================================================================

# --- Ortak test verisi ---
# Tek alanlı veri (Site sütunu yok)
df_single <- data.frame(
  Species   = c("sp1", "sp2", "sp3", "sp4"),
  Genus     = c("G1", "G1", "G2", "G2"),
  Family    = c("F1", "F1", "F1", "F2"),
  Order     = c("O1", "O1", "O1", "O1"),
  Abundance = c(10, 20, 15, 5),
  stringsAsFactors = FALSE
)

# Çoklu alanlı veri (Site sütunu var)
df_multi <- data.frame(
  Site      = c("A", "A", "A", "A", "B", "B", "B", "B"),
  Species   = c("sp1", "sp2", "sp3", "sp4", "sp1", "sp2", "sp3", "sp4"),
  Genus     = c("G1", "G1", "G2", "G2", "G1", "G1", "G2", "G2"),
  Family    = c("F1", "F1", "F1", "F2", "F1", "F1", "F1", "F2"),
  Order     = c("O1", "O1", "O1", "O1", "O1", "O1", "O1", "O1"),
  Abundance = c(10, 20, 15, 5, 5, 5, 5, 5),
  stringsAsFactors = FALSE
)

# Taksonomik agac (compare_indices testleriyle uyumlu)
tax <- data.frame(
  Species = c("sp1", "sp2", "sp3", "sp4"),
  Genus   = c("G1", "G1", "G2", "G2"),
  Family  = c("F1", "F1", "F1", "F2"),
  Order   = c("O1", "O1", "O1", "O1"),
  stringsAsFactors = FALSE
)


# ---- Test 1: Tek alan — Site sutunu olmadan ----
test_that("tek alan verisi dogru sonuc dondurur", {
  result <- batch_analysis(df_single)

  # Data frame mi?
  expect_true(is.data.frame(result))

  # 1 satir olmali (tek alan)
  expect_equal(nrow(result), 1)

  # 16 sutun: Site + N_Species + 14 indeks (6 klasik + 4 PTO + 4 PTO max)
  expect_equal(ncol(result), 16)

  # Site sutunu "All" olmali (otomatik)
  expect_equal(result$Site, "All")

  # N_Species olmali
  expect_equal(result$N_Species, 4L)

  # Tum indeks sutunlari sayisal olmali
  index_cols <- setdiff(names(result), "Site")
  for (col in index_cols) {
    expect_true(is.numeric(result[[col]]),
                info = paste(col, "sayisal olmali"))
  }
})


# ---- Test 2: Coklu alan — Site sutunu ile ----
test_that("coklu alan verisi dogru sonuc dondurur", {
  result <- batch_analysis(df_multi)

  # 2 satir olmali
  expect_equal(nrow(result), 2)

  # Site isimleri korunmali
  expect_equal(result$Site, c("A", "B"))

  # Her satirin 16 sutunu olmali
  expect_equal(ncol(result), 16)

  # Her site icin N_Species
  expect_equal(result$N_Species, c(4L, 4L))
})


# ---- Test 3: compare_indices ile tutarlilik ----
# batch_analysis sonuclari, ayni veri ile compare_indices sonuclariyla eslesmeli
test_that("tek alan sonuclari compare_indices ile eslesir", {
  batch_result <- batch_analysis(df_single)

  # Ayni veriyle compare_indices calistir
  comm <- c(sp1 = 10, sp2 = 20, sp3 = 15, sp4 = 5)
  ci_result <- compare_indices(comm, tax)

  # Shannon eslesir mi?
  expect_equal(batch_result$Shannon, ci_result$Shannon)

  # Simpson eslesir mi?
  expect_equal(batch_result$Simpson, ci_result$Simpson)

  # Delta eslesir mi?
  expect_equal(batch_result$Delta, ci_result$Delta)

  # Delta_star eslesir mi?
  expect_equal(batch_result$Delta_star, ci_result$Delta_star)

  # AvTD eslesir mi?
  expect_equal(batch_result$AvTD, ci_result$AvTD)

  # VarTD eslesir mi?
  expect_equal(batch_result$VarTD, ci_result$VarTD)

  # pTO bilesenleri
  expect_equal(batch_result$uTO, ci_result$uTO)
  expect_equal(batch_result$TO, ci_result$TO)
  expect_equal(batch_result$uTO_plus, ci_result$uTO_plus)
  expect_equal(batch_result$TO_plus, ci_result$TO_plus)

  # pTO max bilesenleri
  expect_equal(batch_result$uTO_max, ci_result$uTO_max)
  expect_equal(batch_result$TO_max, ci_result$TO_max)
  expect_equal(batch_result$uTO_plus_max, ci_result$uTO_plus_max)
  expect_equal(batch_result$TO_plus_max, ci_result$TO_plus_max)
})


# ---- Test 4: Coklu alan — her alan icin compare_indices ile eslestir ----
test_that("coklu alan sonuclari compare_indices ile eslesir", {
  batch_result <- batch_analysis(df_multi)

  # Site A: comm_A = c(sp1=10, sp2=20, sp3=15, sp4=5)
  comm_A <- c(sp1 = 10, sp2 = 20, sp3 = 15, sp4 = 5)
  ci_A <- compare_indices(comm_A, tax)
  expect_equal(batch_result$Shannon[1], ci_A$Shannon)
  expect_equal(batch_result$Simpson[1], ci_A$Simpson)

  # Site B: comm_B = c(sp1=5, sp2=5, sp3=5, sp4=5)
  comm_B <- c(sp1 = 5, sp2 = 5, sp3 = 5, sp4 = 5)
  ci_B <- compare_indices(comm_B, tax)
  expect_equal(batch_result$Shannon[2], ci_B$Shannon)
  expect_equal(batch_result$Simpson[2], ci_B$Simpson)
})


# ---- Test 5: Site sutunu otomatik algilama ----
# "Site", "Alan", "Plot" isimleri otomatik algilanmali
test_that("Site sutunu otomatik algilanir", {
  # "Site" ismiyle
  result_site <- batch_analysis(df_multi)
  expect_equal(nrow(result_site), 2)

  # "Alan" ismiyle
  df_alan <- df_multi
  names(df_alan)[1] <- "Alan"
  result_alan <- batch_analysis(df_alan)
  expect_equal(nrow(result_alan), 2)

  # "Plot" ismiyle
  df_plot <- df_multi
  names(df_plot)[1] <- "Plot"
  result_plot <- batch_analysis(df_plot)
  expect_equal(nrow(result_plot), 2)
})


# ---- Test 6: site_column parametresi ile ----
test_that("site_column parametresi ile belirtme calisiyor", {
  df_custom <- df_multi
  names(df_custom)[1] <- "Lokasyon"
  result <- batch_analysis(df_custom, site_column = "Lokasyon")
  expect_equal(nrow(result), 2)
  expect_equal(result$Site, c("A", "B"))
})


# ---- Test 7: Bos Site sutunu — tek alan gibi calismali ----
test_that("bos Site sutunu tek alan olarak calisir", {
  df_empty_site <- df_single
  df_empty_site$Site <- ""
  result <- batch_analysis(df_empty_site)
  expect_equal(nrow(result), 1)
  expect_equal(result$Site, "All")
})


# ---- Test 8: NA Site sutunu — tek alan gibi calismali ----
test_that("NA Site sutunu tek alan olarak calisir", {
  df_na_site <- df_single
  df_na_site$Site <- NA
  result <- batch_analysis(df_na_site)
  expect_equal(nrow(result), 1)
  expect_equal(result$Site, "All")
})


# ---- Test 9: Case-insensitive abundance sutunu ----
test_that("abundance sutunu case-insensitive eslesiyor", {
  df_lower <- df_single
  names(df_lower)[names(df_lower) == "Abundance"] <- "abundance"
  result <- batch_analysis(df_lower)
  expect_equal(nrow(result), 1)
})


# ---- Test 10: Hatali girdi kontrolleri ----
test_that("hatali girdi durumunda hata firlatir", {
  # data.frame degil
  expect_error(batch_analysis("not_a_df"), "data.*must be a data frame")

  # Bos data.frame
  expect_error(batch_analysis(data.frame()), "no rows")

  # Abundance sutunu yok
  df_no_abd <- df_single[, -5]
  expect_error(batch_analysis(df_no_abd), "Abundance.*not found")

  # Taksonomik sutunlar yetersiz
  df_no_tax <- data.frame(
    Species   = c("sp1", "sp2"),
    Abundance = c(10, 20),
    stringsAsFactors = FALSE
  )
  expect_error(batch_analysis(df_no_tax), "auto-detect")

  # Belirtilen site_column bulunamiyor
  expect_error(batch_analysis(df_single, site_column = "Nonexistent"),
               "not found")
})


# ---- Test 11: tax_columns parametresi ile ----
test_that("tax_columns parametresi ile calisiyor", {
  df_custom_names <- data.frame(
    Tur     = c("sp1", "sp2", "sp3", "sp4"),
    Cins    = c("G1", "G1", "G2", "G2"),
    Familya = c("F1", "F1", "F1", "F2"),
    Takim   = c("O1", "O1", "O1", "O1"),
    Bolluk  = c(10, 20, 15, 5),
    stringsAsFactors = FALSE
  )
  result <- batch_analysis(df_custom_names,
                           tax_columns = c("Tur", "Cins", "Familya", "Takim"),
                           abundance_column = "Bolluk")
  expect_equal(nrow(result), 1)
  expect_true(is.numeric(result$Shannon))
})


# ---- Test 12: Esit dagilimda cesitlilik daha yuksek ----
test_that("esit dagilimda cesitlilik daha yuksek", {
  result <- batch_analysis(df_multi)

  # Site B (esit dagilim) Shannon daha yuksek olmali
  expect_gt(result$Shannon[result$Site == "B"],
            result$Shannon[result$Site == "A"])
})


# ---- Test 13: Turkce alan ismiyle otomatik algilama ----
test_that("alan ismiyle otomatik algilama calisiyor", {
  df_alan <- df_multi
  names(df_alan)[1] <- "alan"  # kucuk harfle
  result <- batch_analysis(df_alan)
  expect_equal(nrow(result), 2)
})


# ---- Test 14: 3 veya daha fazla site ----
test_that("3 site ile calisiyor", {
  df_three <- rbind(
    data.frame(Site = "X", df_single, stringsAsFactors = FALSE),
    data.frame(Site = "Y", df_single, stringsAsFactors = FALSE),
    data.frame(Site = "Z", df_single, stringsAsFactors = FALSE)
  )
  # Y sitesinde farkli abundance
  df_three$Abundance[df_three$Site == "Y"] <- c(5, 5, 5, 5)
  df_three$Abundance[df_three$Site == "Z"] <- c(50, 1, 1, 1)

  result <- batch_analysis(df_three)
  expect_equal(nrow(result), 3)
  expect_equal(result$Site, c("X", "Y", "Z"))
})


# ---- Test 15: Sifir abundance turler filtrelenmeli ----
test_that("sifir abundance turler filtreleniyor", {
  df_zeros <- data.frame(
    Species   = c("sp1", "sp2", "sp3", "sp4", "sp5"),
    Genus     = c("G1", "G1", "G2", "G2", "G3"),
    Family    = c("F1", "F1", "F1", "F2", "F2"),
    Order     = c("O1", "O1", "O1", "O1", "O1"),
    Abundance = c(10, 20, 15, 5, 0),
    stringsAsFactors = FALSE
  )
  # Hata firlatmamali — sp5 (abundance = 0) atlanmali
  result <- batch_analysis(df_zeros)
  expect_equal(nrow(result), 1)
  expect_true(is.numeric(result$Shannon))
})
