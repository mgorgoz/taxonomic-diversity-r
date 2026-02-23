# taxdiv Ornek Ciktilar ve Yontem Rehberi

Bu klasor, `taxdiv` paketinin urettigi ornek gorselleri ve her yontemin
ne anlama geldigini aciklar. 10 agac turlu hayali bir Akdeniz ormani
toplulugu kullanilmistir. Taksonomik hiyerarsi 7 seviyedir:
**Tur > Cins > Familya > Takim > Sinif > Bolum > Alem**

---

## Icindekiler

1. [Ornek Topluluklar](#ornek-topluluklar)
2. [Shannon Indeksi](#1-shannon-indeksi)
3. [Simpson Indeksi](#2-simpson-indeksi)
4. [Pielou Esitlik Indeksi](#3-pielou-esitlik-indeksi)
5. [Taksonomik Mesafe Matrisi](#4-taksonomik-mesafe-matrisi)
6. [Clarke & Warwick Indeksleri](#5-clarke--warwick-indeksleri)
7. [Deng Entropisi](#6-deng-entropisi)
8. [Ozkan pTO](#7-ozkan-pto-taksonomik-cesitlilik)
9. [Rarefaction (Seyreltme Egrisi)](#8-rarefaction-seyreltme-egrisi)
10. [Gorsellestirmeler](#9-gorsellestirmeler)
11. [Kaynaklar](#kaynaklar)

---

## Ornek Topluluklar

Orneklerde iki farkli senaryo kullanilir. Amac, ayni turler oldugu halde
bolluk dagilimi degisince indekslerin nasil tepki verdigini gostermektir.

### "Diverse" (Dengeli Dagilim)

Butun turlerin bollugu birbirine yakin. Kimse cok baskin degil.

```
Pinus_brutia = 30, Quercus_coccifera = 25, Fagus_orientalis = 20,
Quercus_infectoria = 18, Cedrus_libani = 15, Pinus_nigra = 12,
Carpinus_betulus = 10, Juniperus_excelsa = 8, Abies_cilicica = 7,
Juniperus_oxycedrus = 5
```

### "Dominant" (Esitsiz Dagilim)

Bir tur (Quercus coccifera = 80 birey) diger herkesi eziyor.
Geri kalanlar 1-5 birey arasinda.

```
Quercus_coccifera = 80, Quercus_infectoria = 5, Pinus_brutia = 3,
Cedrus_libani = 3, Pinus_nigra = 2, Juniperus_excelsa = 2,
Fagus_orientalis = 2, Juniperus_oxycedrus = 1, Abies_cilicica = 1,
Carpinus_betulus = 1
```

> **Not:** "Diverse" ve "Dominant" resmi ekoloji terimleri degildir.
> Bolluk dagiliminin esitligini (evenness) betimleyen etiketlerdir.

---

## 1. Shannon Indeksi

**Fonksiyon:** `shannon(community)`

**Ne olcuyor?** Topluluktaki belirsizligi. Rastgele bir birey secsen,
hangi turden olacagini ne kadar zor tahmin edersin?

**Formul:**

```
H' = -sum(p_i * ln(p_i))
```

`p_i` = i. turun oransal bolluğu (o turun birey sayisi / toplam birey sayisi)

**Nasil yorumlanir?**
- H' = 0 ise: Tek bir tur var. Tahmin etmek cok kolay.
- H' yuksekse: Cok tur var ve bolluklar dengeli. Tahmin etmek zor = yuksek cesitlilik.

**Ornek:**
- 3 tur esit bollukta: H' = ln(3) = 1.099
- 10 tur esit bollukta: H' = ln(10) = 2.303
- 10 tur ama biri baskin: H' < 2.303 (baskin tur belirsizligi azaltir)

**Kaynak:** Shannon, C.E. (1948). A mathematical theory of communication.
*Bell System Technical Journal*, 27, 379-423.

---

## 2. Simpson Indeksi

**Fonksiyon:** `simpson(community)`

**Ne olcuyor?** Rastgele secilen iki bireyin farkli turlerden olma olasiligini.

**Formul (Gini-Simpson):**

```
1 - D = 1 - sum(p_i^2)
```

**Nasil yorumlanir?**
- 0'a yakinsa: Bir tur her yerde. Iki birey secsen buyuk ihtimal ayni tur.
- 1'e yakinsa: Turler dengeli dagilmis. Iki birey secsen buyuk ihtimal farkli turler.

**Shannon'dan farki ne?** Simpson baskin turlere daha duyarli.
Shannon nadir turlere daha fazla agirlik verir. Ikisi birlikte kullanilir.

**Kaynak:** Simpson, E.H. (1949). Measurement of diversity.
*Nature*, 163, 688.

---

## 3. Pielou Esitlik Indeksi

**Fonksiyon:** Dogrudan bir fonksiyon yok, ama `shannon()` ile hesaplanir.

**Ne olcuyor?** Bolluklar turler arasinda ne kadar esit dagilmis?

**Formul:**

```
J = H' / ln(S)
```

`H'` = Shannon indeksi, `S` = tur sayisi

**Nasil yorumlanir?**
- J = 1: Butun turler esit bollukta. Mukemmel esitlik.
- J = 0: Tek bir tur topluluga tamamen hakim.
- J = 0.8: Gayet dengeli bir topluluk.
- J = 0.3: Ciddi baskinlik var.

**Kaynak:** Pielou, E.C. (1966). The measurement of diversity in
different types of biological collections. *Journal of Theoretical
Biology*, 13, 131-144.

---

## 4. Taksonomik Mesafe Matrisi

**Fonksiyon:** `tax_distance_matrix(tax_tree)`

**Ne olcuyor?** Iki tur arasindaki taksonomik uzakligi. Yani "bu iki tur
taksonomik olarak ne kadar farkli?"

**Nasil hesaplaniyor?** Iki turun taksonomik agacta ilk birlestigi seviyeye
bakilir:

```
Ornek: Pinus nigra ve Quercus cerris
- Ayni cins mi? Hayir (Pinus vs Quercus) -> 1 adim
- Ayni familya mi? Hayir (Pinaceae vs Fagaceae) -> 2 adim
- Ayni takim mi? Hayir (Pinales vs Fagales) -> 3 adim
- Ayni sinif mi? Hayir (Pinopsida vs Magnoliopsida) -> 4 adim
- Ayni bolum mu? Hayir (Pinophyta vs Magnoliophyta) -> 5 adim
- Ayni alem mi? Evet (Plantae) -> Mesafe = 5
```

```
Ornek: Pinus nigra ve Abies nordmanniana
- Ayni cins mi? Hayir -> 1 adim
- Ayni familya mi? Evet (Pinaceae) -> Mesafe = 1
```

Yani ne kadar uzakta birlesiyorlarsa, taksonomik mesafe o kadar buyuk.

**Ilgili gorsel:** `heatmap.png`

**Kaynak:** Clarke, K.R. & Warwick, R.M. (1998). A taxonomic distinctness
index and its statistical properties. *Journal of Applied Ecology*,
35, 523-531.

---

## 5. Clarke & Warwick Indeksleri

### 5a. Delta (Ortalama Taksonomik Mesafe)

**Fonksiyon:** `delta(community, tax_tree)`

**Ne olcuyor?** Topluluktan rastgele secilen iki birey arasindaki
beklenen taksonomik mesafe. Bolluk bilgisini kullanir.

**Formul:**

```
Delta = [sum_i sum_j (w_ij * n_i * n_j)] / [N * (N-1) / 2]
```

`w_ij` = i ve j turleri arasindaki taksonomik mesafe,
`n_i` = i. turun birey sayisi, `N` = toplam birey sayisi

**Yorum:** Yuksek Delta = turler birbirinden taksonomik olarak uzak
VE bolluklar dengeli.

### 5b. Delta* (Taksonomik Mesafe - Bolluk Normalizeli)

**Fonksiyon:** `delta_star(community, tax_tree)`

**Ne olcuyor?** Delta ile ayni ama taksonomik yapiyi bolluk
dagiliminin etkisinden ayirmaya calisir.

### 5c. AvTD / Delta+ (Ortalama Taksonomik Farklilik)

**Fonksiyon:** `avtd(species_names, tax_tree)`

**Ne olcuyor?** Sadece tur listesine (varlik/yokluk) dayali ortalama
taksonomik mesafe. Bolluk bilgisini KULLANMAZ.

**Neden onemli?** Farkli orneklem buyuklugundeki alanlari karsilastirabilirsin.
Cunku bolluktan bagimsiz.

**Formul:**

```
Delta+ = [sum_i<j w_ij] / [S * (S-1) / 2]
```

`S` = tur sayisi, `w_ij` = turler arasi taksonomik mesafe

**Yorum:** Yuksek AvTD = turler taksonomik olarak birbirine uzak.
Dusuk AvTD = turler birbirine yakin (mesela hep ayni familyadan).

### 5d. VarTD / Lambda+ (Taksonomik Farklilik Varyasyonu)

**Fonksiyon:** `vartd(species_names, tax_tree)`

**Ne olcuyor?** Taksonomik mesafelerin ne kadar degisken oldugu.
Bazi tur ciftleri cok yakin, bazilari cok uzak mi?

**Yorum:** Yuksek VarTD = taksonomik yapi duzensiz (bazi turler
cok yakin, bazilari cok uzak). Dusuk VarTD = turler arasi mesafeler
homojen.

**Kaynak (5a-5d icin):** Clarke, K.R. & Warwick, R.M. (1998). A
taxonomic distinctness index and its statistical properties. *Journal
of Applied Ecology*, 35, 523-531.

Clarke, K.R. & Warwick, R.M. (2001). A further biodiversity index
applicable to species lists: variation in taxonomic distinctness.
*Marine Ecology Progress Series*, 216, 265-278.

---

## 6. Deng Entropisi

**Fonksiyon:** `deng_entropy_level(counts, group_sizes)`

**Ne olcuyor?** Klasik Shannon entropisi "her tur bireydir" der.
Deng entropisi ise "turler gruplara ait, gruplar da ust gruplara ait"
der. Yani taksonomik hiyerarsiyi hesaba katar.

**Temel fikir:** Shannon entropisinde her olay (tur) atomiktir.
Deng entropisinde olaylar (turler) kume yapisi icerir.
Dempster-Shafer kanit teorisi cercevesinde calisir.

**Formul:**

```
Ed = -sum(m_i * ln(m_i / (2^|F_i| - 1)))
```

`m_i` = gruptaki turlerin oransal agirligi,
`|F_i|` = gruptaki tur sayisi

**Yorum:** Eger bir familyada 5 tur varsa ve baska familyada 1 tur varsa,
Deng entropisi bu kume yapisini dikkate alir. Shannon bunu goremez.

**Kaynak:** Deng, Y. (2016). Deng entropy. *Chaos, Solitons & Fractals*,
91, 549-553.

---

## 7. Ozkan pTO (Taksonomik Cesitlilik)

**Fonksiyon:** `ozkan_pto(community, tax_tree)` ve `pto_components(community, tax_tree)`

**Ne olcuyor?** Deng entropisini kullanarak taksonomik cesitliligi ve
taksonomik mesafeyi birlikte olcer. 4 bilesen uretir.

### 4 Bilesen:

| Bilesen | Bolluk kullanir mi? | Ne olcuyor? |
|---------|---------------------|-------------|
| **uTO** | Evet (Islem 1+2+3) | Agirliksiz taksonomik cesitlilik |
| **TO** | Evet (Islem 1+2+3) | Agirlikli taksonomik cesitlilik |
| **uTO+** | Hayir (sadece Islem 1) | Agirliksiz taksonomik mesafe |
| **TO+** | Hayir (sadece Islem 1) | Agirlikli taksonomik mesafe |

### Islem 1 (Run 1 - Deterministik):

Tum turler dahil, hepsine esit agirlik verilir. Deng entropisi her
taksonomik seviyede hesaplanir. Sonuc: temel deger.

### Islem 2 (Run 2 - Stokastik Yeniden Ornekleme):

Her tur %50 olasilikla dahil edilir veya cikarilir. Bu islem 101+ kez
tekrarlanir. Her seferinde pTO yeniden hesaplanir. Maksimum degerler
alinir.

**Fonksiyon:** `ozkan_pto_resample(community, tax_tree, n_iter = 101)`

### Islem 3 (Run 3 - Duyarlilik Analizi):

Islem 2'de bulunan maksimum degerleri referans alarak, %5'lik
hassasiyet esigi ile kararliligi test eder. Kac iterasyonun referans
degere yakin olduguna bakar.

**Fonksiyon:** `ozkan_pto_sensitivity(community, tax_tree, n_iter = 101)`

### Neden TO+ sadece Islem 1?

Taksonomik mesafe (TO+) turlerin birbirine ne kadar uzak oldugunu olcer.
Bu, bolluklarla ilgisizdir. 5 bireyi olan Pinus ile 500 bireyi olan
Pinus hala ayni familyadadir. Bu yuzden dilim prosedurune (Islem 2/3)
ihtiyac duymaz.

Taksonomik cesitlilik (TO) ise bollugu da hesaba katmalidir. Dilim
proseduru, bol turlerin daha cok dilimde hayatta kalmasini saglayarak
bolluğu sisteme dolayli olarak sokar.

**Kaynak:** Ozkan, K. (2018). Yeni bir taksonomik cesitlilik olcusu
onerisi. *Turkish Journal of Forestry*, 19(4), 336-346.
DOI: 10.18182/tjf.441061

---

## 8. Rarefaction (Seyreltme Egrisi)

**Fonksiyon:** `rarefaction_taxonomic(community, tax_tree, index = "shannon")`

**Ne olcuyor?** "Orneklemem yeterli mi?" sorusuna cevap verir.

### Problem:

Diyelim ki iki alanda calisan var:
- Alan A: 500 birey sayilmis, 30 tur bulunmus
- Alan B: 100 birey sayilmis, 15 tur bulunmus

Alan A daha cesitli mi? Belki. Ama belki de sadece daha cok birey
sayildigi icin daha cok tur bulunmustur. Adil karsilastirma yapmak
icin ayni orneklem buyuklugune "seyreltmek" gerekir.

### Nasil calisiyor?

1. Toplam N bireyden rastgele n tanesini sec (n < N)
2. Secilen bireylerdeki tur bolluklarini say
3. Secilen cesitlilik indeksini hesapla
4. Bunu 100-200 kez tekrarla, ortalamasini al
5. Farkli n degerleri icin tekrarla (mesela 10, 20, 50, 100, ...)

Sonucta bir egri elde edersin:

```
Cesitlilik
    |        ___________  <-- egri duzlesti = orneklemen yeterli
    |       /
    |      /
    |     /    <-- hizli artis = henuz yetersiz
    |    /
    |   /
    |--/
    +--------------------- Orneklem Buyuklugu
```

### Guven araligi (CI) ne?

Egrinin etrafindaki gri/mavi bant. "Gercek deger %95 olasilikla bu
araliktadir" demek. Bant darsa = sonuc guvenilir. Genisse = belirsiz.

### Desteklenen indeksler:

| index | Aciklama |
|-------|----------|
| `"species"` | Tur zenginligi (S) |
| `"shannon"` | Shannon H' |
| `"simpson"` | Gini-Simpson (1-D) |
| `"uTO"` | Ozkan agirliksiz taksonomik cesitlilik |
| `"TO"` | Ozkan agirlikli taksonomik cesitlilik |
| `"uTO_plus"` | Ozkan agirliksiz taksonomik mesafe |
| `"TO_plus"` | Ozkan agirlikli taksonomik mesafe |
| `"avtd"` | Clarke & Warwick AvTD |

**Gorsellestirme:** `plot_rarefaction(rare_result)`

**Ilgili gorseller:**
- `rarefaction_shannon.png` - Shannon seyreltme egrisi
- `rarefaction_species.png` - Tur zenginligi seyreltme egrisi
- `rarefaction_uTO.png` - uTO seyreltme egrisi
- `rarefaction_avtd.png` - AvTD seyreltme egrisi

**Kaynaklar:**

Gotelli, N.J. & Colwell, R.K. (2001). Quantifying biodiversity:
procedures and pitfalls in the measurement and comparison of species
richness. *Ecology Letters*, 4, 379-391.

Hurlbert, S.H. (1971). The nonconcept of species diversity: a critique
and alternative parameters. *Ecology*, 52, 577-586.

Sanders, H.L. (1968). Marine benthic diversity: a comparative study.
*The American Naturalist*, 102, 243-282.

---

## 9. Gorsellestirmeler

### 9a. Dendogram (Taksonomik Agac)

**Fonksiyon:** `plot_taxonomic_tree(tax_tree, community)`

Turleri taksonomik mesafelerine gore bir agac seklinde gosterir.
Yakin turler yakin dallarda, uzak turler uzak dallarda yer alir.
Renklendirme taksonomik seviyeye gore yapilir (familya, sinif vb.).

**Ilgili gorseller:**
- `dendrogram_family.png` - Familyaya gore renklendirilmis
- `dendrogram_class.png` - Sinifa gore renklendirilmis
- `dendrogram_no_abundance.png` - Bolluksuz, Takima gore renkli

### 9b. Heatmap (Isi Haritasi)

**Fonksiyon:** `plot_heatmap(tax_tree)`

Taksonomik mesafe matrisini renkli kareler olarak gosterir.
Koyu renk = iki tur birbirine yakin. Acik renk = uzak.
Turler hclust ile siralanir, boylece benzer turler yan yana gelir.

**Ilgili gorsel:** `heatmap.png`

### 9c. Karsilastirma Grafigi

**Fonksiyon:** `compare_indices(communities, tax_tree, plot = TRUE)`

Birden fazla toplulugu tum indeksler uzerinden karsilastirir.
Gruplu bar plot olarak gosterir.

**Onemli gozlem:** Bolluga bagimli indeksler (Shannon, Simpson, Delta,
uTO, TO) Diverse ve Dominant arasinda degisirken, bolluktan bagimsiz
indeksler (AvTD, VarTD, uTO+, TO+) ayni kalir. Cunku ikinci gruptakiler
sadece tur listesine bakar, bolluga degil.

**Ilgili gorsel:** `compare_indices_barplot.png`

### 9d. Iterasyon Grafigi (Islem 2/3)

**Fonksiyon:** `plot_iteration(resample_result)`

Islem 2'deki her tekrarin (iterasyon) urettigi pTO degerini gosterir.
Kirmizi kesikli cizgi = Islem 1 (deterministik) sonucu.
Mavi kesikli cizgi = maksimum deger.

**Ilgili gorseller:**
- `iteration_run2_TO.png` - TO bileseninin iterasyonlari
- `iteration_run2_uTO_plus.png` - uTO+ bileseninin iterasyonlari

### 9e. Bubble (Balon) Grafigi

**Fonksiyon:** `plot_bubble(community, tax_tree)`

Her tur bir balon. X ekseni = bolluk, Y ekseni = ortalama taksonomik
mesafe, balon buyuklugu = turun topluluga katkisi (bolluk x mesafe).

Ust sag kosedeki buyuk balonlar = hem bol hem taksonomik olarak farkli
turler. Bunlar toplulugun "onemli" turleridir.

**Ilgili gorsel:** `bubble_family.png`

### 9f. Radar (Orumcek) Grafigi

**Fonksiyon:** `plot_radar(communities, tax_tree)`

Birden fazla toplulugu tum indeksler uzerinden orumcek agi seklinde
karsilastirir. Her eksen bir indeks. Deger buyukse o eksen disari
cikar. Topluluklarin genel profili bir bakista gorulur.

**Ilgili gorsel:** `radar_comparison.png`

### 9g. Rarefaction (Seyreltme) Grafigi

**Fonksiyon:** `plot_rarefaction(rare_result)`

Seyreltme egrisini guven araliklari ile birlikte gosterir.
Kirmizi kesikli cizgi = toplam orneklem buyuklugu (N).
Mavi bant = %95 guven araligi.

**Ilgili gorseller:**
- `rarefaction_shannon.png`
- `rarefaction_species.png`
- `rarefaction_uTO.png`
- `rarefaction_avtd.png`

---

## Dosya Listesi

| Dosya | Fonksiyon | Aciklama |
|-------|-----------|----------|
| `dendrogram_family.png` | `plot_taxonomic_tree()` | Familyaya gore renkli dendogram |
| `dendrogram_class.png` | `plot_taxonomic_tree()` | Sinifa gore renkli dendogram |
| `dendrogram_no_abundance.png` | `plot_taxonomic_tree()` | Bolluksuz dendogram |
| `compare_indices_barplot.png` | `compare_indices()` | 10 indeks karsilastirma grafigi |
| `heatmap.png` | `plot_heatmap()` | Taksonomik mesafe isi haritasi |
| `iteration_run2_TO.png` | `plot_iteration()` | TO iterasyon grafigi |
| `iteration_run2_uTO_plus.png` | `plot_iteration()` | uTO+ iterasyon grafigi |
| `bubble_family.png` | `plot_bubble()` | Balon grafigi |
| `radar_comparison.png` | `plot_radar()` | Orumcek grafigi |
| `rarefaction_shannon.png` | `plot_rarefaction()` | Shannon seyreltme egrisi |
| `rarefaction_species.png` | `plot_rarefaction()` | Tur zenginligi seyreltme egrisi |
| `rarefaction_uTO.png` | `plot_rarefaction()` | uTO seyreltme egrisi |
| `rarefaction_avtd.png` | `plot_rarefaction()` | AvTD seyreltme egrisi |

---

## Kaynaklar

### Klasik Cesitlilik Indeksleri
- Shannon, C.E. (1948). A mathematical theory of communication.
  *Bell System Technical Journal*, 27, 379-423.
- Simpson, E.H. (1949). Measurement of diversity. *Nature*, 163, 688.
- Pielou, E.C. (1966). The measurement of diversity in different types
  of biological collections. *Journal of Theoretical Biology*, 13, 131-144.

### Taksonomik Farklilik (Clarke & Warwick)
- Clarke, K.R. & Warwick, R.M. (1998). A taxonomic distinctness index
  and its statistical properties. *Journal of Applied Ecology*, 35, 523-531.
- Clarke, K.R. & Warwick, R.M. (2001). A further biodiversity index
  applicable to species lists: variation in taxonomic distinctness.
  *Marine Ecology Progress Series*, 216, 265-278.

### Deng Entropisi ve Dempster-Shafer Teorisi
- Deng, Y. (2016). Deng entropy. *Chaos, Solitons & Fractals*, 91, 549-553.
- Dempster, A.P. (1967). Upper and lower probabilities induced by a
  multivalued mapping. *The Annals of Mathematical Statistics*, 38, 325-339.
- Shafer, G. (1976). *A Mathematical Theory of Evidence*. Princeton
  University Press.

### Ozkan pTO Indeksi
- Ozkan, K. (2018). Yeni bir taksonomik cesitlilik olcusu onerisi.
  *Turkish Journal of Forestry*, 19(4), 336-346. DOI: 10.18182/tjf.441061

### Rarefaction (Seyreltme)
- Gotelli, N.J. & Colwell, R.K. (2001). Quantifying biodiversity:
  procedures and pitfalls in the measurement and comparison of species
  richness. *Ecology Letters*, 4, 379-391.
- Hurlbert, S.H. (1971). The nonconcept of species diversity: a critique
  and alternative parameters. *Ecology*, 52, 577-586.
- Sanders, H.L. (1968). Marine benthic diversity: a comparative study.
  *The American Naturalist*, 102, 243-282.

### Genel Ekoloji ve Biyocesitlilik
- Magurran, A.E. (2004). *Measuring Biological Diversity*. Blackwell
  Publishing, Oxford.
- Whittaker, R.H. (1972). Evolution and measurement of species diversity.
  *Taxon*, 21, 213-251.
