# Rakip Paket Analizi / Competitive Analysis

**Tarih:** Mart 2026
**Soru:** R ekosisteminde taksonomik cesitlilik hesaplayan baska paketler var mi?
taxdiv ne fark yaratiyor?

---

## Rakip Paketler

### 1. vegan (en yaygin ekoloji paketi)

**Ne yapar:** Topluluk ekolojisi icin genel amacli paket. Ordinasyon,
permutasyon testleri, cesitlilik indeksleri.

**Taksonomik cesitlilik fonksiyonlari:**
- `taxondive()` — Delta, Delta*, Delta+, Lambda+ hesaplar
- `taxa2dist()` — Taksonomik siniflandirmayi mesafe matrisine cevirir

**Eksikleri:**
- Deng entropisi YOK
- Ozkan pTO YOK
- Dilim proseduru (slicing) YOK
- Gorsellesirme sinirli (temel R grafikleri)
- Stokastik yeniden ornekleme (Islem 2/3) YOK

**taxdiv farki:** vegan sadece Clarke & Warwick indekslerini hesaplar.
taxdiv bunlara ek olarak Deng entropisi tabanli olculeri, dilim
prosedurunu ve tam Islem 1/2/3 pipeline'ini sunar.

---

### 2. hillR (Hill sayilari ile cesitlilik)

**Ne yapar:** Taksonomik, fonksiyonel ve filogenetik cesitliligi Hill
sayilari cercevesinde hesaplar (Chao, Chiu & Jost 2014).

**Fonksiyonlari:**
- `hill_taxa()` — Alfa taksonomik cesitlilik (Hill sayilari)
- `hill_taxa_parti()` — Cesitlilik parcalama (alfa/beta/gamma)
- q parametresi: q=0 (zenginlik), q=1 (Shannon), q=2 (Simpson)

**Eksikleri:**
- Taksonomik hiyerarsiyi KULLANMAZ (sadece tur bolluklari)
- Deng entropisi YOK
- Taksonomik mesafe matrisi YOK
- Clarke & Warwick indeksleri YOK

**taxdiv farki:** hillR tur bolluklarina dayali Hill sayilari hesaplar.
Taksonomik yapi bilgisini hic kullanmaz. taxdiv ise taksonomik
hiyerarsiyi temel alan olculer sunar.

---

### 3. entropart (Entropi tabanli cesitlilik)

**Ne yapar:** HCDT (Tsallis) entropisi ve benzerlik tabanli entropi ile
cesitlilik hesaplar. Alfa/beta/gamma parcalama destekler.

**Fonksiyonlari:**
- `Tsallis()`, `Shannon()`, `Simpson()` — Entropi hesaplamalari
- `Diversity()` — Hill sayilari
- `Dqz()` — Benzerlik matrisli cesitlilik
- `PhyloDiversity()` — Filogenetik entropi
- Bias-corrected (yanlilik duzeltmeli) tahminleyiciler

**Eksikleri:**
- Deng entropisi YOK (Tsallis/HCDT cercevesi kullanir)
- Ozkan pTO YOK
- Clarke & Warwick taksonomik indeksleri YOK
- Taksonomik hiyerarsi tabanli hesaplama YOK

**taxdiv farki:** entropart sofistike bir entropi paketidir ama
Dempster-Shafer kanit teorisi cercevesinde calismaz. Deng entropisi
ile Tsallis entropisi farkli matematiksel temellere dayanir.
taxdiv, Deng entropisinin taksonomik gruplardaki fokal eleman
boyutunu (|F_i|) kullanma avantajini sunar.

---

### 4. adiv (Cesitlilik Analizi)

**Ne yapar:** Tur, fonksiyonel ve filogenetik cesitlilik hesaplar.
Tur orijinalitesi ve benzemezlik olculeri sunar.

**Fonksiyonlari:**
- `EqRao()` — Rao kuadratik entropisi
- `distinctDis()`, `distinctTree()` — Tur ayirt ediciligi
- `dsimTax()` — Taksonomik benzerlik/benzemezlik

**Eksikleri:**
- Deng entropisi YOK
- Ozkan pTO YOK
- Dilim proseduru YOK

**taxdiv farki:** adiv taksonomik benzemezlik hesaplayabilir ama
Deng entropisi tabanli olculer sunmaz.

---

### 5. BiodiversityR (Topluluk Ekolojisi GUI)

**Ne yapar:** vegan'in uzerine GUI (R-Commander) ekler. Cesitlilik
indekslerini hesaplar, grafik olusturur.

**Fonksiyonlari:**
- `diversityvariables()` — Shannon, Simpson, Fisher alfa, evenness
- `diversityresult()` — Alt kume bazli cesitlilik
- vegan fonksiyonlarini sarar

**Eksikleri:**
- Kendi taksonomik cesitlilik hesaplamasi YOK
- vegan'in taxondive'ina bagimli
- Deng entropisi YOK

**taxdiv farki:** BiodiversityR bir arayuz paketidir. taxdiv ise
kendi matematiksel implementasyonlarini sunar.

---

### 6. iNEXT (Rarefaction ve Ekstrapolasyon)

**Ne yapar:** Hill sayilari icin rarefaction ve ekstrapolasyon.
Orneklem buyuklugu ve kapsama tabanli egrilr.

**Fonksiyonlari:**
- `iNEXT()` — Interpolasyon/ekstrapolasyon
- `ggiNEXT()` — ggplot2 grafikleri
- `estimateD()` — Belirli orneklem boyutunda cesitlilik tahmini

**Eksikleri:**
- Sadece Hill sayilari (q=0, q=1, q=2) icin rarefaction
- Taksonomik indeksler icin rarefaction YOK
- Deng entropisi YOK

**taxdiv farki:** taxdiv'in `rarefaction_taxonomic()` fonksiyonu
8 farkli indeks icin (Shannon, Simpson, uTO, TO, uTO+, TO+, AvTD,
tur zenginligi) rarefaction hesaplayabilir. iNEXT sadece Hill
sayilari icin calisir.

---

### 7. taxize (Taksonomik Veri Erisimi)

**Ne yapar:** GBIF, ITIS, NCBI gibi veritabanlarindan taksonomik
siniflandirma bilgisi cekmek icin kullanilir. Cesitlilik hesaplamaz.

**Fonksiyonlari:**
- `classification()` — Taksonomik hiyerarsi cek
- `tax_name()` — Belirli seviyede isim al
- `synonyms()` — Sinonim listesi

**taxdiv ile iliski:** Rakip degil, tamamlayici. taxize ile taksonomik
agaci otomatik olusturup taxdiv'e girdi olarak verebilirsin.
(Hatirlatma listesinde taxize entegrasyonu var.)

---

### 8. picante (Filogenetik Cesitlilik)

**Ne yapar:** Filogenetik agac tabanli cesitlilik. Faith PD, MPD, MNTD.

**Eksikleri:** Taksonomik hiyerarsi degil filogenetik agac kullanir.
Deng entropisi YOK.

---

### 9. metacoder (Taksonomik Gorsellestirme)

**Ne yapar:** "Heat tree" gorsellestirmeleri. Metagenomik veriler icin.

**Eksikleri:** Cesitlilik indeksi hesaplamaz. Gorsellestirme odakli.

---

## Ozet Karsilastirma Tablosu

| Ozellik | taxdiv | vegan | hillR | entropart | adiv | iNEXT |
|---------|--------|-------|-------|-----------|------|-------|
| Shannon/Simpson | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ |
| Clarke & Warwick (Delta, AvTD) | ✅ | ✅ | ❌ | ❌ | ❌ | ❌ |
| Deng Entropisi | ✅ | ❌ | ❌ | ❌ | ❌ | ❌ |
| Ozkan pTO (uTO, TO, uTO+, TO+) | ✅ | ❌ | ❌ | ❌ | ❌ | ❌ |
| Dilim proseduru (Islem 1/2/3) | ✅ | ❌ | ❌ | ❌ | ❌ | ❌ |
| Taksonomik rarefaction | ✅ | ❌ | ❌ | ❌ | ❌ | ❌* |
| Taksonomik mesafe matrisi | ✅ | ✅ | ❌ | ❌ | ✅ | ❌ |
| Hill sayilari | ❌ | ❌ | ✅ | ✅ | ✅ | ✅ |
| Filogenetik cesitlilik | ❌ | ❌ | ✅ | ✅ | ✅ | ❌ |
| Bias-corrected tahminleyiciler | ❌ | ❌ | ❌ | ✅ | ❌ | ❌ |
| ggplot2 gorsellestirme | ✅ | ❌ | ❌ | ✅ | ❌ | ✅ |

*iNEXT sadece Hill sayilari icin rarefaction yapar, taksonomik indeksler icin degil.

---

## taxdiv'in Benzersiz Katkilari

### 1. Deng Entropisi Implementasyonu
CRAN'daki hicbir R paketi Deng entropisini implement etmez. taxdiv,
Dempster-Shafer kanit teorisine dayali bu entropi olcusunu R ekosistemine
tasiyan **ilk ve tek** pakettir.

### 2. Ozkan pTO Pipeline (Islem 1/2/3)
Ozkan (2018) makalesindeki tam analiz hatti (deterministik hesaplama +
stokastik yeniden ornekleme + duyarlilik analizi) sadece
"Macrotakdivozkan" adli bagimsiz bir yazilimda mevcuttur
(http://www.kantitatifekoloji.net/takdivozkan). taxdiv bunu R ekosistemine
tasir ve reproducible research icin uygun hale getirir.

### 3. Taksonomik Rarefaction
Mevcut rarefaction paketleri (iNEXT) sadece Hill sayilari icin calisir.
taxdiv, taksonomik cesitlilik indeksleri (pTO, AvTD vb.) icin de
rarefaction hesaplayabilen tek pakettir.

### 4. Entegre Gorsellestirme
vegan temel R grafikleri kullanir. taxdiv ggplot2 tabanli 7 farkli
gorsellestirme fonksiyonu sunar (dendogram, heatmap, bubble, radar,
iterasyon, karsilastirma, rarefaction).

### 5. Hem Klasik Hem Modern
taxdiv tek pakette Shannon/Simpson + Clarke & Warwick + Deng/Ozkan
indekslerini bir arada sunar. Kullanici ayri paketler yuklemek zorunda
kalmaz.

---

## Sonuc

taxdiv, R ekosisteminde **Deng entropisi tabanli taksonomik cesitlilik**
alaninda doldurulan bir bosluktur. Rakip paketlerin hicbiri bu olculeri
sunmaz. En yakin rakip vegan'dir (Clarke & Warwick indeksleri), ancak
vegan Deng entropisini ve Ozkan pTO'yu implement etmez.

Paketin akademik degeri: Ozkan (2018) makalesinin hesaplamalarini
**tekrarlanabilir** (reproducible) ve **acik kaynak** bir R paketi
olarak sunmak, taksonomik cesitlilik arastirmalarinda standart bir
arac haline gelmesini saglar.
