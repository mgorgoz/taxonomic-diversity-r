# taxdiv Rarefaction Yöntemi İnceleme Raporu

## 1. Amaç

Bu rapor, `taxdiv` paketindeki `rarefaction_taxonomic()` fonksiyonunun hangi
akademik rarefaction yaklaşımına karşılık geldiğini, kullanıcının paylaştığı
Chao et al. (2014) Hill sayıları çerçevesiyle ne kadar uyumlu olduğunu ve
pakette yapılması gereken olası düzeltme/güncellemeleri özetler.

İncelenen ana dosyalar:

- `R/rarefaction_taxonomic.R`
- `R/ozkan_pto.R`
- `R/deng_entropy.R`
- `R/classical_indices.R`
- `tests/testthat/test-rarefaction.R`
- `man/rarefaction_taxonomic.Rd`
- `inst/examples/README.md`

## 2. Kısa Sonuç

Mevcut `rarefaction_taxonomic()` fonksiyonu, Chao et al. (2014) tarafından
sunulan Hill sayıları temelli rarefaction/extrapolation yöntemini uygulamıyor.

Fonksiyonun mevcut hali en doğru şu şekilde tanımlanabilir:

> Monte Carlo ile uygulanan, birey temelli, örneklem büyüklüğüne dayalı
> rarefaction. Yöntem klasik Sanders (1968), Hurlbert (1971) ve Gotelli &
> Colwell (2001) rarefaction mantığına dayanır. Tür zenginliği için klasik
> rarefaction'a yakındır; Shannon, Simpson, AvTD ve Özkan pTO bileşenleri için
> ise aynı alt örnekleme motorunun pratik bir genellemesidir.

Bu nedenle pakette mevcut yöntemin korunması mantıklıdır; ancak dokümantasyonda
metodolojik iddia daha net yazılmalıdır. Chao-Hill yaklaşımı isteniyorsa ayrı
bir fonksiyon veya ayrı bir mod olarak eklenmelidir.

## 3. Mevcut Fonksiyon Ne Yapıyor?

`rarefaction_taxonomic()` şu işlemleri yapar:

1. Kullanıcıdan isimli bir bolluk vektörü alır.
2. Sıfır bolluklu türleri atar.
3. Bollukları tam sayıya yuvarlar.
4. Bolluk vektörünü birey düzeyine genişletir.
5. Farklı örneklem büyüklükleri belirler.
6. Her örneklem büyüklüğünde `n_boot` kez rastgele alt örneklem alır.
7. Alt örneklem alma işlemini `replace = FALSE` ile yapar.
8. Her alt örneklem için seçilen çeşitlilik indeksini yeniden hesaplar.
9. Ortalama, alt güven sınırı, üst güven sınırı ve standart sapma döndürür.

Basit örnek:

```r
community <- c(sp1 = 5, sp2 = 3, sp3 = 1, sp4 = 2)
```

Fonksiyon bunu kavramsal olarak şu birey listesine çevirir:

```text
sp1 sp1 sp1 sp1 sp1 sp2 sp2 sp2 sp3 sp4 sp4
```

Sonra örneğin 4 bireylik alt örnekler çeker:

```text
sp1 sp2 sp2 sp4
sp1 sp1 sp3 sp4
sp2 sp2 sp3 sp4
...
```

Her çekilişte tür bolluklarını yeniden sayar ve seçilen indeksi hesaplar.
Bu işlem birçok kez tekrarlandığı için tek bir rastgele çekilişe bağlı kalınmaz;
ortalama eğri elde edilir.

## 4. Akademik Olarak Hangi Yönteme Karşılık Geliyor?

### 4.1. Sanders (1968)

Sanders, rarefaction yaklaşımını ekolojik çeşitlilik karşılaştırmalarında
kullanılan erken kaynaklardan biridir. Temel fikir, farklı büyüklükteki
örnekleri aynı örneklem büyüklüğüne indirerek karşılaştırmaktır.

Bu tarihsel arka plan paketteki yaklaşım ile uyumludur: farklı örneklem
büyüklüklerinde beklenen çeşitlilik değeri hesaplanır.

### 4.2. Hurlbert (1971)

Hurlbert, rarefaction'ı tür zenginliği için daha formel hale getirir. Ana fikir:
N bireyden oluşan bir koleksiyondan n birey rastgele ve yerine koymadan
seçilirse, beklenen tür sayısı nedir?

`taxdiv` bunu tür zenginliği için analitik formülle değil, Monte Carlo
alt örnekleme ile yaklaşık hesaplıyor. Yani her n için beklenen değeri doğrudan
formülden almak yerine birçok rastgele alt örneklem çekip ortalamasını alıyor.

### 4.3. Gotelli & Colwell (2001)

Gotelli & Colwell (2001), birey temelli rarefaction ile örnek temelli
rarefaction ayrımını, örneklem büyüklüğü standardizasyonunu ve rarefaction
eğrilerinin doğru yorumlanmasını açıklar. Paket en çok bu çerçeveye yakındır.

`rarefaction_taxonomic()` birey temelli rarefaction yapar:

- Veri birey düzeyine genişletilir.
- Alt örneklem bireyler arasından seçilir.
- Her alt örneklemde indeks yeniden hesaplanır.
- Eğri örneklem büyüklüğü ekseninde çizilir.

## 5. Chao et al. (2014) Hill Sayıları Yaklaşımıyla Farkı

Kullanıcının paylaştığı metinde Chao et al. (2014) Hill sayıları çerçevesinin
Özkan pTO indekslerine teorik olarak uygulanabileceği belirtilmişti. Bu önemli
bir fikir, fakat mevcut paket bu yöntemi uygulamıyor.

Chao et al. (2014) yaklaşımının temel özellikleri:

- Hill sayıları kullanır.
- `q = 0`, `q = 1`, `q = 2` gibi çeşitlilik mertebeleri vardır.
- `q = 0` tür zenginliğine karşılık gelir.
- `q = 1` Shannon entropisinin üstel dönüşümüdür: `exp(H)`.
- `q = 2` genellikle inverse Simpson çeşitliliğine karşılık gelir.
- Sample-size rarefaction yanında coverage-based rarefaction da sunar.
- Extrapolation, yani gözlenen örneklem büyüklüğünün ötesine tahmin yapabilir.
- Gözlenmeyen türler ve örneklem tamamlanmamışlığı için tahmin edici yapı içerir.

Mevcut `taxdiv` fonksiyonunda bunlar yoktur:

- `q` parametresi yok.
- Shannon `exp(H)` ile etkin tür sayısına çevrilmiyor.
- Simpson Hill sayısı ölçeğinde raporlanmıyor.
- Sample coverage hesabı yok.
- Coverage standardizasyonu yok.
- Extrapolation yok.
- Gözlenmeyen tür tahmini yok.
- pTO için matematiksel bir Hill-pTO dönüşümü yok.

Bu nedenle mevcut fonksiyon için "Chao et al. (2014) Hill-number
rarefaction/extrapolation uygulanmıştır" demek doğru olmaz.

## 6. Özkan pTO İçin Mevcut Yaklaşım Ne Anlama Geliyor?

Paket pTO bileşenleri için özel bir analitik rarefaction formülü kullanmıyor.
Onun yerine şunu yapıyor:

1. Alt örneklem bireyleri seçiliyor.
2. Alt örneklemden yeni bir topluluk vektörü oluşturuluyor.
3. Bu vektör `ozkan_pto()` fonksiyonuna veriliyor.
4. `uTO`, `TO`, `uTO_plus` veya `TO_plus` değeri alınıyor.
5. Bu işlem birçok kez tekrarlanıp ortalama eğri çıkarılıyor.

Bu, pTO için pratik ve anlaşılır bir yaklaşımdır. Çünkü pTO algoritmasının kendi
iç yapısı değiştirilmez; sadece daha küçük örneklemlerde pTO'nun nasıl
davrandığı görülür.

Ancak bu, literatürde özel olarak türetilmiş bir "pTO rarefaction formülü"
değildir. Bu nedenle dokümantasyonda "practical subsampling extension" veya
"Monte Carlo subsampling-based rarefaction" gibi daha dikkatli ifadeler
kullanılmalıdır.

## 7. Mevcut Yöntemin Güçlü Yanları

- Uygulaması basit ve anlaşılırdır.
- pTO dahil farklı indekslere aynı mantıkla uygulanabilir.
- Kullanıcıya örneklem büyüklüğü arttıkça indeksin nasıl değiştiğini gösterir.
- Tam örneklem büyüklüğünde değer, gerçek indeks değerine yaklaşır.
- Taksonomik indeksler için pratik olarak çalışır.
- Paket bağımlılıklarını artırmaz.
- Chao-Hill gibi daha teorik bir dönüşüm gerektirmez.

## 8. Mevcut Yöntemin Sınırlılıkları

### 8.1. "Bootstrap" terimi tam doğru değil

Kodda ve dokümantasyonda yöntem "bootstrap" olarak adlandırılıyor. Ancak alt
örnekleme `replace = FALSE` ile yapılıyor. Klasik bootstrap genellikle yerine
koyarak örnekleme yapar (`replace = TRUE`).

Bu nedenle "bootstrap rarefaction" yerine şu ifadeler daha doğru olur:

- "Monte Carlo subsampling"
- "repeated random subsampling without replacement"
- "individual-based sample-size rarefaction"

Güven aralıkları için "bootstrap CI" yerine "simulation interval" veya
"Monte Carlo interval" daha temkinli bir ifade olur.

### 8.2. Bollukların yuvarlanması riskli

Fonksiyon bollukları `round()` ile tam sayıya çeviriyor. Bu birey sayısı
verileri için sorun değildir. Ancak paket dokümantasyonu Özkan pTO için
Westhoff-Maarel gibi kaplama-bolluk ölçeklerini de kabul ettiğini belirtiyor.
Bu tür ordinal veya yarı-kantitatif değerleri birey sayısı gibi yuvarlayıp
rarefaction yapmak ekolojik olarak tartışmalıdır.

Öneri:

- Rarefaction için bollukların non-negative integer count olması gerektiği
  açıkça belirtilmeli.
- Non-integer bolluk varsa uyarı verilmeli.
- `round()` sessizce yapılmamalı veya kullanıcıya açıkça bildirilmeli.

### 8.3. Örneklem büyüklüğü başlangıcı tartışmalı

Kodda minimum örneklem büyüklüğü:

```r
min_n <- max(2, min(community[community > 0]))
```

Bu, eğrinin en küçük pozitif bolluktan başlamasına neden olur. Klasik
rarefaction eğrileri genellikle 1 veya 2 bireyden başlayabilir. Bu seçim
dokümantasyonda açıklanmamış.

Öneri:

- `min_size` parametresi eklenebilir.
- Varsayılan `min_size = 1` veya indekslerin anlamlılığı için `2` yapılabilir.
- Mevcut davranış korunacaksa gerekçesi yazılmalıdır.

### 8.4. Chao-Hill yöntemiyle karıştırılabilir

Paket dokümantasyonunda rarefaction genel olarak anlatılıyor, ancak Chao-Hill
çerçevesinden açıkça ayrılmıyor. Kullanıcılar "rarefaction" deyince iNEXT veya
Chao et al. (2014) tarzı Hill number rarefaction/extrapolation bekleyebilir.

Öneri:

- Dokümantasyona şu tür bir not eklenmeli:

```text
This function performs empirical Monte Carlo sample-size rarefaction by
subsampling observed individuals without replacement. It does not implement
Hill-number rarefaction/extrapolation or coverage-based standardization as in
Chao et al. (2014).
```

## 9. Önerilen Paket Güncellemeleri

### 9.1. Kısa vadeli düzeltmeler

Bu düzeltmeler mevcut yöntemi değiştirmeden paketi metodolojik olarak daha
doğru hale getirir.

1. `rarefaction_taxonomic()` dokümantasyonunda yöntemin adı netleştirilmeli.
2. "bootstrap" ifadesi "Monte Carlo subsampling" olarak düzeltilmeli.
3. `replace = FALSE` olduğu özellikle belirtilmeli.
4. Chao et al. (2014) Hill-number rarefaction olmadığı açıkça yazılmalı.
5. Non-integer abundance durumunda uyarı eklenmeli.
6. `round()` davranışı ya kaldırılmalı ya da kullanıcıya bildirilmeli.
7. `min_n` başlangıç mantığı dokümante edilmeli veya parametre haline getirilmeli.
8. Testlere non-integer abundance uyarısı ve `min_size` davranışı eklenmeli.

### 9.2. Orta vadeli geliştirmeler

1. `method` parametresi eklenebilir:

```r
method = c("monte_carlo", "analytical")
```

İlk aşamada analitik yöntem sadece `index = "species"` için desteklenebilir.
Bu, Hurlbert tarzı beklenen tür zenginliği formülünü doğrudan verir.

2. `replace` parametresi eklenebilir, fakat varsayılan rarefaction için
`replace = FALSE` kalmalıdır.

3. Çıktı özniteliklerine yöntem bilgisi eklenebilir:

```r
attr(results, "method") <- "monte_carlo_without_replacement"
```

### 9.3. Uzun vadeli geliştirmeler

Chao et al. (2014) çerçevesi ayrı bir fonksiyon olarak eklenmeli:

```r
hill_rarefaction()
```

veya

```r
rarefaction_hill()
```

Bu fonksiyon şu özellikleri içermeli:

- `q = c(0, 1, 2)` desteği
- sample-size rarefaction
- coverage-based rarefaction
- extrapolation
- effective diversity çıktısı
- uygun referanslar ve metodolojik açıklama

Özkan pTO için Hill dönüşümü yapılacaksa bu ayrı ve dikkatli tanımlanmalıdır.
Sadece `exp(pTO)` almak otomatik olarak geçerli bir Hill sayısı üretmez.
Önce pTO'nun hangi entropy ailesine, hangi q mertebesine ve hangi etkin
çeşitlilik birimine dönüştüğü matematiksel olarak tanımlanmalıdır.

## 10. Önerilen Dokümantasyon Metni

`R/rarefaction_taxonomic.R` roxygen açıklamasına şu yönde bir metin eklenebilir:

```text
This function implements empirical individual-based sample-size rarefaction by
Monte Carlo subsampling of observed individuals without replacement. For each
target sample size, it repeatedly draws a subsample, recalculates the selected
index, and summarizes the resulting distribution.

For species richness, this approximates the classical individual-based
rarefaction expectation described by Sanders (1968), Hurlbert (1971), and
Gotelli & Colwell (2001). For Shannon, Simpson, AvTD, and Ozkan pTO components,
the same subsampling framework is used as a practical extension.

The function does not implement Hill-number rarefaction/extrapolation,
coverage-based standardization, or unseen-species estimation as in Chao et al.
(2014).
```

## 11. Önerilen Test Güncellemeleri

Mevcut `tests/testthat/test-rarefaction.R` dosyası temel davranışı test ediyor
ve çalışıyor. İnceleme sırasında bu test dosyası çalıştırıldı; 24 test geçti.

Eklenmesi önerilen testler:

1. Non-integer abundance verildiğinde uyarı bekleyen test.
2. `replace = FALSE` mantığının tam örneklemde deterministik sonuç verdiğini
   kontrol eden test.
3. `method` veya `min_size` parametresi eklenirse bunların davranış testleri.
4. `index = "species"` için analitik Hurlbert sonucu ile Monte Carlo sonucunun
   yeterli tekrar sayısında yakınsadığını kontrol eden test.
5. Dokümantasyonda belirtilen Chao-Hill özelliklerinin mevcut fonksiyonda
   iddia edilmediğini garanti eden snapshot veya metin testi.

## 12. Arkadaşa Verilecek Net İş Listesi

1. `rarefaction_taxonomic()` fonksiyonunun metodolojik açıklamasını düzelt.
2. "bootstrap" kelimesini kullanım bağlamına göre "Monte Carlo subsampling"
   veya "repeated subsampling without replacement" olarak değiştir.
3. Fonksiyonun Chao et al. (2014) Hill-number rarefaction yapmadığını
   dokümantasyonda açıkça belirt.
4. Non-integer bolluklar için uyarı ekle; sessiz `round()` davranışını gözden
   geçir.
5. `min_n` başlangıç değerini parametre haline getir veya gerekçesini yaz.
6. İleride Chao-Hill yaklaşımı isteniyorsa mevcut fonksiyonu bozma; ayrı
   `hill_rarefaction()` fonksiyonu olarak geliştir.
7. pTO için Hill sayısı dönüşümü eklenecekse önce matematiksel tanımı
   netleştir; bunu mevcut pTO rarefaction'ın yerine koyma.

## 13. Kullanılacak Doğru Kısa Tanım

Paket dokümantasyonu veya makale metni için önerilen kısa tanım:

> `rarefaction_taxonomic()` uses Monte Carlo individual-based sample-size
> rarefaction. It repeatedly subsamples observed individuals without replacement
> at increasing sample sizes and recalculates the selected diversity index. For
> species richness this follows the classical rarefaction logic of Sanders
> (1968), Hurlbert (1971), and Gotelli & Colwell (2001). For Shannon, Simpson,
> AvTD, and Ozkan pTO components, the same subsampling framework is used as a
> practical extension. The function does not implement Hill-number
> rarefaction/extrapolation or coverage-based standardization.

Türkçe kısa tanım:

> `rarefaction_taxonomic()` gözlenen bireyler üzerinden yerine koymadan tekrarlı
> alt örnekleme yapan, birey temelli örneklem-büyüklüğü rarefaction yöntemidir.
> Tür zenginliği için klasik rarefaction mantığını yaklaşıklar; Shannon,
> Simpson, AvTD ve Özkan pTO bileşenleri için aynı alt örnekleme çerçevesini
> pratik olarak genişletir. Chao et al. (2014) tarzı Hill sayıları,
> coverage-based standardizasyon veya extrapolation uygulamaz.

## 14. Kaynaklar

- Sanders, H. L. (1968). Marine benthic diversity: a comparative study.
  The American Naturalist, 102, 243-282.
- Hurlbert, S. H. (1971). The nonconcept of species diversity: a critique and
  alternative parameters. Ecology, 52, 577-586.
- Gotelli, N. J. & Colwell, R. K. (2001). Quantifying biodiversity: procedures
  and pitfalls in the measurement and comparison of species richness.
  Ecology Letters, 4, 379-391.
  https://www.uvm.edu/~ngotelli/manuscriptpdfs/EcologyLetters4p379.pdf
- Chao, A. et al. (2014). Rarefaction and extrapolation with Hill numbers: a
  framework for sampling and estimation in species diversity studies.
  Ecological Monographs, 84, 45-67.
  https://harvardforest.fas.harvard.edu/sites/default/files/ellison-pubs/2014/Chao_etal_2014_EM.pdf
- Deng, Y. (2016). Deng entropy. Chaos, Solitons & Fractals, 91, 549-553.

