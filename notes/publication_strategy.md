# Yayin Stratejisi / Publication Strategy

## Karsilastirma: rOpenSci vs JOSS vs CRAN

| Ozellik | rOpenSci | JOSS | CRAN |
|---------|----------|------|------|
| **Nedir** | Topluluk kod incelemesi | Akademik yazilim dergisi | Paket dagitim deposu |
| **Cikti** | rOpenSci ekosisteminde yer | DOI'li atif yapilabilir makale | CRAN'da paket |
| **Hakemli mi** | Evet (1-3 hakem, kod + dokumantasyon) | Evet (2 hakem, kod + makale) | Hayir (otomatik kontrol) |
| **Sure** | ~2-3 ay | ~1-3 ay | Gunler - haftalar |
| **Maliyet** | Ucretsiz | Ucretsiz | Ucretsiz |
| **Makale gerekli mi** | Hayir | Evet (750-1750 kelime, Markdown) | Hayir |
| **Akademik kredi** | Orta (topluluk taninirlik) | Yuksek (atif yapilabilir yayin) | Dusuk (sadece dagitim) |
| **taxdiv icin uygun mu** | Belirsiz - kapsam disi olabilir | Cok uygun | Her zaman uygun |

---

## JOSS (Journal of Open Source Software)

### Neden taxdiv icin en uygun?

1. **Arastirma uygulamasi acik:** Taksonomik cesitlilik analizi = net bilimsel kullanim
2. **Yeterince kapsamli:** 21+ fonksiyon, vignette'ler, testler = trivial degil
3. **Ucretsiz:** Hicbir odeme yok (diamond open access)
4. **Kisa makale:** 750-1750 kelime, Markdown formatinda
5. **DOI alinir:** Atif yapilabilir akademik yayin

### JOSS Makale Formati

Repoya `paper.md` ve `paper.bib` dosyalari eklenir:

```
paper.md icerigi:
- Summary (paket ne yapar)
- Statement of Need (neden gerekli)
- State of the Field (rakip paketlerle karsilastirma)
- Software Design (mimari)
- Research Impact
- AI Usage Disclosure
- References
```

### JOSS Gereksinimleri

- Acik kaynak (OSI lisansi) -> MIT var ✅
- Acik depo -> GitHub var ✅
- Test altyapisi -> testthat var ✅
- CI/CD -> GitHub Actions var ✅
- Dokumantasyon -> Vignette'ler var ✅
- Onemli arastirma uygulamasi -> Evet ✅

---

## rOpenSci

### Sorun: Kapsam uygunlugu belirsiz

rOpenSci'nin ana odagi **veri yasam dongusu** (veri erisimi, temizleme,
dogrulama, depolama). Taksonomi kategorisindeki paketleri (taxize gibi)
veritabanlarina erisim icin, cesitlilik hesaplama icin degil.

Istatistiksel yazilim incelemesi izinde (EDA/summary statistics) kabul
edilebilir ama garanti degil.

### Oneri: On-basvuru sorgusu yap

rOpenSci'ye tam basvurudan once `ropensci/software-review` deposunda
pre-submission inquiry ac. "taxdiv kapsaminizda mi?" diye sor.

---

## CRAN

### Her durumda yapilmali

CRAN basvurusu JOSS veya rOpenSci'den bagimsiz. `install.packages("taxdiv")`
ile kurulabilir hale getirmek icin CRAN'a girmek sart.

### CRAN Gereksinimleri

- `R CMD check --as-cran` hatasiz gecmeli (0 error, 0 warning)
- Tek sorumlu bakim yapici (maintainer)
- Lisans acik olmali
- Kaynak tarball < 10MB
- Platformlar arasi uyumluluk

---

## Onerilen Eylem Plani

### Adim 1: CRAN Basvurusu (Oncelikli)
- `devtools::check()` -> 0/0/0
- `rhub::check_for_cran()` ile ek kontrol
- https://CRAN.R-project.org/submit.html uzerinden gonder

### Adim 2: JOSS Basvurusu (Esanlamli veya hemen sonra)
- `paper.md` ve `paper.bib` hazirla
- https://joss.theoj.org uzerinden gonder
- ~1-3 ay icinde hakemli yayin

### Adim 3: rOpenSci (Opsiyonel)
- On-basvuru sorgusu ac
- Kabul edilirse JOSS'a hizli gecis (fast-track)
- Reddedilirse JOSS'a dogrudan devam et

---

## Sonuc

**JOSS + CRAN** = taxdiv icin en iyi strateji.
- CRAN: Kolay erisim ve dagitim
- JOSS: Akademik atif ve hakemli yayin
- Toplam maliyet: 0 TL
- Toplam sure: ~2-4 ay
