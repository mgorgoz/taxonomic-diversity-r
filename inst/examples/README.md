# Example Outputs / Ornek Ciktilar

This folder contains example visualizations produced by the `taxdiv` package
using a hypothetical Mediterranean forest community dataset with 10 tree
species and 7 taxonomic levels (Species, Genus, Family, Order, Class, Phylum,
Kingdom).

Bu klasor, `taxdiv` paketinin 10 agac turu ve 7 taksonomik seviyeli (Tur,
Cins, Familya, Takim, Sinif, Bolum, Alem) hayali bir Akdeniz ormani toplulugu
verisiyle urettigi ornek gorselleri icerir.

---

## Example Communities / Ornek Topluluklar

Two community scenarios are used throughout the examples to illustrate how
diversity indices respond to different abundance distributions:

Orneklerde iki farkli topluluk senaryosu kullanilmistir. Amac, cesitlilik
indekslerinin farkli bolluk dagilimlarinda nasil degistigini gostermektir:

### "Diverse" (Cesitli / Dengeli Dagilim)

Species abundances are relatively even across all 10 species. No single
species strongly dominates the community.

Tur bolluklari 10 tur arasinda gorece dengeli dagilmistir. Hicbir tur
topluluga belirgin sekilde hakim degildir.

```
Pinus_brutia = 30, Quercus_coccifera = 25, Fagus_orientalis = 20,
Quercus_infectoria = 18, Cedrus_libani = 15, Pinus_nigra = 12,
Carpinus_betulus = 10, Juniperus_excelsa = 8, Abies_cilicica = 7,
Juniperus_oxycedrus = 5
```

### "Dominant" (Baskin Turlu / Esitsiz Dagilim)

One species (Quercus coccifera) has disproportionately high abundance (80),
while all other species have very low abundances (1-5). This represents a
community with low evenness.

Bir tur (Quercus coccifera) orantisiz olarak yuksek bolluga sahiptir (80),
diger turlerin hepsi cok dusuk bolluktadir (1-5). Bu, dusuk esitlige
(evenness) sahip bir toplulugu temsil eder.

```
Quercus_coccifera = 80, Quercus_infectoria = 5, Pinus_brutia = 3,
Cedrus_libani = 3, Pinus_nigra = 2, Juniperus_excelsa = 2,
Fagus_orientalis = 2, Juniperus_oxycedrus = 1, Abies_cilicica = 1,
Carpinus_betulus = 1
```

### Why these names? / Neden bu isimler?

"Diverse" and "Dominant" are **descriptive labels**, not formal ecological
classification terms. They describe the **evenness** (esitlik) of the
abundance distribution:

"Diverse" ve "Dominant" resmi ekolojik siniflandirma terimleri **degildir**.
Bolluk dagiliminin **esitligini** (evenness) betimleyen etiketlerdir:

- **High evenness** = abundances spread relatively equally among species
  ("Diverse" scenario)
- **Low evenness** = one or few species dominate in abundance
  ("Dominant" scenario)

The concept of evenness is formalized by Pielou's J index:

  J = H' / ln(S)

where H' is Shannon entropy and S is the number of species. J ranges from
0 (one species completely dominates) to 1 (all species equally abundant).

Esitlik kavrami Pielou'nun J indeksi ile formalize edilmistir:

  J = H' / ln(S)

H' Shannon entropisi, S tur sayisidir. J degeri 0 (bir tur tamamen baskin)
ile 1 (tum turler esit bollukta) arasinda degisir.

---

## File Descriptions / Dosya Aciklamalari

| File | Function | Description |
|------|----------|-------------|
| `dendrogram_family.png` | `plot_taxonomic_tree()` | Dendrogram colored by Family, with abundance values |
| `dendrogram_class.png` | `plot_taxonomic_tree()` | Dendrogram colored by Class (Magnoliopsida vs Pinopsida) |
| `dendrogram_no_abundance.png` | `plot_taxonomic_tree()` | Dendrogram without abundance, colored by Order |
| `compare_indices_barplot.png` | `compare_indices()` | Bar plot comparing Diverse vs Dominant across all 10 indices |

---

## Key Observation / Onemli Gozlem

The bar plot comparison reveals an important ecological distinction:

Bar plot karsilastirmasi onemli bir ekolojik ayrimi ortaya koyar:

**Abundance-dependent indices** (bolluga bagimli indeksler):
Shannon, Simpson, Delta, uTO, TO — these change between Diverse and
Dominant scenarios because they incorporate species abundance weights.

**Abundance-independent indices** (bolluktan bagimsiz indeksler):
AvTD, VarTD, uTO+, TO+ — these remain identical in both scenarios
because they depend only on the species list (presence/absence), not
on how many individuals each species has.

This is a fundamental property of these indices documented by
Clarke & Warwick (1998) and reflected in Ozkan's (2018) pTO framework.

## References / Kaynaklar

- Clarke, K.R. & Warwick, R.M. (1998). A taxonomic distinctness index and
  its statistical properties. *Journal of Applied Ecology*, 35, 523-531.
- Ozkan, K. (2018). A new taxonomic diversity index. *Turkish Journal of
  Forestry*, 19(4), 340-346.
- Pielou, E.C. (1966). The measurement of diversity in different types of
  biological collections. *Journal of Theoretical Biology*, 13, 131-144.
- Deng, Y. (2016). Deng entropy. *Chaos, Solitons & Fractals*, 91, 549-553.
