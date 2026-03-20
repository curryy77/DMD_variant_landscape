<p align="center">
  <a href="README_en.md">🇬🇧 English</a> |
  <a href="README_ja.md">🇯🇵 日本語</a> |
  <a href="README_fr.md">🇫🇷 Français</a> |
  <a href="README_ru.md">🇷🇺 Русский</a> 
</p>

# DMD Variant Landscape Analysis

![License](https://img.shields.io/badge/license-MIT-green)
![Python](https://img.shields.io/badge/python-3.10%2B-blue)
![Bioinformatics](https://img.shields.io/badge/field-bioinformatics-orange)
![Data](https://img.shields.io/badge/data-ClinVar%20%7C%20gnomAD%20%7C%20Ensembl-purple)

---
### About
Exploratory analysis of pathogenic and non-pathogenic variants in the **DMD gene** using annotated variant datasets. The goal is to investigate structural and functional patterns of pathogenicity across **exons**, **protein domains**, and **reading frame status**.

---

### Datasets
The analysis uses an annotated dataset of **DMD variants** (merged from **ClinVar**, **Ensembl** and **gnomAD** databases) containing genomic position, exon number, protein domain, pathogenicity classification, phenotype class (DMD or BMD), reading frame status (in-frame or out-of-frame). Variants are preprocessed and stored in:
```commandline
data/processed/DMD_variants_annotated.csv
```

---
### Analyses performed (so far)

#### 1. Exon-level distribution of variants
First, we analyzed the distribution of pathogenic and non-pathogenic variants across all **79 exons of the DMD gene**.

![Exon distribution](figures/0XX_exon_distribution.png)

Two complementary visualizations are produced: **variant counts per exon**, and **fraction of pathogenic variants per exon**.

What we can observe is a strong concentration of variants in exons 42-56, corresponding to deletion hotspot region in DMD, described in multiple studies (Monaco et al., 1988; Aartsma-Rus et al., 2006; Bladen et al., 2015).

This region overlaps with the central rod domain of dystrophin, which contains multiple spectrin-like repeats.

#### 2. Domain level pathogenicity
Variants were mapped (mapping provided by UniProt) to functional **protein domains** of dystrophin. For each domain we computed total number of variants and fraction of pathogenic variants.

![Domain distribution](figures/1XX_domain_distribution.png)

Domains with the highest pathogenic fraction include **actin-binding domain**, **SNTB1 binding region**, **interaction with syntrophin**. These domains correspond to **critical interaction interfaces** required for dystrophin's structural role in the dystrophin-glycoprotein complex. Mutations in such regions are expected to more frequently disrupt protein function and there show higher pathogenicity (Ervasti & Campbell, 1993; Koenig et al., 1988; Blake et al., 2002).

#### 3. Reading frame rule analysis
The classical reading frame rule, described in (Monaco et al., 1988) in dystrophinopathies was tested. The rule states, that the out-of-frame mutations tend to lead to severe, DMD case, while in-frame mutations are resulting in more milder, Becker's muscular dystrophy case.

![Reading frame rule](figures/2XX_reading_frame_rule.png)

Variant counts by phenotype and frame status alongside with proportion plot showing phenotype fractions was generated. 

The results follow the expected trend: out-of-frame mutations are strongly enriched for DMD, while in-frame mutations show a higher proportion of BMD. However, exceptions are observed. Such exceptions may arise due to disruption of critical functional domains, alternative splicing effects and structural instability of truncated dystrophin.

---
### Biological interpretation

The analyses highlight several important properties of the DMD mutational landscape: **mutation hotspots correspond to structurally repetitive regions** of the gene, **functional binding domains** show **elevated pathogenicity rates**, the **Monaco's reading frame rule** is **largely supported**, but not absolute.

---
### TODO

Planned extensions of the analysis include: **statistical enrichment tests**, **odds-ratio analysis** for domain pathogenicity, exon-level mutation density **modeling**, **integration with structural data** of dystrophin domains, comparison with **population variant frequencies**.

---
### Data Sources

The analysis integrates variant information from several public databases:

- **ClinVar** - clinical variant annotations
- **gnomAD** - population allele frequencies
- **Ensembl** - genomic annotation and transcript structure
- **GEO (GSE38417)** - gene expression dataset

We thank the maintainers of the upper-mentioned databases for providing open access to genomic and clinical data.

---
### License

This project is released under the MIT License.

See the LICENSE file for details.

---

### Project Status: 🟨 🟨 🟨 
This repository contains an exploratory analysis pipeline for the DMD mutational landscape. 

**The project is currently under active development!** This README.md file is temporary and will be changed after the project is complete.

---

© 2026, Mikhail Kolesnikov \
Moscow, Higher School of Economics, Faculty of Computer Science, BSc 

MIT License\
GitHub: https://github.com/curryy77