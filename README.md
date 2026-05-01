<p align="center">
  <a href="docs/README_en.md">🇬🇧 English</a> |
  <a href="docs/README_ja.md">🇯🇵 日本語</a> |
  <a href="docs/README_fr.md">🇫🇷 Français</a> |
  <a href="docs/README_ru.md">🇷🇺 Русский</a> 
</p>

# Analysis of the DMD Gene Variant Landscape: an Exploratory Study Based on Merged Data from ClinVar, Ensembl, and gnomAD Sources

![License](https://img.shields.io/badge/license-MIT-green)
![Python](https://img.shields.io/badge/python-3.10%2B-blue)
![Bioinformatics](https://img.shields.io/badge/field-bioinformatics-orange)
![Data](https://img.shields.io/badge/data-ClinVar%20%7C%20gnomAD%20%7C%20Ensembl-purple)

---
### Abstract

This repository contains a research analysis of pathogenic and non-pathogenic variants in the human **DMD gene** using merged data from bioinformatics sources. The goal of the study is to examine structural and functional pathogenicity patterns within **exons**, **protein domains**, **reading-frame status**, etc.

Cleaning and merging of the original datasets were performed; based on the resulting dataset, an **exploratory analysis** of variants was conducted; two **prediction models** of pathogenicity for DMD gene variants were also built based on **logistic regression** (LogisticRegression) and **random forest** (RandomForest). Additionally, a **gene expression analysis** was performed for healthy patients and patients with **Duchenne or Becker muscular dystrophy**, in order to identify the effect of disease on **DMD gene** expression and **overall gene expression** in the organism.

---

### Contents

- [Project task](#project-task)
- [Main points and project pipeline](#main-points-and-project-pipeline)
- [Datasets](#datasets)
- [Methods](#methods)
- [Obtained model metrics](#obtained-model-metrics)
- [EDA: consistency with the literature](#eda-consistency-with-the-literature)
- [Research-gap findings (??)](#research-gap-findings-)
- [Project structure](#project-structure)
- [Reproducibility](#reproducibility)
- [Next steps](#next-steps)
- [License](#license)
- [Project status](#project-status----)

---

### Project task

The **DMD gene** is the largest known gene in the human genome, located on chromosome X, encoding the **dystrophin** protein. This protein links the **structural framework** of muscle cells with the **surrounding matrix**. Pathogenic variants in DMD lead to **disruption of dystrophin synthesis**, resulting in a group of **severe inherited diseases**, including **Duchenne muscular dystrophy** (DMD), a severe form leading to rapid progression of muscle weakness, loss of ambulation and cardiomyopathy, and **Becker muscular dystrophy** (BMD), a milder form characterized by slowly progressive muscle weakness.

The project task is to **study various features** and patterns of variants in the gene that define their **pathogenicity** and **non-pathogenicity**, including the subsequent possibility of **predicting** based on them the **pathogenicity of later discovered variants**.

---

### Main points and project pipeline

An image showing the **pipeline** of this project:

![Project Pipeline2](assets/pipeline/project_pipeline_ru.png)

At the **data preparation** stage, initially **three tables** on DMD gene variants were collected from **ClinVar**, **Ensembl**, and **gnomAD**. Next, keys were **unified** (for example, ``var_id``, ``rsid``), coordinate and categorical field formats were **brought** to a single standard, **technical duplicates** and empty categories were **removed**, and numerical features were **normalized** with explicit missing-value handling. After that, step-by-step **"merging"** was applied (left-join from the main variant table), at each step **controlling** row loss, growth of missingness, and field consistency. Code for all of the above is located in ``src/data_preparation.py``.

Next, variants are **annotated**: in ``src/annotate_variants.py`` the script loads the **variant table** and two **reference tables** (boundaries of dystrophin protein domains and DMD exon coordinates), then assigns a domain by **amino-acid position** and an exon by **genomic coordinates** for each variant. After that, derived features are formed, such as **mutation type** and **reading-frame status**: mutation type is aggregated from consequence/variant-type columns with a separate rule for large events when a variant occurs at the level of **several nucleotides** and its length is greater than 50, and reading-frame status is inferred from the **mutation class**. The output is ``DMD_variants_annotated.csv``, the dataset used for **further work**.

Exploratory analysis starts with ``notebooks/00_data_preview.ipynb`` and ``notebooks/01_annotation_check.ipynb``, and then blocks in ``notebooks/02_EDA`` are processed **systematically**. At initial stages, **uniqueness** of ``var_id`` and ``rsid``, duplicates, missing share, consistency of ``chr`` and ``pos``, coverage by domains, exons, and reading-frame status were checked. As a result, it became clear that data are overall suitable for statistics, however missingness is distributed unevenly, because it is more often concentrated in more **complex** variants and separate mutation classes, therefore missingness could not be interpreted as random noise.

Further, tests provided the **main picture of the spectrum** by phenotype, mutation type, reading-frame status, domain, ``consequence``, interval length, ``REVEL/meta_lr`` metrics, allele frequency, and also by exons and amino-acid positions. In most cross-tab and rank-based comparisons, the structure turned out to be **non-random**: clinical groups really differ by mutation types and molecular context, pathogenic variants are more often associated with more “radical” ``consequence``, and population features behave as expected toward lower frequency in clinically significant variants.

Checking of ``⭐`` hypotheses (results from scientific literature) overall **confirmed known literature supports**: reading-frame rule (relationship between reading-frame status and DMD/BMD phenotype), exon hotspot architecture (including distal zone 45-55), **clinical relevance** of isoform and distal regions (dp140, dp71), and **therapeutic focus** of skip regions (including skip51). Domain-structural observations and functional ``score`` metrics REVEL and MetaLR as pathogenicity markers are also consistent in a similar way.

``📖`` blocks showed rare ``consequence`` enrichments, entropy profiles, **analysis of exceptions** to the reading-frame status rule, and domain-positional anomalies. The main result is **identification of concrete zones** where annotation diverges or remains **under-described**. These places require **external validation** on an independent set and deeper revision, and may be either a real new signal or an artifact of this source.

In building prediction models, the task of **binary classification** of DMD gene variants into ``pathogenic`` and ``benign`` was solved using annotated features from the common pipeline (``exon/domain/mutation/frame`` and other fields). Data preparation was done with emphasis on filtering noisy categories, duplicate control, feature engineering, and also robustness checking with cross-validation. Two interpretable baseline models were trained — **LogisticRegression** and **RandomForest**, saved as ``.pkl``. Evaluation was built on a set of metrics such as ``accuracy``, ``AUROC``, ``confusion matrix``.

Quality results turned out **strong**: a high level of accuracy was reached (for LogisticRegression: accuracy ``~0.982``, AUROC ``~0.995``). High AUROC was considered together with threshold-dependent metrics and leakage risks, which makes conclusions reproducible. Nevertheless, **external validation is required** on an independent dataset.

![auroc_1](assets/auroc.png)

In the **expression block**, we used the **Gene Expression Omnibus** dataset ``GSE38417``, where the **expression matrix** and **phenotypes** were extracted, and **DMD** and **control** groups were formed. After that, a **differential analysis** was performed to see how ``DMD`` and related markers behave in patients relative to control.

The results are as follows: differential signals highlight **muscle regeneration genes** (including ``MYH3/MYH8``), which agrees with chronic damage and **regenerative response** in Duchenne muscular dystrophy.

![expression](assets/expression.png)

---

### Datasets

The project uses an integrated dataset of **DMD gene** variants, collected from open sources **ClinVar** (clinical interpretation), **Ensembl** (genomic-transcript annotation), and **gnomAD** (population frequencies). The final working table is stored in:
```commandline
data/processed/DMD_variants_annotated.csv
```

The table contains key fields for **clinical-genetic** analysis: variant coordinates, event **type**, consequence, exon, domain, reading-frame **status**, functional **scores** REVEL, MetaLR, and clinical pathogenicity/phenotype labels.

Final dataset size: **11308 variants across 29 features**. For clinical classes in `clinvar_class_simple` there are pathogenic variants (`pathogenic`, 2858), likely pathogenic variants (`likely_pathogenic`, 560), healthy (`benign`, 640), likely healthy (`likely_benign`, 3062), variants with unknown consequences (`vus`, 3251), and other (`other`, 937); by phenotype in `phenotype_group`: Duchenne muscular dystrophy (`DMD`, 7807), Becker muscular dystrophy (`BMD`, 1023), and others (`other`, 2478). For expression, a separate dataset from **Gene Expression Omnibus** is used — **GSE38417**.

Dataset limitations are related to the nature of source databases, as they are aggregated observations with non-uniform depth of clinical annotation, and some fields may be incomplete or ambiguous (for example, `vus`, mixed phenotype and patient condition wording).

Key bias sources: **ClinVar** predominantly publishes clinically interesting variants, **label-noise** in likely/vus categories, population heterogeneity of **gnomAD** frequencies, and dependence of a number of features on quality of primary annotation.

### I thank the administrators of the above-mentioned databases for providing open access to genomic and clinical data.

**Please note** that this repository includes only processed data. Source data remain subject to source-database license terms (see `data/raw/raw_info.md`).

---
### Methods

For ``02_EDA`` we use a single reproducible statistical stack: preliminary normalization and cleaning of categorical and numerical fields (including removal of placeholder categories), then a combination of descriptive analytics (frequencies, shares, cross-tables, CI, entropy) and hypothesis tests selected for data type — Fisher Exact for 2x2 and rare events, chi-square for categorical associations, Mann-Whitney and Kruskal-Wallis tests for non-parametric distribution comparisons, Spearman for monotonic associations, and KS for comparing distribution shapes; results are visualized in Plotly (bar, heatmap, box, scatter, density, lollipop).

In the ML part of the project, a binary classification task is solved for DMD variants: ``pathogenic`` versus ``benign``. The target label is formed from ``clinvar_class_simple`` by a ``likely-inclusive`` scheme: ``pathogenic + likely_pathogenic`` is ``1``, ``benign + likely_benign`` is ``0``, variants of other classes are not included in training. This allows using more labeled examples while preserving the clinical meaning of the target. Before training, deduplication and removal of conflicting duplicates are performed (when the same variant appears with different labels).

As features, structural and functional characteristics of a variant are used: exon number, protein domain, aggregated mutation type, reading-frame status, consequence/variant-type, and also numerical features (``interval_length``, ``aa_pos``, ``revel``, ``meta_lr``, ``allele_freq`` etc). Categorical features are encoded via `OneHotEncoder`, numerical features are scaled and, if necessary, imputed inside the pipeline. It is important to note that fields directly involved in target formation are not used as input features.

Two models are trained: ``LogisticRegression`` (as an interpretable linear basis) and ``RandomForestClassifier`` (as a non-linear ensemble). Logistic regression gives stable interpretation of feature contribution via coefficients, while RandomForest better captures non-linear interactions between molecular characteristics. Both models are trained in a single preprocessing pipeline, making comparison correct and reproducible.

For validation, a group-aware scheme is used so that close variants do not leak between train and test. At cross-validation stage, ``StratifiedGroupKFold`` is applied, and in final evaluation holdout metrics are additionally fixed. Quality is assessed by a set of indicators (``ROC-AUC``, ``accuracy``, ``balanced accuracy``, ``F1``, ``confusion matrix``), and also via sanity-checks (including permutation checks) for leakage and overfitting. Both trained models are saved in ``models/`` as separate ``.pkl`` artifacts.

---
### Obtained model metrics

In the project, quality assessment is split into two levels: **holdout** and **group-aware cross-validation**. For holdout, we calculate four key metrics: ``accuracy``, ``AUROC``, ``F1``, ``balanced accuracy``, so as not to rely on a single number. Logistic regression was the best.

For **CV** (5-fold StratifiedGroupKFold), we compute the same four metrics and look at mean and spread across folds.

Key conclusion for the metrics section: holdout and cross-validation provide a consistent picture without a sharp collapse outside training, and fold spread is small, which indicates pipeline stability. At the same time, the main limitation remains: **lack of external validation**.


| Model | Holdout Accuracy | Holdout ROC-AUC | Holdout F1 | Holdout Balanced Acc | CV Accuracy (mean ± std) | CV ROC-AUC (mean ± std) | CV F1 (mean ± std) | CV Balanced Acc (mean ± std) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| LogReg_tuned | 0.9833 | 0.9950 | 0.9825 | 0.9829 | 0.9815 ± 0.0041 | 0.9965 ± 0.0009 | 0.9813 ± 0.0043 | 0.9815 ± 0.0041 |
| RandomForest_tuned | 0.9822 | 0.9936 | 0.9813 | 0.9817 | 0.9806 ± 0.0030 | 0.9948 ± 0.0008 | 0.9804 ± 0.0031 | 0.9806 ± 0.0029 |

---
### EDA: consistency with the literature

#### 1. Reading-frame rule: confirmed at cohort level
![lc1](assets/lc/lc1.png)
In the analysis, the relation between `frame_status` and phenotype (`DMD/BMD`) reproduces the classical observation: `out-of-frame` variants are much more often associated with a more severe clinical profile, whereas `in-frame` with a milder one. This is consistent with the basic dystrophinopathies model, but we do not make an overclaim at the individual prognosis level for each variant.

**Paper:**  
Koenig et al., 1989 — [PubMed](https://pubmed.ncbi.nlm.nih.gov/2491009/)  
Aartsma-Rus et al., 2006 — [PubMed](https://pubmed.ncbi.nlm.nih.gov/16770791/)


#### 2. Exon hotspot architecture is generally consistent with registries
![lc1](assets/lc/lc2.png)
Distribution of events by exons in our dataset demonstrates clustering in hotspot regions, which is consistent with large DMD databases. This confirms that the variant structure in the project is biologically plausible and does not look like a random sampling artifact.

**Paper:**  
Bladen et al., 2015 (TREAT-NMD, >7000 mutations) — [PubMed](https://pubmed.ncbi.nlm.nih.gov/25604253/)  
Tuffery-Giraud et al., 2009 — [PubMed](https://pubmed.ncbi.nlm.nih.gov/19367636/)


#### 3. Therapeutic skip regions: effect direction matches the literature
In block `02N` (skip45/51/53), we see expected clinical heterogeneity by therapeutically relevant regions, which corresponds to the logic of literature on exon-skipping applicability (including priority of exon 51 for a significant share of patients). Interpretation is careful: this is a population signal, not an evaluation of efficacy of a specific therapy.

**Paper:**  
Aartsma-Rus et al., 2009 — [PubMed](https://pubmed.ncbi.nlm.nih.gov/19156838/)

#### 4. Distal isoforms (Dp140/Dp71): trend agrees with clinical works
In `02M`, observations on `dp140/dp71` and phenotype agree in direction with studies where distal regions and brain isoforms are linked to a more severe/complicated clinical profile. At the same time, we do not transfer conclusions directly to neurocognitive outcomes, because our dataset is not specialized for deep cognitive phenotyping.

**Paper:**  
Moizard et al., 1998 — [PubMed](https://pubmed.ncbi.nlm.nih.gov/9800909/)  
Bardoni et al., 2000 — [PubMed](https://pubmed.ncbi.nlm.nih.gov/10734267/)  
Milic Rasic et al., 2015 — [PubMed](https://pubmed.ncbi.nlm.nih.gov/25937795/)

#### 5. Functional/population features: consistent as supporting evidence, but not as standalone
In `02J/02K`, signals of `REVEL/meta_lr` and population frequencies (`in_gnomad`, `allele_freq`) behave in the expected direction for benign vs pathogenic separation. This fits well with clinical interpretation practice: in-silico and population criteria increase confidence, but by themselves do not replace integrated expert classification.

**Paper:**  
REVEL (Ioannidis et al., 2016) — [DOI](https://doi.org/10.1016/j.ajhg.2016.08.016)  
dbNSFP v3.0 (MetaLR) — [PubMed](https://pubmed.ncbi.nlm.nih.gov/26555599/)  
ACMG/AMP guideline — [PubMed](https://pubmed.ncbi.nlm.nih.gov/25741868/)

---

### Research-gap findings (??)

⚠️ **Important! All hypotheses are applied exclusively to this dataset, and their consequences are not implied as reality (in other words, this project makes no serious claims 😊): verification on other datasets is required.**

#### 1. Rare consequence classes show strong enrichment for pathogenicity
![rg1](assets/rg/rg1.png)

In the consequence block, it is visible that rare consequence categories are associated with a much higher share of pathogenic variants: `pathogenic_fraction = 0.734` versus `0.193` for non-rare, `OR = 11.54` (95% CI `10.37-12.84`, `p < 1e-6`). This is a strong statistical signal, but it should be interpreted correctly as a hypothesis about rare high-risk classes, because ClinVar may have sampling bias toward clinically noticeable variants.

#### 2. condition_raw semantics may hide biologically different subgroups

![rg2](assets/rg/rg2.png)

In the condition block, it was found that not all clinical formulations behave equally with respect to hotspot regions: for example, for one aggregated “neuromuscular disease …” category `OR ≈ 3.00`, `p ≈ 1.0e-9`, whereas for “duchenne muscular dystrophy” the effect is non-significant (`OR ≈ 0.92`, `p ≈ 0.26`). This points to a potential research gap in clinical ontology: coarse-grained textual condition labels may mix subtypes with different genetic architecture.  

#### 3. There is a link between exon position and event length, but the effect is very small

![rg3](assets/rg/rg3.png)

For `exon_num` vs `interval_length`, a statistically significant but weak monotonic relation was obtained: `Spearman rho = 0.072`, `p = 0.000484`, `n = 2328`. This means that a positional gradient exists, but by itself explains a small part of event-length variability; region-specific models and stratification by mutation classes are needed.  

#### 4. Domain scale and pathogenicity: there is a trend, but still insufficient power

![rg4](assets/rg/rg4.png)

At domain level, a positive trend is observed between domain size and share of pathogenic variants (`rho = 0.276`), but there is no statistical significance (`p = 0.133`, `n_domains = 31`). At the same time, domain distribution entropy is high (`4.85 bits`), which confirms structural heterogeneity. The signal may be real, but current power is insufficient for a confident conclusion.  

#### 5. Exceptions to reading-frame rule and meta-inconsistency require a separate QC layer

![rg5](assets/rg/rg5.png)

In exception analysis, enrichment of exceptions in rod domain was found (`OR = 1.25`, 95% CI `1.03–1.51`, `p = 0.0259`), and in the meta-consistency block `1304` potentially inconsistent variants were noted (`mismatch_any`). Part of biology really goes beyond the simple rule framework, and part may be a consequence of annotation noise. External validation is needed.

---
### Project structure

```commandline
DMD_variant_landscape/
├── data/                          # Project data (raw + processed)
│   ├── raw/                       # Raw exports (ClinVar/Ensembl/gnomAD/GEO etc.)
│   │   ├── annotation/            # Raw tables for exons/domains
│   │   └── expression/            # Expression files (e.g. GSE38417 series matrix)
│   └── processed/                 # Prepared tables for analysis/modeling
│       ├── annotation/            # Normalized annotations (exons/domains)
│       ├── DMD_variants_master.csv
│       └── DMD_variants_annotated.csv
│
├── notebooks/                     # Jupyter notebooks by project stages
│   ├── 00_data_preview.ipynb      # Primary QC and data-structure overview
│   ├── 01_annotation_check.ipynb  # Annotation correctness checking
│   ├── 02_EDA/                    # Main EDA block (02A–02O)
│   │   ├── 02A_dataset_integrity.ipynb
│   │   ├── 02B_clinical_structure.ipynb
│   │   ├── ...
│   │   └── 02O_meta_consistency.ipynb
│   ├── 03_model_training.ipynb    # Training/evaluation of ML models
│   └── 04_expression_analysis.ipynb # Expression analysis (GEO GSE38417)
│
├── src/                           # Scripted project pipeline
│   ├── data_preparation.py        # Cleaning/merge of source datasets
│   ├── annotate_variants.py       # Annotation exon/domain/mutation/frame
│   ├── exploratory.py             # Scripted EDA plots/tables
│   ├── modeling.py                # ML pipeline (LogReg + RandomForest, validation)
│   └── utils.py                   # Common functions (cleaning, tests, QC, helper metrics)
│
├── models/                        # Saved models
│   ├── model_logreg.pkl           # Logistic regression
│   ├── model_random_forest.pkl    # RandomForest
│   └── model.pkl                  # Best/current model alias
│
├── figures/                       # Plots and visualizations generated by scripts
│   └── ...
│ 
├── assets/                        # Nice things for README.md
│   └── ...
│ 
├── requirements.txt               # Python dependencies
├── README.md                      # Project description and results
└── LICENSE                        # Repository license
```
---

### Reproducibility

Install Python `3.11+` (recommended `3.12`), create a virtual environment and install dependencies:

```commandline
python -m venv .venv && source .venv/bin/activate && pip install -r requirements.txt
```

Then check that source exports are in `data/raw/` (ClinVar/gnomAD/Ensembl; for expression module — GEO `GSE38417`). After that, run the script pipeline step by step from the repository root:
```commandline
python src/data_preparation.py
python src/annotate_variants.py
python src/exploratory.py
python src/modeling.py
```
If you already downloaded a prepared dataset (annotated), then simply run:

```commandline
python src/exploratory.py
python src/modeling.py
```

For interactive analysis, open Jupyter: `jupyter lab`, then run sequentially `notebooks/00_data_preview.ipynb`, `01_annotation_check.ipynb`, block `02_EDA/`, `03_model_training.ipynb`, and `04_expression_analysis.ipynb`.

The project fixes key reproducibility sources: single `RANDOM_STATE=42` in the ML pipeline, group-aware validation (`StratifiedGroupKFold`), and deterministic class balancing with fixed seed. All dependencies are listed in [requirements.txt](requirements.txt), major artifacts are saved in `models/` (`model_logreg.pkl`, `model_random_forest.pkl`, `model.pkl`) and `figures/` (ROC/feature importance/confusion matrix and EDA/expression plots).

---

### Next steps

The next priority step is **external validation** on an independent set (not ClinVar-centric), to check model transferability beyond the current source and assess how metrics change with another distribution of classes and annotations. Practically, this means repeating the full infer pipeline on an external dataset with the same feature engineering, fixed `.pkl` models, and a separate report on `accuracy / balanced accuracy / F1 / ROC-AUC`, including error analysis by mutation classes and exon regions. This step will close the main risk of dataset-specific performance and provide a more honest assessment of clinical generalizability.

In parallel, it is worth adding **prospective validation** and **probability calibration**. Prospective scenario: the model is tested on “new in time” variants (temporal split), not only on random holdout; this is closer to real-world use in the flow of new interpretations.

---
### References

#### DMD genotype-phenotype, hotspots, domains, isoforms

1. Koenig M, Beggs AH, Moyer M, et al. **The molecular basis for Duchenne versus Becker muscular dystrophy: correlation of severity with type of deletion**. *Am J Hum Genet*. 1989;45(4):498-506. PMID: 2491009.  
   https://pubmed.ncbi.nlm.nih.gov/2491009/

2. Aartsma-Rus A, van Deutekom JCT, Fokkema IFJ, van Ommen GJB, den Dunnen JT. **Entries in the Leiden Duchenne muscular dystrophy mutation database: an overview of mutation types and paradoxical cases that confirm the reading-frame rule**. *Muscle Nerve*. 2006;34(2):135-144. PMID: 16770791.  
   https://pubmed.ncbi.nlm.nih.gov/16770791/

3. Tuffery-Giraud S, Béroud C, Leturcq F, et al. **Genotype-phenotype analysis in 2,405 patients with a dystrophinopathy using the UMD-DMD database**. *Hum Mutat*. 2009;30(6):934-945. doi:10.1002/humu.20976. PMID: 19367636.  
   https://pubmed.ncbi.nlm.nih.gov/19367636/

4. Bladen CL, Salgado D, Monges S, et al. **The TREAT-NMD DMD Global Database: analysis of more than 7,000 Duchenne muscular dystrophy mutations**. *Hum Mutat*. 2015;36(4):395-402. doi:10.1002/humu.22758. PMID: 25604253.  
   https://pubmed.ncbi.nlm.nih.gov/25604253/

5. Matsumura K, Nonaka I, Tomé FMS, et al. **Mild deficiency of dystrophin-associated proteins in Becker muscular dystrophy patients having in-frame deletions in the rod domain of dystrophin**. *Am J Hum Genet*. 1993;53(2):409-416. PMID: 8328458.  
   https://pubmed.ncbi.nlm.nih.gov/8328458/

6. Moizard MP, Billard C, Toutain A, et al. **Are Dp71 and Dp140 brain dystrophin isoforms related to cognitive impairment in Duchenne muscular dystrophy?** *Am J Med Genet*. 1998;80(1):32-41. PMID: 9800909.  
   https://pubmed.ncbi.nlm.nih.gov/9800909/

7. Bardoni A, Felisari G, Sironi M, et al. **Loss of Dp140 regulatory sequences is associated with cognitive impairment in dystrophinopathies**. *Neuromuscul Disord*. 2000;10(3):194-199. doi:10.1016/S0960-8966(99)00108-X. PMID: 10734267.  
   https://pubmed.ncbi.nlm.nih.gov/10734267/

8. Milic Rasic V, Vojinovic D, Pesovic J, et al. **Intellectual ability in the Duchenne muscular dystrophy and dystrophin gene mutation location**. *Balkan J Med Genet*. 2015;17(2):25-35. doi:10.2478/bjmg-2014-0071. PMID: 25937795.  
   https://pubmed.ncbi.nlm.nih.gov/25937795/

#### Exon-skipping and therapy-relevant regions

9. Aartsma-Rus A, Fokkema I, Verschuuren J, et al. **Theoretic applicability of antisense-mediated exon skipping for Duchenne muscular dystrophy mutations**. *Hum Mutat*. 2009;30(3):293-299. doi:10.1002/humu.20918. PMID: 19156838.  
   https://pubmed.ncbi.nlm.nih.gov/19156838/

10. Waldrop MA, Ben Yaou R, Lucas KK, et al. **Clinical Phenotypes of DMD Exon 51 Skip Equivalent Deletions: A Systematic Review**. *J Neuromuscul Dis*. 2020;7(3):217-229. doi:10.3233/JND-200483. PMID: 32417793.  
   https://pubmed.ncbi.nlm.nih.gov/32417793/

#### Variant interpretation, functional predictors, population criteria

11. Richards S, Aziz N, Bale S, et al. **Standards and guidelines for the interpretation of sequence variants (ACMG/AMP)**. *Genet Med*. 2015;17(5):405-424. doi:10.1038/gim.2015.30. PMID: 25741868.  
   https://pubmed.ncbi.nlm.nih.gov/25741868/

12. Ghosh R, Harrison SM, Rehm HL, Plon SE, Biesecker LG. **Updated Recommendation for the Benign Stand Alone ACMG/AMP Criterion (BA1)**. *Hum Mutat*. 2018;39(11):1525-1530. doi:10.1002/humu.23642. PMID: 30311383.  
   https://pubmed.ncbi.nlm.nih.gov/30311383/

13. Ioannidis NM, Rothstein JH, Pejaver V, et al. **REVEL: An Ensemble Method for Predicting the Pathogenicity of Rare Missense Variants**. *Am J Hum Genet*. 2016;99(4):877-885. doi:10.1016/j.ajhg.2016.08.016. PMID: 27666373.  
   https://pubmed.ncbi.nlm.nih.gov/27666373/

14. Liu X, Wu C, Li C, Boerwinkle E. **dbNSFP v3.0: A One-Stop Database of Functional Predictions and Annotations for Human Nonsynonymous and Splice-Site SNVs**. *Hum Mutat*. 2016;37(3):235-241. doi:10.1002/humu.22932. PMID: 26555599.  
   https://pubmed.ncbi.nlm.nih.gov/26555599/

#### Data resources used in the project

15. Landrum MJ, Lee JM, Benson M, et al. **ClinVar: improving access to variant interpretations and supporting evidence**. *Nucleic Acids Res*. 2018;46(D1):D1062-D1067. doi:10.1093/nar/gkx1153.  
   https://academic.oup.com/nar/article/46/D1/D1062/4641904

16. Karczewski KJ, Francioli LC, Tiao G, et al. **The mutational constraint spectrum quantified from variation in 141,456 humans (gnomAD)**. *Nature*. 2020;581:434-443. doi:10.1038/s41586-020-2308-7. PMID: 32461654.  
   https://pubmed.ncbi.nlm.nih.gov/32461654/

17. Harrison PW, Amode MR, Austine-Orimoloye O, et al. **Ensembl 2024**. *Nucleic Acids Res*. 2024;52(D1):D891-D899. doi:10.1093/nar/gkad1049. PMID: 37953337.  
   https://pubmed.ncbi.nlm.nih.gov/37953337/

18. NCBI GEO. **GSE38417: Gene expression data from Duchenne muscular dystrophy patients versus controls**.  
   https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE38417

19. Darras BT, Urion DK, Ghosh PS. **Dystrophinopathies (GeneReviews®)**. Last revision: 2022.  
   https://www.ncbi.nlm.nih.gov/books/NBK1119/

---

### License

This project is distributed under the **MIT** license.

For details, see the LICENSE file.

---

### Project status: 🟩 🟨 🟨 

The technical part of the project is complete! 🎉 

**The project is currently in the results validation stage.**

---

© 2026, Mikhail Kolesnikov (Mikhail Kolesnikov) \
Moscow, Higher School of Economics, Faculty of Computer Science, BSc 

MIT License\
GitHub: https://github.com/curryy77
