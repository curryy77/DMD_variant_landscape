<p align="center">
  <a href="README_en.md">🇬🇧 English</a> |
  <a href="README_ja.md">🇯🇵 日本語</a> |
  <a href="README_fr.md">🇫🇷 Français</a> |
  <a href="README_ru.md">🇷🇺 Русский</a> 
</p>

# Analyse des variantes de la DMD

![License](https://img.shields.io/badge/license-MIT-green)
![Python](https://img.shields.io/badge/python-3.10%2B-blue)
![Bioinformatics](https://img.shields.io/badge/field-bioinformatics-orange)
![Data](https://img.shields.io/badge/data-ClinVar%20%7C%20gnomAD%20%7C%20Ensembl-purple)

---
### À propos

Analyse exploratoire des variants pathogènes et non pathogènes du **gène DMD** à l'aide d'ensembles de données annotés sur les variants. L'objectif est d'étudier les schémas structurels et fonctionnels de la pathogénicité au niveau des **exons**, des **domaines protéiques** et du **statut du cadre de lecture**.

---

### Ensembles de données
L'analyse utilise un ensemble de données annotées de **variants du gène DMD** (issu de la fusion des bases de données **ClinVar**, **Ensembl** et **gnomAD**) contenant la position génomique, le numéro d'exon, le domaine protéique, la classification de pathogénicité, la classe de phénotype (DMD ou BMD) et le statut du cadre de lecture (en phase ou hors phase). Les variants sont prétraités et stockés dans :
```commandline
data/processed/DMD_variants_annotated.csv
```

---
### Analyses réalisées (à ce jour)

#### 1. Répartition des variants au niveau des exons
Nous avons tout d'abord analysé la répartition des variants pathogènes et non pathogènes sur l'ensemble des **79 exons du gène DMD**.

![Exon distribution](../figures/0XX_exon_distribution.png)

Deux visualisations complémentaires sont générées : **le nombre de variants par exon** et **la proportion de variants pathogènes par exon**.

On constate une forte concentration de variants dans les exons 42 à 56, qui correspondent à la région de délétion fréquente du gène DMD, décrite dans de nombreuses études (Monaco et al., 1988 ; Aartsma-Rus et al., 2006 ; Bladen et al., 2015).

Cette région recouvre le domaine central en forme de tige de la dystrophine, qui contient de multiples motifs répétitifs de type spectrine.

#### 2. Pathogénicité au niveau du domaine
Les variants ont été cartographiés (cartographie fournie par UniProt) par rapport aux **domaines protéiques** fonctionnels de la dystrophine. Pour chaque domaine, nous avons calculé le nombre total de variants ainsi que la proportion de variants pathogènes.

![Domain distribution](../figures/1XX_domain_distribution.png)

Les domaines présentant la plus forte fraction pathogène comprennent le **domaine de liaison à l'actine**, la **région de liaison à la SNTB1** et l'**interaction avec la syntrophine**. Ces domaines correspondent à des **interfaces d'interaction critiques** nécessaires au rôle structurel de la dystrophine au sein du complexe glycoprotéique de la dystrophine. On s'attend à ce que les mutations dans ces régions perturbent plus fréquemment la fonction protéique et présentent une pathogénicité plus élevée (Ervasti & Campbell, 1993 ; Koenig et al., 1988 ; Blake et al., 2002).

#### 3. Analyse de la règle du cadre de lecture

La règle classique du cadre de lecture, décrite dans (Monaco et al., 1988) dans le contexte des dystrophinopathies, a été vérifiée. Selon cette règle, les mutations hors cadre ont tendance à entraîner des cas graves de DMD, tandis que les mutations dans le cadre se traduisent par des cas plus bénins de dystrophie musculaire de Becker.

![Reading frame rule](../figures/2XX_reading_frame_rule.png)

Un graphique présentant le nombre de variants par phénotype et par statut de cadre de lecture, ainsi qu'un graphique illustrant la répartition des fractions phénotypiques, a été généré. 

Les résultats confirment la tendance attendue : les mutations hors cadre sont fortement surreprésentées dans le DMD, tandis que les mutations dans le cadre présentent une proportion plus élevée de BMD. On observe toutefois des exceptions. Ces exceptions peuvent s'expliquer par la perturbation de domaines fonctionnels essentiels, par des effets d'épissage alternatif et par l'instabilité structurelle de la dystrophine tronquée.

---
### Interprétation biologique

Les analyses mettent en évidence plusieurs caractéristiques importantes du paysage mutationnel de la DMD : **les zones à forte concentration de mutations correspondent à des régions structurellement répétitives** du gène, les **domaines de liaison fonctionnels** présentent **des taux de pathogénicité élevés**, la **règle de cadre de lecture de Monaco** est **largement confirmée**, mais n'est pas absolue.

---
### TODO

Les extensions prévues de l'analyse comprennent : **des tests d'enrichissement statistique**, **une analyse du rapport de cotes** pour la pathogénicité des domaines, **la modélisation** de la densité des mutations au niveau des exons, **l'intégration avec les données structurelles** des domaines de la dystrophine, ainsi qu'une comparaison avec **les fréquences des variants dans la population**.

---
### Data Sources

L'analyse intègre des informations sur les variants provenant de plusieurs bases de données publiques :

- **ClinVar** - annotations cliniques des variants
- **gnomAD** - fréquences alléliques dans la population
- **Ensembl** - annotation génomique et structure des transcrits
- **GEO (GSE38417)** - ensemble de données sur l'expression génique

Nous remercions les responsables des bases de données susmentionnées de nous avoir donné libre accès aux données génomiques et cliniques.

---
### Licence

Ce projet est publié sous licence MIT.

Consultez le fichier LICENSE pour plus de détails.

---

### État du projet : 🟨 🟨 🟨 
Ce référentiel contient un pipeline d'analyse exploratoire du paysage mutationnel de la DMD. 

**Le projet est actuellement en cours de développement !** Ce fichier README.md est provisoire et sera modifié une fois le projet terminé.

---

© 2026, Mikhail Kolesnikov (Mikhaïl Kolesnikov)\
Moscow, Higher School of Economics, Faculty of Computer Science, BSc 

MIT License\
GitHub: https://github.com/curryy77