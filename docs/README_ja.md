<p align="center">
  <a href="README_en.md">🇬🇧 English</a> |
  <a href="README_ja.md">🇯🇵 日本語</a> |
  <a href="README_fr.md">🇫🇷 Français</a> |
  <a href="README_ru.md">🇷🇺 Русский</a> 
</p>

# DMD遺伝子バリアントスペクトラム解析：ClinVar・Ensembl・gnomAD由来の統合データに基づく探索的研究

![License](https://img.shields.io/badge/license-MIT-green)
![Python](https://img.shields.io/badge/python-3.10%2B-blue)
![Bioinformatics](https://img.shields.io/badge/field-bioinformatics-orange)
![Data](https://img.shields.io/badge/data-ClinVar%20%7C%20gnomAD%20%7C%20Ensembl-purple)

---
### 要旨

このリポジトリは、バイオインフォマティクス情報源から統合したデータを用いた、ヒト **DMD遺伝子** における病原性および非病原性バリアントの研究分析を含みます。研究の目的は、**エクソン**、**タンパク質ドメイン**、**リーディングフレーム状態** などの枠組みで病原性の構造的・機能的な規則性を調べることです。

元データセットのクリーニングと統合を実施し、得られたデータセットに基づいてバリアントの **探索的解析** を行いました。また、**ロジスティック回帰**（LogisticRegression）と **ランダムフォレスト**（RandomForest）に基づいて、DMD遺伝子バリアントの病原性を予測する2つの **予測モデル** を構築しました。さらに、**デュシェンヌ型またはベッカー型筋ジストロフィー** の患者と健常患者における **遺伝子発現解析** も実施し、疾患が **DMD遺伝子** の発現および生体内の **全体的な遺伝子発現** に与える影響を明らかにすることを目的としました。

---

### 目次

- [プロジェクト課題](#プロジェクト課題)
- [主要ポイントとプロジェクトパイプライン](#主要ポイントとプロジェクトパイプライン)
- [データセット](#データセット)
- [方法](#方法)
- [得られたモデル指標](#得られたモデル指標)
- [EDA：文献との整合性](#eda文献との整合性)
- [Research-gap の発見 (??)](#research-gap-の発見-)
- [プロジェクト構造](#プロジェクト構造)
- [再現可能性](#再現可能性)
- [次のステップ](#次のステップ)
- [ライセンス](#ライセンス)
- [プロジェクト状況](#プロジェクト状況----)

---

### プロジェクト課題

**DMD遺伝子** はヒトゲノムで既知最大の遺伝子であり、X染色体上に位置し、**ジストロフィン** タンパク質をコードします。このタンパク質は、筋細胞の **構造的骨格** と **周辺マトリクス** を連結します。DMD遺伝子の病原性バリアントは **ジストロフィン合成の障害** を引き起こし、その結果として **重篤な遺伝性疾患** の群が発症します。その中には、筋力低下・歩行不能化・心筋症の急速進行をもたらす重症型の **デュシェンヌ型筋ジストロフィー**（DMD）と、より軽症で緩徐進行性の筋力低下を特徴とする **ベッカー型筋ジストロフィー**（BMD）があります。

本プロジェクトの課題は、遺伝子内バリアントの **さまざまな特徴** とパターンを研究し、それらが規定する **病原性** と **非病原性** を明らかにすることです。これには、将来的にその情報に基づいて **後続で発見されるバリアントの病原性を予測** する可能性も含まれます。

---

### 主要ポイントとプロジェクトパイプライン

本プロジェクトの **パイプライン** を示す図：

![Project Pipeline2](../assets/pipeline/project_pipeline_ru.png)

**データ準備** 段階では、まず DMD遺伝子バリアントに関する **3つのテーブル** を **ClinVar**・**Ensembl**・**gnomAD** から収集しました。次に、キー（例：``var_id``, ``rsid``）を **統一** し、座標およびカテゴリフィールドの形式を単一標準に **整合** し、**技術的重複** と空カテゴリを **削除** し、欠損値を明示的に処理した上で数値特徴を **正規化** しました。その後、主バリアントテーブルからの left-join による **段階的「結合」** を適用し、各段階で行損失、欠損増加、フィールド整合性を **制御** しました。以上のコードは ``src/data_preparation.py`` にあります。

次にバリアントを **アノテーション** します：``src/annotate_variants.py`` では、スクリプトが **バリアントテーブル** と2つの **参照テーブル**（ジストロフィン蛋白ドメイン境界および DMDエクソン座標）を読み込み、各バリアントに対して **アミノ酸位置** に基づくドメイン、**ゲノム座標** に基づくエクソンを割り当てます。続いて **変異タイプ** と **リーディングフレーム状態** などの派生特徴を構築します：変異タイプは consequence/variant-type 列から集約され、バリアントが **複数ヌクレオチド** レベルで起こり長さが50超である大規模イベントには別ルールを適用し、リーディングフレーム状態は **変異クラス** から導出されます。出力は ``DMD_variants_annotated.csv`` で、これを用いて **後続作業** を行います。

探索的解析は ``notebooks/00_data_preview.ipynb`` と ``notebooks/01_annotation_check.ipynb`` から開始し、その後 ``notebooks/02_EDA`` の各ブロックを **体系的** に実行します。初期段階では ``var_id`` と ``rsid`` の **一意性**、重複、欠損率、``chr`` と ``pos`` の整合、ドメイン・エクソン・リーディングフレーム状態のカバレッジを確認しました。その結果、データは全体として統計解析に適することが分かりましたが、欠損は不均一に分布し、より **複雑な** バリアントや一部の変異クラスに集中するため、欠損を単なるランダムノイズとして扱うことはできませんでした。

続くテストにより、表現型、変異タイプ、リーディングフレーム状態、ドメイン、``consequence``、区間長、``REVEL/meta_lr`` 指標、アレル頻度、さらにエクソンおよびアミノ酸位置に関する **スペクトラムの主要像** が得られました。多くのクロス表・順位比較で構造は **非ランダム** であることが示されました：臨床群は実際に変異タイプと分子文脈で異なり、病原性バリアントはより「急進的」な ``consequence`` と結びつく傾向が強く、集団特徴も臨床的に重要なバリアントで低頻度側に向かうという期待通りの挙動を示しました。

``⭐`` 仮説（文献由来の結果）の検証は、全体として **既知の文献的支柱を確認** しました：reading-frame rule（リーディングフレーム状態と DMD/BMD 表現型の関係）、エクソン hotspot 構造（遠位 45-55 を含む）、同位体/遠位領域（dp140, dp71）の **臨床的関連性**、および skip 領域（skip51を含む）の **治療的フォーカス** です。ドメイン構造的観察と、病原性マーカーとしての REVEL・MetaLR の機能 ``score`` 指標も同様に整合しました。

``📖`` ブロックでは、稀な ``consequence`` 濃縮、エントロピープロファイル、リーディングフレーム規則からの **例外解析**、ドメイン位置異常を示しました。主結果は、アノテーションが不一致または **記述不足** のまま残る **具体的ゾーンの抽出** です。これらの領域は独立データでの **外部バリデーション** と深い再検証を要し、真の新規シグナルである可能性も、この情報源特有のアーティファクトである可能性もあります。

予測モデル構築では、共通パイプラインの注釈特徴（``exon/domain/mutation/frame`` ほか）に基づき、DMD遺伝子バリアントを ``pathogenic`` と ``benign`` に分類する **二値分類** 課題を解きました。データ準備は、ノイズカテゴリのフィルタリング、重複管理、特徴量エンジニアリング、クロスバリデーションでの頑健性確認を重視しました。解釈可能な基礎モデルとして **LogisticRegression** と **RandomForest** の2つを学習し、``.pkl`` に保存しました。評価は ``accuracy``、``AUROC``、``confusion matrix`` などの指標セットで行いました。

品質結果は **強い** ものでした：高い精度水準に到達しました（LogisticRegression で accuracy ``~0.982``、AUROC ``~0.995``）。高い AUROC は、閾値依存指標およびリークリスクと併せて解釈しており、結論の再現可能性を担保します。それでも、独立データセットでの **外部バリデーションが必要** です。

![auroc_1](../assets/auroc.png)

**発現ブロック** では、**Gene Expression Omnibus** の ``GSE38417`` を使用し、**発現マトリクス** と **表現型** を抽出して **DMD** 群と **対照** 群を構築しました。その後、患者群で ``DMD`` と関連マーカーが対照に対してどう振る舞うかを調べるため、**差次解析** を実施しました。

得られた結果として、差次シグナルでは **筋再生遺伝子**（``MYH3/MYH8`` を含む）が際立っており、これはデュシェンヌ型筋ジストロフィーにおける慢性的障害と **再生応答** と整合します。

![expression](../assets/expression.png)

---

### データセット

本プロジェクトでは、公開情報源 **ClinVar**（臨床解釈）、**Ensembl**（ゲノム・転写産物注釈）、**gnomAD**（集団頻度）から収集した **DMD遺伝子** バリアントの統合データセットを使用します。最終作業テーブルは次に保存されています：
```commandline
data/processed/DMD_variants_annotated.csv
```

このテーブルには **臨床遺伝学的** 解析の主要フィールドが含まれます：バリアント座標、イベント **タイプ**、consequence、エクソン、ドメイン、リーディングフレーム **状態**、機能 **スコア** REVEL・MetaLR、および病原性/表現型の臨床ラベル。

最終データサイズは **29特徴に対して 11308 バリアント** です。`clinvar_class_simple` の臨床クラスには、病原性（`pathogenic`, 2858）、おそらく病原性（`likely_pathogenic`, 560）、健常（`benign`, 640）、おそらく健常（`likely_benign`, 3062）、影響不明（`vus`, 3251）、その他（`other`, 937）が含まれます。`phenotype_group` では、デュシェンヌ型筋ジストロフィー（`DMD`, 7807）、ベッカー型筋ジストロフィー（`BMD`, 1023）、その他（`other`, 2478）です。発現については **Gene Expression Omnibus** の **GSE38417** を別途使用します。

このデータセットの制約は元データベースの性質に由来します。すなわち、これらは臨床注釈の深さが不均一な集約観測であり、一部フィールドは不完全または曖昧になり得ます（例：`vus`、表現型と患者状態の混在表記）。

主なバイアス源：**ClinVar** には主に臨床的に注目されるバリアントが掲載されること、likely/vus カテゴリの **label-noise**、**gnomAD** 頻度の集団的不均一性、そして複数特徴が一次注釈品質に依存することです。

### 上記データベース管理者の皆様がゲノム・臨床データをオープンアクセスで提供していることに感謝します。

**注意**：このリポジトリには処理済みデータのみが含まれます。元データには引き続き元データベースのライセンス条件が適用されます（`data/raw/raw_info.md` 参照）。

---
### 方法

``02_EDA`` では、単一の再現可能な統計スタックを使用します：カテゴリ/数値フィールドの前処理正規化とクリーニング（placeholderカテゴリ削除を含む）を行い、その後、記述解析（頻度、割合、クロステーブル、CI、エントロピー）とデータ型に合わせた仮説検定を組み合わせます。2x2 と稀イベントには Fisher Exact、カテゴリ関連にはカイ二乗、分布のノンパラ比較には Mann-Whitney と Kruskal-Wallis、単調関連には Spearman、分布形比較には KS を用い、結果は Plotly（bar, heatmap, box, scatter, density, lollipop）で可視化します。

本プロジェクトの ML 部分では、DMD バリアントの二値分類課題（``pathogenic`` 対 ``benign``）を解きます。ターゲットラベルは ``clinvar_class_simple`` から ``likely-inclusive`` スキームで作成します：``pathogenic + likely_pathogenic`` を ``1``、``benign + likely_benign`` を ``0`` とし、他クラスは学習に含めません。これにより、ターゲットの臨床的意味を保ちながら、より多くのラベル付き例を利用できます。学習前に重複除去と競合重複（同一バリアントに異なるラベル）除去を実施します。

特徴量として、エクソン番号、タンパク質ドメイン、集約変異タイプ、リーディングフレーム状態、consequence/variant-type、および数値特徴（``interval_length``, ``aa_pos``, ``revel``, ``meta_lr``, ``allele_freq`` など）を使用します。カテゴリ特徴は `OneHotEncoder` で符号化し、数値特徴はパイプライン内でスケーリングし必要に応じて補完します。重要な点として、ターゲット形成に直接関与するフィールドは入力特徴として使用しません。

2つのモデルを学習します：``LogisticRegression``（解釈可能な線形ベース）と ``RandomForestClassifier``（非線形アンサンブル）。ロジスティック回帰は係数を通じて特徴寄与を安定に解釈でき、RandomForest は分子特徴間の非線形相互作用をより良く捉えます。両モデルは共通の前処理パイプラインで学習されるため、比較は妥当かつ再現可能です。

バリデーションでは、近いバリアントが train/test 間で漏れないよう、グループを考慮したスキームを使用します。クロスバリデーション段階で ``StratifiedGroupKFold`` を適用し、最終評価では holdout 指標も追加で固定します。品質は ``ROC-AUC``, ``accuracy``, ``balanced accuracy``, ``F1``, ``confusion matrix`` の指標セットで評価し、リークや過学習の sanity-check（置換検定を含む）も行います。学習済み2モデルは ``models/`` に別個の ``.pkl`` アーティファクトとして保存されます。

---
### 得られたモデル指標

本プロジェクトでは品質評価を2レベルに分けます：**holdout**（分離テスト集合）と **group-aware クロスバリデーション**。holdout では ``accuracy``, ``AUROC``, ``F1``, ``balanced accuracy`` の4つの主要指標を算出し、単一数値に依存しません。最良はロジスティック回帰でした。

**CV**（5-fold StratifiedGroupKFold）でも同じ4指標を算出し、fold ごとの平均とばらつきを確認します。

この指標セクションの主要結論は、holdout とクロスバリデーションが学習外での急激な崩壊なしに整合した像を与え、fold 間ばらつきも小さいため、パイプラインが安定していることです。同時に、主な制約は依然として **外部バリデーションの欠如** です。


| Model | Holdout Accuracy | Holdout ROC-AUC | Holdout F1 | Holdout Balanced Acc | CV Accuracy (mean ± std) | CV ROC-AUC (mean ± std) | CV F1 (mean ± std) | CV Balanced Acc (mean ± std) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| LogReg_tuned | 0.9833 | 0.9950 | 0.9825 | 0.9829 | 0.9815 ± 0.0041 | 0.9965 ± 0.0009 | 0.9813 ± 0.0043 | 0.9815 ± 0.0041 |
| RandomForest_tuned | 0.9822 | 0.9936 | 0.9813 | 0.9817 | 0.9806 ± 0.0030 | 0.9948 ± 0.0008 | 0.9804 ± 0.0031 | 0.9806 ± 0.0029 |

---
### EDA：文献との整合性

#### 1. Reading-frame rule：コホートレベルで確認
![lc1](../assets/lc/lc1.png)
解析では、`frame_status` と表現型（`DMD/BMD`）の関係が古典的観察を再現します：`out-of-frame` バリアントはより重い臨床プロファイルと有意に関連し、`in-frame` はより軽いプロファイルと関連します。これは dystrophinopathies の基本モデルと整合しますが、各バリアントの個別予後に対する overclaim は行いません。

**論文：**  
Koenig et al., 1989 — [PubMed](https://pubmed.ncbi.nlm.nih.gov/2491009/)  
Aartsma-Rus et al., 2006 — [PubMed](https://pubmed.ncbi.nlm.nih.gov/16770791/)


#### 2. エクソン hotspot 構造は全体としてレジストリと整合
![lc1](../assets/lc/lc2.png)
本データセットのエクソン別イベント分布はホット領域でのクラスタリングを示し、大規模 DMD データベースと一致します。これは本プロジェクトのバリアント構造が生物学的に妥当であり、ランダムサンプリングのアーティファクトには見えないことを裏付けます。

**論文：**  
Bladen et al., 2015 (TREAT-NMD, >7000 mutations) — [PubMed](https://pubmed.ncbi.nlm.nih.gov/25604253/)  
Tuffery-Giraud et al., 2009 — [PubMed](https://pubmed.ncbi.nlm.nih.gov/19367636/)


#### 3. 治療関連 skip 領域：効果方向は文献と一致
`02N`（skip45/51/53）ブロックでは、治療的に関連のある領域ごとの期待される臨床的不均一性が見られ、exon-skipping applicability 文献のロジック（有意な患者割合に対する exon 51 の優先性を含む）と一致します。解釈は慎重で、これは集団シグナルであり、特定治療の有効性評価ではありません。

**論文：**  
Aartsma-Rus et al., 2009 — [PubMed](https://pubmed.ncbi.nlm.nih.gov/19156838/)

#### 4. 遠位アイソフォーム（Dp140/Dp71）：傾向は臨床研究と一致
`02M` における `dp140/dp71` と phenotype の観察は、遠位領域や脳アイソフォームがより重症/複雑な臨床像と関連する研究と方向一致します。同時に、本データセットは深い認知表現型化に特化していないため、神経認知アウトカムへ直接結論を移しません。

**論文：**  
Moizard et al., 1998 — [PubMed](https://pubmed.ncbi.nlm.nih.gov/9800909/)  
Bardoni et al., 2000 — [PubMed](https://pubmed.ncbi.nlm.nih.gov/10734267/)  
Milic Rasic et al., 2015 — [PubMed](https://pubmed.ncbi.nlm.nih.gov/25937795/)

#### 5. Functional/population 特徴：supporting evidence として整合、ただし standalone ではない
`02J/02K` では、`REVEL/meta_lr` と集団頻度（`in_gnomad`, `allele_freq`）のシグナルが benign vs pathogenic 分離で期待される方向に振る舞います。これは臨床解釈の実務に適合します：in-silico と population 基準は確信を高めますが、それ単独では統合的専門分類を置き換えません。

**論文：**  
REVEL (Ioannidis et al., 2016) — [DOI](https://doi.org/10.1016/j.ajhg.2016.08.016)  
dbNSFP v3.0 (MetaLR) — [PubMed](https://pubmed.ncbi.nlm.nih.gov/26555599/)  
ACMG/AMP guideline — [PubMed](https://pubmed.ncbi.nlm.nih.gov/25741868/)

---

### Research-gap の発見 (??)

⚠️ **重要！ すべての仮説はこのデータセットにのみ適用され、その帰結を現実へ外挿するものではありません（言い換えると、本プロジェクトは強い断定を行いません 😊）：他データセットでの検証が必要です。**

#### 1. 稀な consequence クラスは病原性で強い濃縮を示す
![rg1](../assets/rg/rg1.png)

consequence ブロックでは、稀な consequence カテゴリが非稀カテゴリより病原性割合がはるかに高いことが見られます：`pathogenic_fraction = 0.734` 対 `0.193`、`OR = 11.54`（95% CI `10.37-12.84`, `p < 1e-6`）。これは強い統計シグナルですが、ClinVar に臨床的に目立つバリアントへの sampling bias があり得るため、稀な高リスククラスに関する仮説として解釈するのが適切です。

#### 2. condition_raw の意味論は生物学的に異なるサブグループを隠し得る

![rg2](../assets/rg/rg2.png)

condition ブロックでは、すべての臨床表現が hotspot 領域に対して同様に振る舞うわけではないことが見つかりました。例として、ある集約 “neuromuscular disease …” カテゴリでは `OR ≈ 3.00`, `p ≈ 1.0e-9`、一方 “duchenne muscular dystrophy” では効果は非有意（`OR ≈ 0.92`, `p ≈ 0.26`）です。これは臨床オントロジーにおける潜在的 research-gap を示します：coarse-grained なテキスト condition ラベルは、遺伝構造の異なるサブタイプを混在させる可能性があります。  

#### 3. エクソン位置とイベント長の関係はあるが、効果は非常に小さい

![rg3](../assets/rg/rg3.png)

`exon_num` vs `interval_length` では、統計的に有意だが弱い単調関係が得られました：`Spearman rho = 0.072`, `p = 0.000484`, `n = 2328`。これは位置勾配が存在することを意味しますが、それ単独ではイベント長の差の小部分しか説明しません；領域特異モデルと変異クラス層別化が必要です。  

#### 4. ドメイン規模と病原性：傾向はあるが、現時点で検出力不足

![rg4](../assets/rg/rg4.png)

ドメインレベルでは、ドメインサイズと病原性割合の間に正の傾向（`rho = 0.276`）が観察されますが、統計的有意性はありません（`p = 0.133`, `n_domains = 31`）。同時に、ドメイン分布エントロピーは高く（`4.85 bits`）、構造的不均一性を確認します。シグナルは実在し得ますが、確信的結論には現在の検出力では不足します。  

#### 5. reading-frame rule の例外とメタ不整合には独立した QC 層が必要

![rg5](../assets/rg/rg5.png)

exception 解析では rod ドメインにおける例外の濃縮が検出されました（`OR = 1.25`, 95% CI `1.03–1.51`, `p = 0.0259`）。また meta-consistency ブロックでは `1304` の潜在的不整合バリアント（`mismatch_any`）が記録されました。生物学の一部は単純規則の枠を実際に超え、別の一部は注釈ノイズの結果である可能性があります。外部バリデーションが必要です。

---
### プロジェクト構造

```commandline
DMD_variant_landscape/
├── data/                          # プロジェクトデータ (raw + processed)
│   ├── raw/                       # 元エクスポート (ClinVar/Ensembl/gnomAD/GEO など)
│   │   ├── annotation/            # エクソン/ドメイン用の生テーブル
│   │   └── expression/            # 発現ファイル (例: GSE38417 series matrix)
│   └── processed/                 # 解析/モデリング用の準備済みテーブル
│       ├── annotation/            # 正規化済み注釈 (エクソン/ドメイン)
│       ├── DMD_variants_master.csv
│       └── DMD_variants_annotated.csv
│
├── notebooks/                     # プロジェクト段階ごとの Jupyter ノートブック
│   ├── 00_data_preview.ipynb      # 初期 QC とデータ構造レビュー
│   ├── 01_annotation_check.ipynb  # 注釈の正確性確認
│   ├── 02_EDA/                    # メイン EDA ブロック (02A–02O)
│   │   ├── 02A_dataset_integrity.ipynb
│   │   ├── 02B_clinical_structure.ipynb
│   │   ├── ...
│   │   └── 02O_meta_consistency.ipynb
│   ├── 03_model_training.ipynb    # ML モデルの学習/評価
│   └── 04_expression_analysis.ipynb # 発現解析 (GEO GSE38417)
│
├── src/                           # スクリプト化パイプライン
│   ├── data_preparation.py        # 元データソースのクリーニング/merge
│   ├── annotate_variants.py       # exon/domain/mutation/frame 注釈
│   ├── exploratory.py             # スクリプト EDA グラフ/テーブル
│   ├── modeling.py                # ML パイプライン (LogReg + RandomForest, validation)
│   └── utils.py                   # 共通関数 (クリーニング, テスト, QC, helper 指標)
│
├── models/                        # 保存済みモデル
│   ├── model_logreg.pkl           # ロジスティック回帰
│   ├── model_random_forest.pkl    # RandomForest
│   └── model.pkl                  # 最良/現行モデルのエイリアス
│
├── figures/                       # スクリプト生成の図表・可視化
│   └── ...
│ 
├── assets/                        # README.md 用アセット
│   └── ...
│ 
├── requirements.txt               # Python 依存関係
├── README.md                      # プロジェクト説明と結果
└── LICENSE                        # リポジトリライセンス
```
---

### 再現可能性

Python `3.11+`（推奨 `3.12`）をインストールし、仮想環境を作成して依存関係をインストールしてください：

```commandline
python -m venv .venv && source .venv/bin/activate && pip install -r requirements.txt
```

次に、`data/raw/` に元エクスポート（ClinVar/gnomAD/Ensembl；expression モジュール用に GEO `GSE38417`）があることを確認してください。その後、リポジトリルートから段階的にスクリプトパイプラインを実行します：
```commandline
python src/data_preparation.py
python src/annotate_variants.py
python src/exploratory.py
python src/modeling.py
```
すでに準備済みデータセット（annotated）を取得済みの場合は、次だけで十分です：

```commandline
python src/exploratory.py
python src/modeling.py
```

対話的解析では Jupyter を開いてください：`jupyter lab`。その後 `notebooks/00_data_preview.ipynb`、`01_annotation_check.ipynb`、`02_EDA/` ブロック、`03_model_training.ipynb`、`04_expression_analysis.ipynb` を順に実行します。

プロジェクトは再現可能性の主要要素を固定しています：ML パイプラインの単一 `RANDOM_STATE=42`、group-aware バリデーション（`StratifiedGroupKFold`）、固定 seed による決定論的クラスバランス。依存関係は [requirements.txt](/Users/franceballin/PycharmProjects/DMD_variant_landscape/requirements.txt) に列挙され、主要アーティファクトは `models/`（`model_logreg.pkl`, `model_random_forest.pkl`, `model.pkl`）と `figures/`（ROC/feature importance/confusion matrix と EDA/expression 図）に保存されます。

---

### 次のステップ

次の優先ステップは、独立セット（ClinVar-centric ではない）での **外部バリデーション** です。これにより、現在のソース外へのモデル移植可能性を検証し、クラス分布や注釈分布が異なるときに指標がどう変化するかを評価します。実務的には、同一 feature engineering、固定 `.pkl` モデル、`accuracy / balanced accuracy / F1 / ROC-AUC` の別レポート（変異クラス別・エクソン領域別の誤り解析を含む）で、外部データセット上で infer パイプライン全体を再実行することを意味します。このステップは dataset-specific performance の主要リスクを閉じ、臨床的一般化可能性のより誠実な評価を与えます。

並行して **prospective validation** と **確率校正** を追加すべきです。Prospective シナリオでは、モデルをランダム holdout だけでなく「時間的に新しい」バリアント（temporal split）で評価します；これは新規解釈が流入する実運用により近いです。

---
### 参考文献

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

### ライセンス

このプロジェクトは **MIT** ライセンスで配布されます。

詳細は LICENSE ファイルを参照してください。

---

### プロジェクト状況：🟩 🟨 🟨 

プロジェクトの技術部分は完了しました！ 🎉 

**プロジェクトは現在、結果検証段階にあります。**

---

© 2026, Mikhail Kolesnikov (ミハイル・コレスニコフ) \
Moscow, Higher School of Economics, Faculty of Computer Science, BSc 

MIT License\
GitHub: https://github.com/curryy77
