import sys
from pathlib import Path

import joblib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.base import clone
from sklearn.compose import ColumnTransformer
from sklearn.ensemble import RandomForestClassifier
from sklearn.dummy import DummyClassifier
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    ConfusionMatrixDisplay,
    accuracy_score,
    balanced_accuracy_score,
    f1_score,
    roc_auc_score,
    roc_curve,
    confusion_matrix,
)
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.utils import ensure_column, normalize_cat


# Constant paths
PROCESSED_DATA = Path("../data/processed")
FIGURES = Path("../figures")
MODELS = Path("../models")

RANDOM_STATE = 42
VALIDATION_PERMUTATIONS = 10
VALIDATION_BOOTSTRAP = 300

CATEGORICAL_COLS = [
    "clinvar_class_simple",
    "mutation_type_group",
    "domain",
    "frame_status",
    "clinvar_consequence",
    "var_type",
    "rsid",
]

NUMERIC_SOURCE_COLS = [
    "exon",
    "interval_length",
    "aa_pos",
    "revel",
    "meta_lr",
    "allele_freq",
    "hemizygote_count",
    "homozygote_count",
    "pos",
]

FEATURES = [
    "exon_num",
    "frame_bin",
    "domain_clean",
    "mutation_type_clean",
    "consequence_clean",
    "var_type_clean",
    "domain_group",
    "exon_bin",
    "interval_length_num",
    "aa_pos_num",
    "revel_num",
    "meta_lr_num",
    "allele_freq_num",
    "hemizygote_count_num",
    "homozygote_count_num",
    "pos_num",
    "log_interval_length",
    "log_allele_freq",
    "log_hemizygote",
    "log_homozygote",
    "is_distal",
    "is_hotspot_distal",
    "is_hotspot_proximal",
    "mut_is_splice",
    "mut_is_nonsense",
    "mut_is_frameshift",
    "mut_is_large",
    "missing_domain",
    "missing_revel",
    "missing_meta_lr",
    "missing_af",
    "missing_aa_pos",
    "rsid_missing",
]

# Load the annotated CSV
def load_annotated(path: Path | str = PROCESSED_DATA / "DMD_variants_annotated.csv") -> pd.DataFrame:
    df = pd.read_csv(path)
    return df

# Preparing likely-inclusive target and base numeric fields
def prepare_modeling_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    work = df.copy()

    for col in CATEGORICAL_COLS:
        ensure_column(work, col)
        work[col] = work[col].apply(normalize_cat)

    for col in NUMERIC_SOURCE_COLS:
        ensure_column(work, col)
        work[col + "_num"] = pd.to_numeric(work[col], errors="coerce")

    work["frame_bin"] = work["frame_status"].map({"in-frame": 0, "out-of-frame": 1})
    work["rsid_missing"] = work["rsid"].isna().astype(float)

    work["target_class"] = pd.Series(np.nan, index=work.index, dtype="object")

    work.loc[
        work["clinvar_class_simple"].isin(["pathogenic", "likely_pathogenic"]),
        "target_class",
    ] = "pathogenic"

    work.loc[
        work["target_class"].isna()
        & work["clinvar_class_simple"].isin(["benign", "likely_benign"]),
        "target_class",
    ] = "benign"

    full_df = work[work["target_class"].notna()].copy()
    full_df["is_pathogenic"] = (full_df["target_class"] == "pathogenic").astype(int)

    return full_df


def _normalized_key_token(series: pd.Series, uppercase: bool = False) -> pd.Series:
    out = series.astype("string").str.strip()
    out = out.replace({"": pd.NA, "nan": pd.NA, "None": pd.NA, "none": pd.NA, "<NA>": pd.NA})
    if uppercase:
        out = out.str.upper()
    return out


def build_variant_group_id(df: pd.DataFrame) -> pd.Series:
    tmp = df.copy()
    for col in ["chr", "pos", "ref", "alt", "var_id", "rsid"]:
        ensure_column(tmp, col)

    chr_key = _normalized_key_token(tmp["chr"], uppercase=True).str.replace("^CHR", "", regex=True)
    pos_key = pd.to_numeric(tmp["pos"], errors="coerce").round().astype("Int64").astype("string")
    ref_key = _normalized_key_token(tmp["ref"], uppercase=True)
    alt_key = _normalized_key_token(tmp["alt"], uppercase=True)
    var_id_key = _normalized_key_token(tmp["var_id"])
    rsid_key = _normalized_key_token(tmp["rsid"], uppercase=True)

    group_id = pd.Series(pd.NA, index=tmp.index, dtype="string")
    has_cpra = chr_key.notna() & pos_key.notna() & ref_key.notna() & alt_key.notna()
    group_id.loc[has_cpra] = (
        "CPRA:"
        + chr_key.loc[has_cpra]
        + ":"
        + pos_key.loc[has_cpra]
        + ":"
        + ref_key.loc[has_cpra]
        + ">"
        + alt_key.loc[has_cpra]
    )

    missing = group_id.isna()
    has_cpos = missing & chr_key.notna() & pos_key.notna()
    group_id.loc[has_cpos] = "CPOS:" + chr_key.loc[has_cpos] + ":" + pos_key.loc[has_cpos]

    missing = group_id.isna()
    group_id.loc[missing & var_id_key.notna()] = "VARID:" + var_id_key.loc[missing & var_id_key.notna()]

    missing = group_id.isna()
    group_id.loc[missing & rsid_key.notna()] = "RSID:" + rsid_key.loc[missing & rsid_key.notna()]

    missing = group_id.isna()
    group_id.loc[missing] = pd.Series(tmp.index[missing], index=tmp.index[missing]).map(lambda x: f"ROW:{x}")

    return group_id.astype(str)


def deduplicate_training_pool(full_df: pd.DataFrame):
    tmp = full_df.copy()
    tmp["group_id"] = build_variant_group_id(tmp)

    group_summary = tmp.groupby("group_id")["is_pathogenic"].agg(["size", "nunique"])
    conflict_groups = group_summary[group_summary["nunique"] > 1].index
    n_conflict_groups = int(len(conflict_groups))
    n_conflict_rows = int(tmp["group_id"].isin(conflict_groups).sum())

    tmp = tmp.loc[~tmp["group_id"].isin(conflict_groups)].copy()
    n_before_dedup = len(tmp)
    tmp = tmp.drop_duplicates(subset=["group_id"], keep="first").reset_index(drop=True)
    n_after_dedup = len(tmp)

    report = {
        "n_conflict_groups_removed": n_conflict_groups,
        "n_conflict_rows_removed": n_conflict_rows,
        "n_same_label_duplicates_collapsed": int(n_before_dedup - n_after_dedup),
    }
    return tmp, report


# Balance classes to 50/50
def balance_binary(full_df: pd.DataFrame) -> pd.DataFrame:
    path_df = full_df[full_df["is_pathogenic"] == 1]
    ben_df = full_df[full_df["is_pathogenic"] == 0]

    min_n = min(len(path_df), len(ben_df))

    if min_n == 0:
        raise ValueError("One class is empty after target assignment")

    model_df = pd.concat(
        [
            path_df.sample(n=min_n, random_state=RANDOM_STATE),
            ben_df.sample(n=min_n, random_state=RANDOM_STATE),
        ],
        ignore_index=True,
    ).sample(frac=1.0, random_state=RANDOM_STATE).reset_index(drop=True)

    return model_df


# Feature engineering
def build_feature_matrix(model_df: pd.DataFrame):
    feat = model_df.copy()

    feat["domain_clean"] = feat["domain"].fillna("Missing_domain")
    feat["mutation_type_clean"] = feat["mutation_type_group"].fillna("Missing_mutation")
    feat["consequence_clean"] = feat["clinvar_consequence"].fillna("Missing_consequence")
    feat["var_type_clean"] = feat["var_type"].fillna("Missing_var_type")

    feat["is_distal"] = feat["exon_num"].ge(45).astype(float)
    feat["is_hotspot_distal"] = feat["exon_num"].between(45, 55, inclusive="both").astype(float)
    feat["is_hotspot_proximal"] = feat["exon_num"].between(3, 9, inclusive="both").astype(float)

    feat["exon_bin"] = pd.cut(feat["exon_num"], bins=[0, 9, 44, 55, 79], include_lowest=True)
    feat["exon_bin"] = feat["exon_bin"].astype(str).replace("nan", "Missing_exon_bin")

    feat["domain_group"] = np.where(
        feat["domain_clean"].str.contains("spectrin", case=False, na=False),
        "rod",
        np.where(
            feat["domain_clean"].str.contains(
                "interaction|actin|binding|ch1|ch2|ww|cysteine",
                case=False,
                regex=True,
                na=False,
            ),
            "binding",
            np.where(
                feat["domain_clean"].str.contains("disordered", case=False, na=False),
                "disordered",
                "other_domain",
            ),
        ),
    )

    feat["mut_is_splice"] = feat["mutation_type_clean"].str.contains("splice", case=False, na=False).astype(float)
    feat["mut_is_nonsense"] = feat["mutation_type_clean"].str.contains("nonsense", case=False, na=False).astype(float)
    feat["mut_is_frameshift"] = feat["mutation_type_clean"].str.contains("frameshift", case=False, na=False).astype(float)
    feat["mut_is_large"] = feat["mutation_type_clean"].str.contains("large", case=False, na=False).astype(float)

    feat["log_interval_length"] = np.log1p(feat["interval_length_num"])
    feat["log_allele_freq"] = np.log1p(feat["allele_freq_num"])
    feat["log_hemizygote"] = np.log1p(feat["hemizygote_count_num"])
    feat["log_homozygote"] = np.log1p(feat["homozygote_count_num"])

    feat["missing_domain"] = feat["domain"].isna().astype(float)
    feat["missing_revel"] = feat["revel_num"].isna().astype(float)
    feat["missing_meta_lr"] = feat["meta_lr_num"].isna().astype(float)
    feat["missing_af"] = feat["allele_freq_num"].isna().astype(float)
    feat["missing_aa_pos"] = feat["aa_pos_num"].isna().astype(float)

    for col in FEATURES:
        ensure_column(feat, col)

    X = feat[FEATURES].copy()
    y = feat["is_pathogenic"].astype(int)

    return feat, X, y


def split_train_test_by_group(X: pd.DataFrame, y: pd.Series, groups: pd.Series):
    splitter = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=RANDOM_STATE)
    train_idx, test_idx = next(splitter.split(X, y, groups=groups))
    X_train = X.iloc[train_idx].copy()
    X_test = X.iloc[test_idx].copy()
    y_train = y.iloc[train_idx].copy()
    y_test = y.iloc[test_idx].copy()
    return X_train, X_test, y_train, y_test


# Build preprocess pipeline and models
def build_preprocess(X: pd.DataFrame) -> ColumnTransformer:
    num_cols = [c for c in X.columns if pd.api.types.is_numeric_dtype(X[c])]
    cat_cols = [c for c in X.columns if c not in num_cols]

    preprocess = ColumnTransformer(
        transformers=[
            (
                "num",
                Pipeline(
                    [
                        ("imputer", SimpleImputer(strategy="median")),
                        ("scaler", StandardScaler()),
                    ]
                ),
                num_cols,
            ),
            (
                "cat",
                Pipeline(
                    [
                        ("imputer", SimpleImputer(strategy="most_frequent")),
                        ("onehot", OneHotEncoder(handle_unknown="ignore")),
                    ]
                ),
                cat_cols,
            ),
        ]
    )

    return preprocess


def build_logreg_model(preprocess: ColumnTransformer) -> Pipeline:
    model = Pipeline(
        [
            ("preprocess", preprocess),
            (
                "model",
                LogisticRegression(
                    max_iter=8000,
                    solver="saga",
                    random_state=RANDOM_STATE,
                    C=4.5705630998,
                    penalty="l1",
                    class_weight=None,
                ),
            ),
        ]
    )
    return model


def build_rf_model(preprocess: ColumnTransformer) -> Pipeline:
    model = Pipeline(
        [
            ("preprocess", preprocess),
            (
                "model",
                RandomForestClassifier(
                    random_state=RANDOM_STATE,
                    n_jobs=-1,
                    n_estimators=300,
                    min_samples_split=8,
                    min_samples_leaf=4,
                    max_features=0.3,
                    max_depth=12,
                    class_weight="balanced_subsample",
                    bootstrap=True,
                ),
            ),
        ]
    )
    return model


# Evaluate two models
def evaluate_models(models: dict, X_test: pd.DataFrame, y_test: pd.Series):
    rows = []
    pred_store = {}

    for name, model in models.items():
        y_proba = model.predict_proba(X_test)[:, 1]
        y_pred = (y_proba >= 0.5).astype(int)
        pred_store[name] = {"proba": y_proba, "pred": y_pred}

        rows.append(
            {
                "model": name,
                "accuracy": accuracy_score(y_test, y_pred),
                "balanced_accuracy": balanced_accuracy_score(y_test, y_pred),
                "f1": f1_score(y_test, y_pred),
                "auc": roc_auc_score(y_test, y_proba),
            }
        )

    results_df = (
        pd.DataFrame(rows)
        .sort_values(["accuracy", "auc"], ascending=False)
        .reset_index(drop=True)
    )

    return results_df, pred_store


def _compute_metrics(y_true: pd.Series, y_proba: np.ndarray) -> dict:
    y_pred = (y_proba >= 0.5).astype(int)
    return {
        "accuracy": accuracy_score(y_true, y_pred),
        "balanced_accuracy": balanced_accuracy_score(y_true, y_pred),
        "f1": f1_score(y_true, y_pred),
        "auc": roc_auc_score(y_true, y_proba),
    }


def evaluate_models_train_test(
    models: dict,
    X_train: pd.DataFrame,
    y_train: pd.Series,
    X_test: pd.DataFrame,
    y_test: pd.Series,
) -> pd.DataFrame:
    rows = []
    for name, model in models.items():
        train_proba = model.predict_proba(X_train)[:, 1]
        test_proba = model.predict_proba(X_test)[:, 1]
        train_metrics = _compute_metrics(y_train, train_proba)
        test_metrics = _compute_metrics(y_test, test_proba)

        rows.append(
            {
                "model": name,
                "split": "train",
                **train_metrics,
            }
        )
        rows.append(
            {
                "model": name,
                "split": "test",
                **test_metrics,
            }
        )

    return pd.DataFrame(rows)


def evaluate_baselines(
    X_train: pd.DataFrame,
    y_train: pd.Series,
    X_test: pd.DataFrame,
    y_test: pd.Series,
) -> pd.DataFrame:
    rows = []
    for strategy in ["most_frequent", "stratified"]:
        baseline = DummyClassifier(strategy=strategy, random_state=RANDOM_STATE)
        baseline.fit(X_train, y_train)
        y_proba = baseline.predict_proba(X_test)[:, 1]
        rows.append(
            {
                "model": f"Dummy_{strategy}",
                **_compute_metrics(y_test, y_proba),
            }
        )

    return pd.DataFrame(rows)


def permutation_test(
    estimator: Pipeline,
    X_train: pd.DataFrame,
    y_train: pd.Series,
    X_test: pd.DataFrame,
    y_test: pd.Series,
    observed_auc: float,
    observed_accuracy: float,
    n_permutations: int = VALIDATION_PERMUTATIONS,
) -> tuple[pd.DataFrame, dict]:
    rng = np.random.default_rng(RANDOM_STATE)
    y_train_arr = y_train.to_numpy()
    rows = []

    for i in range(n_permutations):
        perm_model = clone(estimator)
        y_perm = rng.permutation(y_train_arr)
        perm_model.fit(X_train, y_perm)
        perm_proba = perm_model.predict_proba(X_test)[:, 1]
        perm_pred = (perm_proba >= 0.5).astype(int)
        rows.append(
            {
                "iteration": i + 1,
                "accuracy": accuracy_score(y_test, perm_pred),
                "auc": roc_auc_score(y_test, perm_proba),
            }
        )

    perm_df = pd.DataFrame(rows)
    p_value_auc = (float((perm_df["auc"] >= observed_auc).sum()) + 1.0) / (n_permutations + 1.0)
    p_value_acc = (float((perm_df["accuracy"] >= observed_accuracy).sum()) + 1.0) / (n_permutations + 1.0)

    summary = {
        "n_permutations": int(n_permutations),
        "perm_auc_mean": float(perm_df["auc"].mean()),
        "perm_auc_std": float(perm_df["auc"].std(ddof=1)),
        "perm_accuracy_mean": float(perm_df["accuracy"].mean()),
        "perm_accuracy_std": float(perm_df["accuracy"].std(ddof=1)),
        "p_value_auc": float(p_value_auc),
        "p_value_accuracy": float(p_value_acc),
    }
    return perm_df, summary


def bootstrap_ci(y_true: pd.Series, y_proba: np.ndarray, n_bootstrap: int = VALIDATION_BOOTSTRAP) -> dict:
    rng = np.random.default_rng(RANDOM_STATE)
    y_true_arr = y_true.to_numpy()
    n = len(y_true_arr)
    auc_vals = []
    acc_vals = []

    for _ in range(n_bootstrap):
        idx = rng.integers(0, n, n)
        y_b = y_true_arr[idx]
        if np.unique(y_b).size < 2:
            continue
        p_b = y_proba[idx]
        pred_b = (p_b >= 0.5).astype(int)
        auc_vals.append(roc_auc_score(y_b, p_b))
        acc_vals.append(accuracy_score(y_b, pred_b))

    auc_arr = np.array(auc_vals, dtype=float)
    acc_arr = np.array(acc_vals, dtype=float)

    if len(auc_arr) == 0 or len(acc_arr) == 0:
        return {
            "n_bootstrap_valid": 0,
            "auc_ci_low": np.nan,
            "auc_ci_high": np.nan,
            "accuracy_ci_low": np.nan,
            "accuracy_ci_high": np.nan,
        }

    return {
        "n_bootstrap_valid": int(len(auc_arr)),
        "auc_ci_low": float(np.quantile(auc_arr, 0.025)),
        "auc_ci_high": float(np.quantile(auc_arr, 0.975)),
        "accuracy_ci_low": float(np.quantile(acc_arr, 0.025)),
        "accuracy_ci_high": float(np.quantile(acc_arr, 0.975)),
    }


def save_permutation_plot(perm_df: pd.DataFrame, observed_auc: float, observed_accuracy: float):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

    axes[0].hist(perm_df["auc"], bins=12, color="#4C78A8", alpha=0.8)
    axes[0].axvline(observed_auc, color="#E45756", linestyle="--", linewidth=2, label=f"Observed AUC={observed_auc:.3f}")
    axes[0].set_title("Permutation Test: AUC")
    axes[0].set_xlabel("AUC")
    axes[0].set_ylabel("Count")
    axes[0].legend()

    axes[1].hist(perm_df["accuracy"], bins=12, color="#72B7B2", alpha=0.8)
    axes[1].axvline(
        observed_accuracy,
        color="#E45756",
        linestyle="--",
        linewidth=2,
        label=f"Observed Acc={observed_accuracy:.3f}",
    )
    axes[1].set_title("Permutation Test: Accuracy")
    axes[1].set_xlabel("Accuracy")
    axes[1].set_ylabel("Count")
    axes[1].legend()

    plt.tight_layout()
    plt.savefig(FIGURES / "model_permutation_test.png", dpi=200)
    plt.close()


# Save plots and artifacts
def save_roc_plot(pred_store: dict, y_test: pd.Series):
    plt.figure(figsize=(8, 6))
    for name, vals in pred_store.items():
        fpr, tpr, _ = roc_curve(y_test, vals["proba"])
        auc = roc_auc_score(y_test, vals["proba"])
        plt.plot(fpr, tpr, linewidth=2, label=f"{name} (AUC={auc:.3f})")

    plt.plot([0, 1], [0, 1], "--", color="gray")
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title("ROC Curve: likely-inclusive label scheme")
    plt.legend(loc="lower right")
    plt.tight_layout()
    plt.savefig(FIGURES / "roc_curve.png", dpi=200)
    plt.close()


def save_feature_importance(best_model: Pipeline, best_model_name: str):
    feature_names = best_model.named_steps["preprocess"].get_feature_names_out()

    if "LogReg" in best_model_name:
        importances = np.abs(best_model.named_steps["model"].coef_[0])
    else:
        importances = best_model.named_steps["model"].feature_importances_

    imp_df = pd.DataFrame(
        {
            "feature": feature_names,
            "importance": importances,
        }
    ).sort_values("importance", ascending=False)

    plot_df = imp_df.head(25).iloc[::-1]

    plt.figure(figsize=(10, 8))
    plt.barh(plot_df["feature"], plot_df["importance"])
    plt.xlabel("Importance")
    plt.ylabel("Feature")
    plt.title(f"Feature importance: {best_model_name}")
    plt.tight_layout()
    plt.savefig(FIGURES / "feature_importance.png", dpi=200)
    plt.close()

    return imp_df


def save_confusion(best_pred, y_test, best_model_name: str):
    cm = confusion_matrix(y_test, best_pred)
    fig, ax = plt.subplots(figsize=(5, 4))
    ConfusionMatrixDisplay(cm, display_labels=["benign", "pathogenic"]).plot(ax=ax, colorbar=False)
    plt.title(f"Confusion matrix: {best_model_name}")
    plt.tight_layout()
    plt.savefig(FIGURES / "confusion_matrix.png", dpi=200)
    plt.close()


# Main modeling pipeline
def run_modeling_pipeline():
    MODELS.mkdir(parents=True, exist_ok=True)

    df = load_annotated()
    full_df = prepare_modeling_dataframe(df)
    dedup_df, dedup_report = deduplicate_training_pool(full_df)
    model_df = balance_binary(dedup_df)

    feat, X, y = build_feature_matrix(model_df)
    groups = feat["group_id"].astype(str)

    X_train, X_test, y_train, y_test = split_train_test_by_group(X, y, groups)

    preprocess = build_preprocess(X)

    logreg_model = build_logreg_model(preprocess)
    rf_model = build_rf_model(preprocess)

    models = {
        "LogReg_tuned": logreg_model,
        "RandomForest_tuned": rf_model,
    }

    for model in models.values():
        model.fit(X_train, y_train)

    results_df, pred_store = evaluate_models(models, X_test, y_test)

    best_model_name = results_df.iloc[0]["model"]
    best_model = models[best_model_name]
    best_pred = pred_store[best_model_name]["pred"]
    best_proba = pred_store[best_model_name]["proba"]
    observed_accuracy = accuracy_score(y_test, best_pred)
    observed_auc = roc_auc_score(y_test, best_proba)

    # Validation artifacts
    train_test_metrics_df = evaluate_models_train_test(models, X_train, y_train, X_test, y_test)
    baseline_df = evaluate_baselines(X_train, y_train, X_test, y_test)
    perm_df, perm_summary = permutation_test(
        best_model,
        X_train,
        y_train,
        X_test,
        y_test,
        observed_auc=observed_auc,
        observed_accuracy=observed_accuracy,
        n_permutations=VALIDATION_PERMUTATIONS,
    )
    ci_summary = bootstrap_ci(y_test, best_proba, n_bootstrap=VALIDATION_BOOTSTRAP)

    train_groups = set(groups.loc[X_train.index].astype(str))
    test_groups = set(groups.loc[X_test.index].astype(str))
    overlap_n = int(len(train_groups & test_groups))

    # Save both models + best alias
    joblib.dump(models["LogReg_tuned"], MODELS / "model_logreg.pkl")
    joblib.dump(models["RandomForest_tuned"], MODELS / "model_random_forest.pkl")
    joblib.dump(best_model, MODELS / "model.pkl")

    print("Full labeled rows:", len(full_df))
    print("After dedup rows:", len(dedup_df))
    print("Balanced rows:", len(model_df))
    print("Train rows:", len(X_train), "Test rows:", len(X_test))
    print("Conflict groups removed:", dedup_report["n_conflict_groups_removed"])
    print("Conflict rows removed:", dedup_report["n_conflict_rows_removed"])
    print("Same-label duplicates collapsed:", dedup_report["n_same_label_duplicates_collapsed"])
    print("Group overlap train/test:", overlap_n)
    print("Best model:", best_model_name)
    print("Best accuracy:", round(observed_accuracy, 4))
    print("Best balanced accuracy:", round(balanced_accuracy_score(y_test, best_pred), 4))
    print("Best AUC:", round(observed_auc, 4))
    print()
    print("Train/test metrics:")
    print(train_test_metrics_df.to_string(index=False))
    print()
    print("Dummy baselines:")
    print(baseline_df.to_string(index=False))
    print("Permutation p-value (AUC):", round(perm_summary["p_value_auc"], 6))
    print("Permutation p-value (accuracy):", round(perm_summary["p_value_accuracy"], 6))
    print("95% CI AUC:", f"[{ci_summary['auc_ci_low']:.4f}, {ci_summary['auc_ci_high']:.4f}]")
    print("95% CI accuracy:", f"[{ci_summary['accuracy_ci_low']:.4f}, {ci_summary['accuracy_ci_high']:.4f}]")

    return results_df


# By launching main, we obtain model artifacts in ../models and ../figures
if __name__ == "__main__":
    run_modeling_pipeline()
