import re
from typing import Iterable

import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, fisher_exact, kruskal, mannwhitneyu, spearmanr


PLACEHOLDER_CATS = {
    "",
    "na",
    "n/a",
    "nan",
    "none",
    "null",
    "unknown",
    "other",
    "not provided",
    "not specified",
    "see cases",
    ".",
    "-",
    "--",
}

MUTATION_TO_CONSEQUENCE_KEYWORD = {
    "missense": "missense",
    "nonsense": "nonsense",
    "frameshift": "frameshift",
    "splice": "splice",
}

FRAME_EXPECTED_BY_MUTATION = {
    "frameshift": "out-of-frame",
    "large_deletion": "out-of-frame",
    "nonsense": "out-of-frame",
    "splice": "out-of-frame",
    "missense": "in-frame",
}

LOW_RISK_CONSEQUENCE_TERMS = ["synonymous", "intron variant", "utr"]
HIGH_RISK_CONSEQUENCE_TERMS = ["nonsense", "frameshift", "splice donor", "splice acceptor"]


def ensure_column(df: pd.DataFrame, col: str, fill_value=np.nan) -> None:
    if col not in df.columns:
        df[col] = fill_value


def normalize_category(value):
    if pd.isna(value):
        return np.nan
    text = str(value).strip()
    if text.lower() in PLACEHOLDER_CATS:
        return np.nan
    return text


def normalize_cat(value):
    return normalize_category(value)


def clean_condition(value):
    if pd.isna(value):
        return np.nan
    text = str(value).lower().strip()
    text = text.replace("|", ";")
    text = re.sub(r"[^a-z0-9; /_\-]+", " ", text)
    text = re.sub(r"\s+", " ", text).strip()
    if text == "":
        return np.nan
    if text in PLACEHOLDER_CATS:
        return np.nan
    return text


def to_bool(value):
    if pd.isna(value):
        return False
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    text = str(value).strip().lower()
    return text in {"1", "true", "t", "yes", "y"}


def odds_ratio_ci(a, b, c, d):
    vals = np.array([a, b, c, d], dtype=float)
    if (vals == 0).any():
        vals = vals + 0.5
    or_value = (vals[0] * vals[3]) / (vals[1] * vals[2])
    se = np.sqrt(np.sum(1 / vals))
    low = np.exp(np.log(or_value) - 1.96 * se)
    high = np.exp(np.log(or_value) + 1.96 * se)
    return or_value, low, high


def fisher_bool(data: pd.DataFrame, exposure_col: str, outcome_col: str):
    tmp = data[[exposure_col, outcome_col]].dropna().copy()
    tmp[exposure_col] = tmp[exposure_col].astype(bool)
    tmp[outcome_col] = tmp[outcome_col].astype(bool)

    tab = pd.crosstab(tmp[exposure_col], tmp[outcome_col]).reindex(
        index=[False, True],
        columns=[False, True],
        fill_value=0,
    )

    odds, p = fisher_exact(tab.values)
    a = int(tab.loc[True, True])
    b = int(tab.loc[True, False])
    c = int(tab.loc[False, True])
    d = int(tab.loc[False, False])
    or_ci = odds_ratio_ci(a, b, c, d)

    return tab, odds, p, or_ci


def chi2_table(data: pd.DataFrame, col1: str, col2: str):
    tmp = data[[col1, col2]].dropna()
    tab = pd.crosstab(tmp[col1], tmp[col2])
    if tab.shape[0] < 2 or tab.shape[1] < 2:
        return tab, np.nan, np.nan, np.nan
    chi2, p, dof, _ = chi2_contingency(tab)
    return tab, chi2, p, dof


def mann_whitney_bool(data: pd.DataFrame, bool_col: str, value_col: str):
    tmp = data[[bool_col, value_col]].dropna()
    g1 = tmp.loc[tmp[bool_col].astype(bool), value_col]
    g0 = tmp.loc[~tmp[bool_col].astype(bool), value_col]
    if len(g1) == 0 or len(g0) == 0:
        return g1, g0, np.nan, np.nan
    stat, p = mannwhitneyu(g1, g0, alternative="two-sided")
    return g1, g0, stat, p


def kruskal_group(data: pd.DataFrame, group_col: str, value_col: str):
    tmp = data[[group_col, value_col]].dropna()
    grouped = [(k, g[value_col].values) for k, g in tmp.groupby(group_col) if len(g) > 1]
    names = [k for k, _ in grouped]
    groups = [v for _, v in grouped]
    if len(groups) < 2:
        return names, groups, np.nan, np.nan
    stat, p = kruskal(*groups)
    return names, groups, stat, p


def spearman_xy(data: pd.DataFrame, x_col: str, y_col: str):
    tmp = data[[x_col, y_col]].dropna()
    if len(tmp) < 3:
        return tmp, np.nan, np.nan
    rho, p = spearmanr(tmp[x_col], tmp[y_col])
    return tmp, rho, p


def bh_adjust(pvalues: Iterable[float]):
    pvalues = np.array(list(pvalues), dtype=float)
    n = len(pvalues)
    if n == 0:
        return np.array([], dtype=float)

    order = np.argsort(pvalues)
    ranked = np.empty(n, dtype=float)
    prev = 1.0

    for i in range(n - 1, -1, -1):
        idx = order[i]
        rank = i + 1
        value = pvalues[idx] * n / rank
        prev = min(prev, value)
        ranked[idx] = prev

    return np.clip(ranked, 0, 1)


def ecdf_xy(values: Iterable[float]):
    arr = np.sort(np.array(list(values), dtype=float))
    if len(arr) == 0:
        return np.array([]), np.array([])
    y = np.arange(1, len(arr) + 1) / len(arr)
    return arr, y


def mut_cons_mismatch(
    row,
    mapping: dict | None = None,
):
    mapping = mapping or MUTATION_TO_CONSEQUENCE_KEYWORD
    mt = row.get("mutation_type")
    cons = str(row.get("consequence", "")).lower()
    if mt in mapping:
        return mapping[mt] not in cons
    return False


def frame_mismatch(
    row,
    expected: dict | None = None,
):
    expected = expected or FRAME_EXPECTED_BY_MUTATION
    mt = row.get("mutation_type")
    fr = row.get("frame")
    if mt in expected:
        return fr != expected[mt]
    return False


def path_cons_mismatch(
    row,
    low_risk_terms: list[str] | None = None,
    high_risk_terms: list[str] | None = None,
):
    low_risk_terms = low_risk_terms or LOW_RISK_CONSEQUENCE_TERMS
    high_risk_terms = high_risk_terms or HIGH_RISK_CONSEQUENCE_TERMS

    cons = str(row.get("consequence", "")).lower()
    pathogenic = bool(row.get("pathogenic", False))

    if pathogenic and any(t in cons for t in low_risk_terms):
        return True
    if (not pathogenic) and any(t in cons for t in high_risk_terms):
        return True
    return False


def prepare_eda_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    d = df.copy()

    ensure_column(d, "is_pathogenic", False)
    ensure_column(d, "phenotype_group")
    ensure_column(d, "mutation_type_group")
    ensure_column(d, "frame_status")
    ensure_column(d, "domain")
    ensure_column(d, "clinvar_consequence")
    ensure_column(d, "exon")
    ensure_column(d, "interval_length")
    ensure_column(d, "aa_pos")
    ensure_column(d, "revel")
    ensure_column(d, "meta_lr")
    ensure_column(d, "allele_freq")
    ensure_column(d, "hemizygote_count")
    ensure_column(d, "homozygote_count")
    ensure_column(d, "pos")
    ensure_column(d, "rsid")
    ensure_column(d, "condition_raw")
    ensure_column(d, "in_gnomad")

    d["pathogenic"] = d["is_pathogenic"].fillna(False).astype(bool)
    d["phenotype"] = d["phenotype_group"].apply(normalize_category)
    d["mutation_type"] = d["mutation_type_group"].apply(normalize_category)
    d["frame"] = d["frame_status"].apply(normalize_category)
    d["domain_clean"] = d["domain"].apply(normalize_category)
    d["consequence"] = d["clinvar_consequence"].apply(normalize_category)

    d["exon_num"] = pd.to_numeric(d["exon"], errors="coerce")
    d["interval_length_num"] = pd.to_numeric(d["interval_length"], errors="coerce")
    d["aa_pos_num"] = pd.to_numeric(d["aa_pos"], errors="coerce")
    d["revel_num"] = pd.to_numeric(d["revel"], errors="coerce")
    d["meta_lr_num"] = pd.to_numeric(d["meta_lr"], errors="coerce")
    d["allele_freq_num"] = pd.to_numeric(d["allele_freq"], errors="coerce")
    d["hemizygote_num"] = pd.to_numeric(d["hemizygote_count"], errors="coerce")
    d["homozygote_num"] = pd.to_numeric(d["homozygote_count"], errors="coerce")
    d["pos_num"] = pd.to_numeric(d["pos"], errors="coerce")

    d["missing_domain"] = d["domain"].isna()
    d["missing_rsid"] = d["rsid"].isna() | d["rsid"].astype(str).str.strip().eq("")

    d["condition_clean"] = d["condition_raw"].apply(clean_condition)
    d["condition_main"] = (
        d["condition_clean"].str.split(";").str[0].str.strip().apply(normalize_category)
    )

    d["is_hotspot_distal"] = d["exon_num"].between(45, 55, inclusive="both")
    d["is_hotspot_proximal"] = d["exon_num"].between(3, 9, inclusive="both")
    d["is_distal_region"] = d["exon_num"] >= 45

    d["dp140_status"] = pd.Series(np.nan, index=d.index, dtype="object")
    d.loc[d["exon_num"] <= 44, "dp140_status"] = "dp140_preserved"
    d.loc[d["exon_num"] >= 51, "dp140_status"] = "dp140_affected"
    d.loc[d["exon_num"].between(45, 50, inclusive="both"), "dp140_status"] = "dp140_indeterminate"

    d["is_dp140_affected"] = d["dp140_status"].eq("dp140_affected")
    d["is_dp71_affected"] = d["exon_num"] >= 63

    d["skip51_region"] = d["exon_num"].eq(51)
    d["skip53_region"] = d["exon_num"].eq(53)
    d["skip45_region"] = d["exon_num"].eq(45)

    d["in_gnomad_bool"] = d["in_gnomad"].apply(to_bool)

    return d


def add_consistency_flags(data: pd.DataFrame) -> pd.DataFrame:
    tmp = data.copy()
    tmp["mismatch_mut_cons"] = tmp.apply(mut_cons_mismatch, axis=1)
    tmp["mismatch_frame_mut"] = tmp.apply(frame_mismatch, axis=1)
    tmp["mismatch_path_cons"] = tmp.apply(path_cons_mismatch, axis=1)
    tmp["mismatch_any"] = tmp[["mismatch_mut_cons", "mismatch_frame_mut", "mismatch_path_cons"]].any(axis=1)
    return tmp
