import pandas as pd
from pathlib import Path

PROCESSED_DATA = Path("../data/processed")

def load_data():
    variants = pd.read_csv(PROCESSED_DATA / "DMD_variants_master.csv")
    domains = pd.read_csv(PROCESSED_DATA / "annotation/dystrophin_domains.csv")
    exons = pd.read_csv(PROCESSED_DATA / "annotation/DMD_exons.csv")
    return variants, domains, exons


def annotate_domain(variants: pd.DataFrame, domains: pd.DataFrame) -> pd.Series:
    if "aa_pos" not in variants.columns:
        return pd.Series(index=variants.index, dtype="object")

    dom_series = pd.Series(index=variants.index, dtype="object")

    for _, row in domains.iterrows():
        mask = variants["aa_pos"].between(row["start_aa"], row["end_aa"], inclusive="both")
        dom_series.loc[mask] = row["domain_name"]

    return dom_series


def annotate_exon(variants: pd.DataFrame, exons: pd.DataFrame) -> pd.Series:
    exon_series = pd.Series(index=variants.index, dtype="float")

    for _, ex in exons.iterrows():
        mask = (
        (variants["chr"].astype(str) == str(ex["chr"])) &
        (variants["pos"].between(ex["start"], ex["end"], inclusive="both"))
        )
        exon_series.loc[mask] = ex["number"]

    return exon_series

# I will continue this later

# def infer_mutation_type_group(row) -> str:
#    cons = (str(row.get("clinvar_consequence", "")) + " " +
#            str(row.get("ensembl_consequence", ""))).lower()
