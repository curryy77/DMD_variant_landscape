import pandas as pd
from pathlib import Path

PROCESSED_DATA = Path("../data/processed")

# Load the master-csv and the annotation data
def load_data():
    variants = pd.read_csv(PROCESSED_DATA / "DMD_variants_master.csv")
    domains = pd.read_csv(PROCESSED_DATA / "annotation/dystrophin_domains.csv")
    exons = pd.read_csv(PROCESSED_DATA / "annotation/DMD_exons.csv")
    return variants, domains, exons

# Annotating dystrophin domains
def annotate_domain(variants: pd.DataFrame, domains: pd.DataFrame) -> pd.Series:
    if "aa_pos" not in variants.columns:
        return pd.Series(index=variants.index, dtype="object")

    dom_series = pd.Series(index=variants.index, dtype="object")

    for _, row in domains.iterrows():
        mask = variants["aa_pos"].between(row["start_aa"], row["end_aa"], inclusive="both")
        dom_series.loc[mask] = row["domain_name"]

    return dom_series

# Annotating DMD exons
def annotate_exon(variants: pd.DataFrame, exons: pd.DataFrame) -> pd.Series:
    exon_series = pd.Series(index=variants.index, dtype="float")

    for _, ex in exons.iterrows():
        mask = (
        (variants["chr"].astype(str) == str(ex["chr"])) &
        (variants["pos"].between(ex["start"], ex["end"], inclusive="both"))
        )
        exon_series.loc[mask] = ex["number"]

    return exon_series

# Classify mutation groups
def infer_mutation_type_group(row) -> str:
    cons = (str(row.get("clinvar_consequence", "")) + " " +
    str(row.get("ensembl_consequence", ""))).lower()
    var_type = str(row.get("var_type", "")).lower()
    interval = row.get("interval_length", pd.NA)

    # Separate large deletions or duplications
    if pd.notna(interval) and interval > 50:
        return "large_deletion"

    if "frameshift" in cons or "frameshift" in var_type:
        return "frameshift"
    if "stop gained" in cons or "nonsense" in var_type:
        return "nonsense"
    if "splice" in cons:
        return "splice"
    if "missense" in cons:
        return "missense"

    return "other"

# Classify frame status
def infer_frame_status(row) -> str:
    mt = row["mutation_type_group"]

    if mt in ["frameshift", "large_deletion", "nonsense", "splice"]:
        return "out-of-frame"
    if mt in ["missense"]:
        return "in-frame"
    return "unknown"

# Main function of annotating
def annotate_variants():
    variants, domains, exons = load_data()

    # turn aa_pos to numeric
    if "aa_pos" in variants.columns:
        variants["aa_pos"] = pd.to_numeric(variants["aa_pos"], errors="coerce")

    variants["domain"] = annotate_domain(variants, domains)
    variants["exon"] = annotate_exon(variants, exons)

    variants["mutation_type_group"] = variants.apply(infer_mutation_type_group, axis=1)
    variants["frame_status"] = variants.apply(infer_frame_status, axis=1)

    out_path = PROCESSED_DATA / "DMD_variants_annotated.csv"
    variants.to_csv(out_path, index=False)

    return variants

# By launching main, we obtain all the csv-data in the ../data/processed
if __name__ == "__main__":
    annotate_variants()


