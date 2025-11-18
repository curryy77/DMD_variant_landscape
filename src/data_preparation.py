import pandas as pd
from pathlib import Path

# Working with ClinVar DMD Database
def load_and_clean_clinvar(path: str | Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")

    # Removing other genes that may appear in the raw data
    df = df[df["Gene(s)"].str.contains("DMD", na=False)]

    df = df.rename(columns= {
        "VariationID": "var_id",
        "Protein change": "protein_change",
        "Variant type": "var_type",
        "GRCh38Chromosome": "chr38",
        "GRCh38Location": "pos38",
        "dbSNP ID": "rsid",
        "Germline classification": "clinvar_class",
        "Molecular consequence": "clinvar_consequence",
        "Condition(s)": "condition_raw",
    })

    # Turning the genomic position into the numeric type, else NaN
    df["pos38"] = pd.to_numeric(df["pos38"], errors="coerce")

    mapping = {
        "Pathogenic" : "pathogenic",
        "Pathogenic/Likely pathogenic" : "pathogenic",
        "Likely pathogenic" : "likely_pathogenic",
        "Benign" : "benign",
        "Likely benign" : "likely_benign",
        "Benign/Likely benign" : "benign",
        "Uncertain significance" : "vus",
    }

    # Mapping into simple ClinVar class, else - put in the other group
    df["clinvar_class_simple"] = df["clinvar_class"].map(mapping).fillna("other")

    # likely_pathogenic group here was considered to be pathogenic
    df["is_pathogenic"] = df["clinvar_class_simple"].isin(
        ["pathogenic", "likely_pathogenic"]
    )

    # Distributing the data into the phenotypes
    cond = df["condition_raw"].fillna("")
    df["phenotype_group"] = "other"
    df.loc[cond.str.contains("Duchenne", case=False), "phenotype_group"] = "DMD"
    df.loc[cond.str.contains("Becker", case=False), "phenotype_group"] = "BMD"

    # Final columns
    cols = [
        "var_id",
        "rsid",
        "chr38",
        "pos38",
        "protein_change",
        "var_type",
        "clinvar_consequence",
        "clinvar_class",
        "clinvar_class_simple",
        "is_pathogenic",
        "phenotype_group",
        "condition_raw",
    ]

    df = df[cols]
    return df
