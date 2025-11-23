import pandas as pd
from pathlib import Path

# Constant paths
RAW_DATA = Path("../data/raw")
PROCESSED_DATA = Path("../data/processed")

# Working with ClinVar DMD Dataset - the main dataset in the analysis
def load_and_clean_clinvar(path: str | Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")

    # Removing other genes that may appear in the raw data
    df = df[df["Gene(s)"].str.contains("DMD", na=False)]

    df = df.rename(columns= {
        "VariationID": "var_id",
        "dbSNP ID": "rsid",
        "GRCh37Chromosome": "chr37",
        "GRCh37Location": "pos37",
        "GRCh38Chromosome": "chr38",
        "GRCh38Location": "pos38",
        "Protein change": "protein_change",
        "Variant type": "var_type",
        "Germline classification": "clinvar_class",
        "Molecular consequence": "clinvar_consequence",
        "Condition(s)": "condition_raw",
    })

    # Fill NaN in chromosomes column with blank values
    chr38 = df["chr38"].fillna("")
    chr37 = df["chr37"].fillna("")

    # Unify chr38 and chr37 (chr38, else chr37), replace artifacts to X
    chr_unified = chr38.where(chr38 != "", chr37)
    chr_unified = chr_unified.replace({"X|X": "X", "X||X": "X"})

    df["chr"] = chr_unified

    # Fill NaN in positions with blank values
    pos38 = df["pos38"].fillna("").astype(str)
    pos37 = df["pos37"].fillna("").astype(str)

    # Unify pos38 and pos37 (pos38, else pos37)
    pos_unified = pos38.where(pos38 != "", pos37)

    # Work with the intervals in pos_unified
    parts = pos_unified.str.split("-", n = 1, expand=True)
    start_str = parts[0]
    # If it is not an interval, then the string is None
    end_str = parts[1] if parts.shape[1] > 1 else None

    # Save the interval length for intervals in the dataframe
    start = pd.to_numeric(start_str, errors="coerce")
    if end_str is not None:
        end = pd.to_numeric(end_str, errors="coerce")
        interval_length = end - start
    else:
        interval_length = pd.Series(index=df.index, dtype="float64") # for points, it will be NaN

    # Save the columns
    df["pos"] = start
    df["interval_length"] = interval_length

    # Drop the remaining NaN values for chr and pos
    df = df.dropna(subset=["chr", "pos"])

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
        "chr",
        "pos",
        "interval_length",
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

# Working with gnomAD DMD Database
def load_and_clean_gnomad(path: str | Path) -> pd.DataFrame:
    df = pd.read_csv(path)

    df = df.rename(columns= {
        "Chromosome": "chr",
        "Position": "pos",
        "Reference": "ref",
        "Alternate": "alt",
        "Allele Frequency": "allele_freq",
        "ClinVar Variation ID": "var_id",
        "Hemizygote Count": "hemizygote_count",
        "Homozygote Count": "homozygote_count",
    })

    # Turning the genomic position in the numeric type, same as in ClinVar
    df["pos"] = pd.to_numeric(df["pos"], errors="coerce")

    # Returns true if the allele frequency is not NaN
    df["in_gnomad"] = ~df["allele_freq"].isna()

    # Checks if the allele frequency is high enough
    df["af_high"] = df["allele_freq"] > 1e-4

    cols = [
        "var_id",
        "chr",
        "pos",
        "ref",
        "alt",
        "allele_freq",
        "in_gnomad",
        "af_high",
        "hemizygote_count",
        "homozygote_count",
    ]

    df = df[cols]
    return df

# Working with Ensembl DMD Database
def load_and_clean_ensembl(path: str | Path) -> pd.DataFrame:
    df = pd.read_csv(path)

    df = df.rename(columns= {
        "Variant ID": "rsid",
        "Location": "location",
        "Conseq. Type": "ensembl_consequence",
        "AA": "aa_change",
        "AA coord": "aa_pos",
        "REVEL": "revel",
        "MetaLR": "meta_lr",
        "Mutation Assessor": "mutation_assessor",
    })

    # Split the location into two columns - chromosome and the position in it
    chr_pos = df["location"].str.split(":",  expand=True)
    df["chr"] = chr_pos[0]
    df["pos"] = pd.to_numeric(chr_pos[1], errors="coerce")

    # Deal with the duplicates on rsid
    df = df.drop_duplicates(subset=["rsid", "pos"])

    cols = [
        "rsid",
        "chr",
        "pos",
        "ensembl_consequence",
        "aa_change",
        "aa_pos",
        "revel",
        "meta_lr",
        "mutation_assessor",
    ]

    df = df[cols]
    return df

# Preparing two CSV-tables from GSE38417 series matrix, one with metadata, another with the expression
def prepare_gse38417():
    series_path = RAW_DATA / "expression/GSE38417_series_matrix.txt"

    # Names for the output csv-tables
    metadata_out = PROCESSED_DATA / "GSE38417_samples.csv"
    expression_out = PROCESSED_DATA / "GSE38417_expression.csv"

    # Preparing the metadata (GSE38417_samples.csv)
    meta_lines = {}
    with series_path.open(encoding="utf-8", errors="ignore") as f:
        for raw_line in f:
            line = raw_line.rstrip("\n")

            # We take the !Sample_ strings for metadata
            if line.startswith("!Sample_"):
                parts = line.split("\t")

                key = parts[0].replace("!Sample_", "")
                values = parts[1:]
                meta_lines[key] = values

            # Signal string
            if line.startswith("!series_matrix_table_begin"):
                break

    # Transform the dictionary into the dataframe
    meta_df = pd.DataFrame(meta_lines)

    # Save the safe columns
    cols = [
        "geo_accession",
        "title",
        "source_name_ch1",
        "characteristics_ch1"
    ]

    meta_df = meta_df[cols]

    # Transpose the dataframe to match the expression
    meta_df = meta_df.T

    # Convert the metadata to csv
    meta_df.to_csv(metadata_out, index=False)

    # Preparing the expression (GSE38417_expression.csv), via pandas
    expr_df = pd.read_csv(
        series_path,
        sep = "\t",
        comment = "!",
        index_col = 0
    )

    # Convert the expression to csv
    expr_df.to_csv(expression_out, index=False)

    return expr_df, meta_df

# Merge the dataframes obtained
def prepare_master():
    clinvar_path = RAW_DATA / "ClinVar_DMD_raw.tsv"
    gnomad_path = RAW_DATA / "gnomAD_DMD_raw.csv"
    ensembl_path = RAW_DATA / "Ensembl_DMD_raw.csv"

    # The output file name
    master_out = PROCESSED_DATA / "DMD_variants_master.csv"

    # Loading dataframes
    cvar = load_and_clean_clinvar(clinvar_path)
    gnom = load_and_clean_gnomad(gnomad_path)
    ensm = load_and_clean_ensembl(ensembl_path)

    # Firstly, we merge ClinVar <= gnomAD by var_id
    merged = cvar.merge(
        gnom,
        on=["var_id", "chr", "pos"],
        how="left",
    )

    # Secondly, we do ClinVar <= Ensembl by rsid
    merged = merged.merge(
        ensm,
        on=["rsid", "chr", "pos"],
        how="left",
    )

    merged.to_csv(master_out, index=False)

    return merged

# By launching main, we obtain all the csv-data in the ../data/processed
if __name__ == "__main__":
    prepare_master()
    prepare_gse38417()
