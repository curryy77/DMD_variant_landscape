from __future__ import annotations

# Some preliminary imports
import pandas as pd
from pathlib import Path
from scipy.stats import chi2_contingency, fisher_exact
import plotly.express as px
import kaleido

# Constant paths
PROCESSED_DATA = Path("../data/processed")
FIGURES = Path("../figures")

# Loads the annotated CSV
def load_annotated(path: Path | str = PROCESSED_DATA / "DMD_variants_annotated.csv") -> pd.DataFrame:
    df = pd.read_csv(path)
    return df

# Preparing function, ensuring exon and exon position is numeric, phenotype class exists
def prepare_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()

    # If there's no phenotype_class
    if "phenotype_class" not in df.columns and "phenotype_group" in df.columns:
        df["phenotype_class"] = df["phenotype_group"].fillna("other")

    # Turn exon number into numeric value
    if "exon" in df.columns:
        df["exon"] = pd.to_numeric(df["exon"], errors="coerce")

    return df

# Main plot #1: Histogram chart of pathogenic vs non-pathogenic variants by exon
def save_exon_distrib_plot(df: pd.DataFrame) -> pd.DataFrame:
    tmp = df.dropna(subset=["exon"]).copy()
    tmp["exon"] = tmp["exon"].astype(int)

    exon_counts = (
        tmp.groupby(["exon", "is_pathogenic"])
        .size()
        .unstack(fill_value=0)
        .rename(columns={True: "pathogenic", False: "non_pathogenic"})
        .sort_index()
        .reset_index()
    )

    exon_long = exon_counts.melt(
        id_vars=["exon"],
        value_vars= [c for c in ["pathogenic", "non_pathogenic"] if c in exon_counts.columns],
        var_name="variant_class",
        value_name="count"
    )

    fig = px.bar(
        exon_long,
        x="exon",
        y="count",
        color="variant_class",
        title="Pathogenic vs non-pathogenic variants by exon",
        labels={
            "exon": "Exon number",
            "count": "Number of variants",
            "variant_class": "Variant class",
        },

        # Order of the legend
        category_orders={
            "variant_class": ["pathogenic", "non_pathogenic"]
        },

        # Custom colors
        color_discrete_map={
            "pathogenic": "#d62728",
            "non_pathogenic": "#1f77b4",
        }
    )

    # X axis from 1 to 79
    fig.update_xaxes(
        tickmode="linear",
        dtick=5,
        range=[0.5, 79.5]
    )

    fig.write_image(FIGURES / "0XX_exon_distribution.png", scale=2)

    return exon_counts

# Main function
def run_exploratory() -> None:
    df = load_annotated()
    df = prepare_dataframe(df)

    exon_table = save_exon_distrib_plot(df)

    print(f"Saved figures to: {FIGURES}")
    print(exon_table)

# By launching main, we obtain all the figures in the ../figures
if __name__ == "__main__":
    run_exploratory()






