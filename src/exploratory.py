from __future__ import annotations

# Some preliminary imports
import pandas as pd
from pathlib import Path
from scipy.stats import chi2_contingency, fisher_exact
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.graph_objects as go
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

# Main plot #1: Histograms of pathogenic vs non-pathogenic variants by exon
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

    # Case if we have only pathogenic/non-pathogenic variants
    if "pathogenic" not in exon_counts:
        exon_counts["pathogenic"] = 0
    if "non_pathogenic" not in exon_counts:
        exon_counts["non_pathogenic"] = 0

    # Define total amount of variants
    exon_counts["total"] = exon_counts["pathogenic"] + exon_counts["non_pathogenic"]

    # Define pathogenic fraction
    exon_counts["pathogenic_fraction"] = exon_counts["pathogenic"] / exon_counts["total"]

    fig = make_subplots(
        rows=1,
        cols=2,
        subplot_titles=[
            "Pathogenic vs non-pathogenic variants by exon",
            "Fraction of pathogenic variants by exon"
        ]
    )

    # Figure 1.1 - pathogenic variants
    fig.add_trace(
        go.Bar(
            x=exon_counts["exon"],
            y=exon_counts["pathogenic"],
            name="pathogenic",
            marker_color="#d62728"
        ),
        row=1,
        col=1
    )

    # Figure 1.2 - non-pathogenic variants
    fig.add_trace(
        go.Bar(
            x=exon_counts["exon"],
            y=exon_counts["non_pathogenic"],
            name="non_pathogenic",
            marker_color="#1f77b4"
        ),
        row=1,
        col=1
    )

    # Figure 2 - pathogenic fraction
    fig.add_trace(
        go.Bar(
            x=exon_counts["exon"],
            y=exon_counts["pathogenic_fraction"],
            name="pathogenic fraction",
            marker_color="#4f4f4f",
            showlegend=False
        ),
        row=1,
        col=2
    )

    fig.update_xaxes(
        tickmode="linear",
        dtick=5,
        range=[0.5, 79.5],
        title="Exon number"
    )

    fig.update_yaxes(title="Number of variants", row=1, col=1)
    fig.update_yaxes(title="Fraction pathogenic", row=1, col=2)

    fig.update_layout(
        barmode="stack",
        width=1400,
        height=500
    )

    fig.write_image(FIGURES / "0XX_exon_distribution.png", scale=2)

    return exon_counts

# Main plot #2: Histograms of pathogenic vs non-pathogenic variants by domain
def save_domain_distrib_plot(df: pd.DataFrame) -> pd.DataFrame:
    tmp = df.dropna(subset=["domain"]).copy()

    domain_counts = (
        tmp.groupby(["domain", "is_pathogenic"])
        .size()
        .unstack(fill_value=0)
        .rename(columns={True: "pathogenic", False: "non_pathogenic"})
        .reset_index()
    )

    # Case if we have only pathogenic/non-pathogenic variants
    if "pathogenic" not in domain_counts.columns:
        domain_counts["pathogenic"] = 0
    if "non_pathogenic" not in domain_counts.columns:
        domain_counts["non_pathogenic"] = 0

    # Define total amount of variants
    domain_counts["total"] = (
        domain_counts["pathogenic"] + domain_counts["non_pathogenic"]
    )

    # Define pathogenic fraction
    domain_counts["pathogenic_fraction"] = (
        domain_counts["pathogenic"] / domain_counts["total"]
    )

    # Sort by pathogenic fraction
    domain_counts = domain_counts.sort_values(
        ["pathogenic_fraction", "total"],
        ascending=[False, False]
    ).reset_index(drop=True)

    # Save only top 25 domains by pathogenicity
    top_domains = domain_counts.head(25).copy()
    x_order = top_domains["domain"].tolist()

    fig = make_subplots(
        rows=1,
        cols=2,
        subplot_titles=[
            "Variants by domain",
            "Pathogenic fraction by domain"
        ],
        horizontal_spacing=0.12
    )

    # Figure 1.1 - pathogenic variants
    fig.add_trace(
        go.Bar(
            x=x_order,
            y=top_domains["pathogenic"],
            name="pathogenic",
            marker_color="#d62728",
            customdata=top_domains[["non_pathogenic", "total", "pathogenic_fraction"]],
        ),
        row=1,
        col=1
    )

    # Figure 1.2 - non-pathogenic variants
    fig.add_trace(
        go.Bar(
            x=x_order,
            y=top_domains["non_pathogenic"],
            name="non_pathogenic",
            marker_color="#1f77b4",
            customdata=top_domains[["pathogenic", "total", "pathogenic_fraction"]],
        ),
        row=1,
        col=1
    )

    # Figure 2 - pathogenic fraction
    fig.add_trace(
        go.Bar(
            x=x_order,
            y=top_domains["pathogenic_fraction"],
            name="pathogenic fraction",
            marker_color="#4f4f4f",
            showlegend=False,
            customdata=top_domains[["pathogenic", "non_pathogenic", "total"]],
        ),
        row=1,
        col=2
    )

    fig.update_yaxes(title_text="Number of variants", row=1, col=1)
    fig.update_yaxes(title_text="Fraction of pathogenic variants", row=1, col=2)

    fig.update_layout(
        barmode="stack",
        width=1400,
        height=500,
    )

    # Domain names
    fig.update_xaxes(
        tickangle=-45,
        title_text="Domain",
        categoryorder="array",
        categoryarray=x_order,
        row=1,
        col=1
    )
    fig.update_xaxes(
        tickangle=-45,
        title_text="Domain",
        categoryorder="array",
        categoryarray=x_order,
        row=1,
        col=2
    )

    fig.write_image(FIGURES / "1XX_domain_distribution.png", scale=2)

    return domain_counts

# Main function
def run_exploratory() -> None:
    df = load_annotated()
    df = prepare_dataframe(df)

    exon_table = save_exon_distrib_plot(df)
    domain_table = save_domain_distrib_plot(df)

    print(f"Saved figures to: {FIGURES}")
    print(exon_table)

# By launching main, we obtain all the figures in the ../figures
if __name__ == "__main__":
    run_exploratory()






