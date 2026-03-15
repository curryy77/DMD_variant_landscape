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
            "Pathogenic fraction by domain",
            "Variants by domain"
        ],
        horizontal_spacing=0.12
    )

    # Figure 1 - pathogenic fraction
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
        col=1
    )

    # Figure 2.1 - pathogenic variants
    fig.add_trace(
        go.Bar(
            x=x_order,
            y=top_domains["pathogenic"],
            name="pathogenic",
            marker_color="#d62728",
            customdata=top_domains[["non_pathogenic", "total", "pathogenic_fraction"]],
        ),
        row=1,
        col=2
    )

    # Figure 2.2 - non-pathogenic variants
    fig.add_trace(
        go.Bar(
            x=x_order,
            y=top_domains["non_pathogenic"],
            name="non_pathogenic",
            marker_color="#1f77b4",
            customdata=top_domains[["pathogenic", "total", "pathogenic_fraction"]],
        ),
        row=1,
        col=2
    )

    fig.update_yaxes(title_text="Fraction of pathogenic variants", row=1, col=1)
    fig.update_yaxes(title_text="Number of variants", row=1, col=2)


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

# Main plot #3: Histograms of frame status vs. phenotype class
def save_reading_frame_plot(df: pd.DataFrame) -> pd.DataFrame:
    tmp = df.dropna(subset=["frame_status", "phenotype_class"]).copy()

    # Contingency table
    frame_counts = (
        tmp.groupby(["frame_status", "phenotype_class"])
        .size()
        .unstack(fill_value=0)
        .reset_index()
    )

    frame_counts["total"] = frame_counts.iloc[:,1:].sum(axis=1)

    # Counting proportions
    for col in frame_counts.columns:
        if col not in ["frame_status", "total"]:
            frame_counts[col + "_prop"] = frame_counts[col] / frame_counts["total"]

    phenotypes = [c for c in frame_counts.columns if c not in ["frame_status","total"] and not c.endswith("_prop")]
    phenotypes_prop = [p + "_prop" for p in phenotypes]

    x = frame_counts["frame_status"]

    fig = make_subplots(
        rows=1,
        cols=2,
        subplot_titles=[
            "Reading frame rule: counts",
            "Reading frame rule: proportions"
        ],
        horizontal_spacing=0.12
    )

    colors = {
        "DMD": "#e74c3c",
        "BMD": "#4c6ef5",
        "other": "#afafaf"
    }


    # Figure 1 - Reading frame rule: counts
    for p in phenotypes:
        fig.add_trace(
            go.Bar(
                x=x,
                y=frame_counts[p],
                name=p,
                marker_color=colors.get(p, None)
            ),
            row=1,
            col=1
        )

    # Figure 2 - Reading frame rule: proportions
    for p, pp in zip(phenotypes, phenotypes_prop):
        fig.add_trace(
            go.Bar(
                x=x,
                y=frame_counts[pp],
                name=p,
                marker_color=colors.get(p, None),
                showlegend=False
            ),
            row=1,
            col=2
        )

    fig.update_layout(
        barmode="stack",
        width=1400,
        height=600,
        title="Reading frame rule: frame status vs phenotype"
    )

    fig.update_yaxes(title="Number of variants", row=1, col=1)
    fig.update_yaxes(title="Proportion of variants", row=1, col=2)

    fig.update_xaxes(title="Frame status")

    fig.write_image(FIGURES / "2XX_reading_frame_rule.png", scale=2)

    return frame_counts

# Main function
def run_exploratory() -> None:
    df = load_annotated()
    df = prepare_dataframe(df)

    exon_table = save_exon_distrib_plot(df)
    domain_table = save_domain_distrib_plot(df)
    reading_frame_table = save_reading_frame_plot(df)

    print(f"Saved figures to: {FIGURES}")
    print(exon_table)

# By launching main, we obtain all the figures in the ../figures
if __name__ == "__main__":
    run_exploratory()






