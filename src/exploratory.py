from __future__ import annotations

# Some preliminary imports
import pandas as pd
from pathlib import Path
from scipy.stats import chi2_contingency, fisher_exact
import matplotlib.pyplot as plt

# Constant paths
PROCESSED_DATA = Path("../data/processed")
FIGURES = Path("../figures")
EDA_TABLES = PROCESSED_DATA / "eda_tables"

# Loads the annotated CSV
def load_annotated(path: Path | str = PROCESSED_DATA / "DMD_variants_annotated.csv") -> pd.DataFrame:
    df = pd.read_csv(path)
    return df