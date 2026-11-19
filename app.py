"""Streamlit version of the GeneScope co-expression heatmap."""

from __future__ import annotations

import io

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import streamlit as st


SAMPLE_DATA = """Gene,Sample_1,Sample_2,Sample_3,Sample_4,Sample_5
BRCA1,12.5,15.2,14.0,20.1,22.4
BRCA2,11.8,14.9,13.8,19.5,21.9
TP53,80.1,65.4,50.2,30.1,15.0
EGFR,45.2,46.1,44.8,45.9,45.0
MYC,10.2,12.1,25.4,30.8,40.2
"""


st.set_page_config(
    page_title="GeneScope | Co-expression lab",
    page_icon="🧬",
    layout="wide",
)

st.markdown(
    """
    <style>
    .stApp { background: #f5f3ef; }
    [data-testid="stMetric"] {
        background: #fffefa;
        border: 1px solid #e5e5df;
        border-radius: 12px;
        padding: 14px 16px;
    }
    div[data-testid="stAlert"] { border-radius: 10px; }
    </style>
    """,
    unsafe_allow_html=True,
)


def parse_expression_data(raw_text: str) -> pd.DataFrame:
    """Read a gene-by-sample CSV or TSV table and validate its values."""
    cleaned = raw_text.lstrip("\ufeff").strip()
    if not cleaned:
        raise ValueError("Add a table or upload a CSV/TSV file to begin.")

    first_line = cleaned.splitlines()[0]
    separator = "\t" if first_line.count("\t") > first_line.count(",") else ","

    try:
        table = pd.read_csv(io.StringIO(cleaned), sep=separator)
    except Exception as error:
        raise ValueError(f"I couldn't read that table: {error}") from error

    if table.shape[1] < 3:
        raise ValueError("Include one gene-name column and at least two sample columns.")
    if table.shape[0] < 2:
        raise ValueError("Include at least two genes to compare.")

    headers = [str(column).strip() for column in table.columns]
    if not headers[0] or headers[0].lower().startswith("unnamed:"):
        raise ValueError("The first column must have a name such as Gene.")
    if any(not name or name.lower().startswith("unnamed:") for name in headers[1:]):
        raise ValueError("Every sample column needs a header.")
    table.columns = headers

    gene_names = table.iloc[:, 0].astype("string").str.strip()
    if gene_names.isna().any() or gene_names.eq("").any():
        raise ValueError("Every row needs a gene name in the first column.")
    if gene_names.str.casefold().duplicated().any():
        raise ValueError("Gene names must be unique.")

    expression = table.iloc[:, 1:].copy()
    expression.columns = headers[1:]
    expression.index = gene_names.astype(str)
    expression = expression.apply(pd.to_numeric, errors="coerce")

    if expression.isna().any().any():
        raise ValueError("Every sample value must be a number; blank or text values were found.")
    if not np.isfinite(expression.to_numpy(dtype=float)).all():
        raise ValueError("Sample values must be finite numbers.")

    return expression.astype(float)


def read_uploaded_file(uploaded_file) -> str:
    try:
        return uploaded_file.getvalue().decode("utf-8-sig")
    except UnicodeDecodeError as error:
        raise ValueError("The file must be saved as UTF-8 text.") from error


def find_gene_groups(correlations: pd.DataFrame, cutoff: float) -> list[list[str]]:
    """Group genes connected by positive correlations at or above cutoff."""
    genes = list(correlations.index)
    neighbors = {gene: set() for gene in genes}

    for left_index, left_gene in enumerate(genes):
        for right_index in range(left_index + 1, len(genes)):
            right_gene = genes[right_index]
            value = correlations.iat[left_index, right_index]
            if pd.notna(value) and value >= cutoff:
                neighbors[left_gene].add(right_gene)
                neighbors[right_gene].add(left_gene)

    visited = set()
    groups = []
    for gene in genes:
        if gene in visited:
            continue
        stack = [gene]
        visited.add(gene)
        group = []
        while stack:
            current = stack.pop()
            group.append(current)
            for neighbor in neighbors[current]:
                if neighbor not in visited:
                    visited.add(neighbor)
                    stack.append(neighbor)
        if len(group) > 1:
            groups.append(sorted(group))

    return groups


st.title("GeneScope")
st.caption("Turn gene expression data into a vivid map of genes moving together.")
st.write(
    "Compare genes across samples with Pearson correlation, find similar patterns, "
    "and download the correlation matrix."
)

with st.sidebar:
    st.header("Your data")
    input_mode = st.radio(
        "Choose an input method",
        ["Example data", "Upload CSV/TSV", "Paste data"],
        index=0,
    )

    raw_data = SAMPLE_DATA
    if input_mode == "Upload CSV/TSV":
        uploaded_file = st.file_uploader("Choose a data file", type=["csv", "tsv"])
        if uploaded_file is None:
            st.info("Upload a CSV or TSV file to continue.")
            st.stop()
        try:
            raw_data = read_uploaded_file(uploaded_file)
        except ValueError as error:
            st.error(str(error))
            st.stop()
    elif input_mode == "Paste data":
        raw_data = st.text_area(
            "Paste CSV or tab-separated data",
            value=SAMPLE_DATA,
            height=220,
            help="First column: gene names. Remaining columns: numeric sample values.",
        )

    st.divider()
    cutoff = st.slider(
        "Positive correlation cutoff for groups",
        min_value=0.70,
        max_value=0.95,
        value=0.85,
        step=0.05,
        format="%.2f",
    )
    st.caption("Groups connect genes with Pearson r at or above the cutoff.")
    st.info(
        "This Streamlit version processes uploaded data on the hosting server. "
        "Use only data you are allowed to upload."
    )


try:
    expression = parse_expression_data(raw_data)
except ValueError as error:
    st.error(str(error))
    st.stop()

correlations = expression.T.corr(method="pearson")
genes = list(expression.index)
sample_count = expression.shape[1]

pairs = []
for left_index, left_gene in enumerate(genes):
    for right_index in range(left_index + 1, len(genes)):
        right_gene = genes[right_index]
        value = correlations.iat[left_index, right_index]
        if pd.notna(value):
            pairs.append((left_gene, right_gene, float(value)))

st.subheader("Dataset summary")
metric_columns = st.columns(4)
metric_columns[0].metric("Genes", len(genes))
metric_columns[1].metric("Samples", sample_count)

if pairs:
    strongest = max(pairs, key=lambda pair: pair[2])
    most_inverse = min(pairs, key=lambda pair: pair[2])
    metric_columns[2].metric(
        "Strongest positive",
        f"{strongest[2]:+.2f}",
        f"{strongest[0]} · {strongest[1]}",
    )
    metric_columns[3].metric(
        "Most inverse",
        f"{most_inverse[2]:+.2f}",
        f"{most_inverse[0]} · {most_inverse[1]}",
    )
else:
    metric_columns[2].metric("Strongest positive", "—", "No variable gene pair")
    metric_columns[3].metric("Most inverse", "—", "No variable gene pair")

st.subheader("Gene-to-gene Pearson correlation")
st.caption(
    f"{len(genes)} genes × {sample_count} samples · each cell shows the correlation across samples."
)

side = max(5.5, min(15.0, len(genes) * 0.72))
figure, axis = plt.subplots(figsize=(side + 1.2, side))
sns.heatmap(
    correlations,
    annot=True,
    fmt=".4f",
    cmap="coolwarm",
    vmin=-1,
    vmax=1,
    center=0,
    square=True,
    linewidths=0.6,
    linecolor="white",
    mask=correlations.isna(),
    cbar_kws={"label": "Pearson r"},
    ax=axis,
)
axis.set_title("Co-expression heatmap", pad=14, weight="bold")
axis.set_xlabel("Gene")
axis.set_ylabel("Gene")
axis.tick_params(axis="x", rotation=45)
axis.tick_params(axis="y", rotation=0)
figure.tight_layout()
st.pyplot(figure, use_container_width=True)
plt.close(figure)

st.caption("Blue indicates negative correlation, pale colors are near zero, and red indicates positive correlation.")

st.subheader("Similar-gene groups")
groups = find_gene_groups(correlations, cutoff)
grouped_genes = {gene for group in groups for gene in group}

if groups:
    for number, group in enumerate(groups, start=1):
        st.text(f"Group {number} · {len(group)} genes: {', '.join(group)}")
    ungrouped = [gene for gene in genes if gene not in grouped_genes]
    if ungrouped:
        with st.expander("Genes below this cutoff"):
            st.text(", ".join(ungrouped))
else:
    st.info(f"No gene pairs meet r ≥ {cutoff:.2f}. Try a lower cutoff.")

st.subheader("Download results")
csv_bytes = correlations.to_csv(index=True, index_label="Gene").encode("utf-8")
st.download_button(
    "Download correlation matrix CSV",
    data=csv_bytes,
    file_name="gene-correlation-matrix.csv",
    mime="text/csv",
)

st.caption(
    "Pearson correlation describes a linear pattern across samples; it does not prove that one gene causes a change in another."
)
