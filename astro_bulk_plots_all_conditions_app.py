import streamlit as st
import pandas as pd
import plotly.express as px
import uuid
import os
import plotly.graph_objects as go

st.set_page_config(layout="wide")
st.title("🔬 Bulk RNA-seq Gene Expression Explorer")

# === Default file paths ===
default_expr_path = "/Users/zeenalpatel/Astrocyte_Integration_SCI_Project/HPC-astro_integration_project/AnalysisRDSfiles/CSV_GSE241628/astro_bulk_with_symbols.xlsx"
default_meta_path = "/Users/zeenalpatel/Astrocyte_Integration_SCI_Project/HPC-astro_integration_project/AnalysisRDSfiles/CSV_GSE241628/astro_bulk_metadata.csv"

# === Load expression matrix from default path ===
if os.path.exists(default_expr_path):
    expr_df = pd.read_excel(default_expr_path)
    st.success("✅ Loaded default expression matrix.")
else:
    st.error("❌ Default expression matrix not found.")
    st.stop()

# === Load metadata from default path ===
if os.path.exists(default_meta_path):
    metadata = pd.read_csv(default_meta_path, index_col=0)
    st.success("✅ Loaded default metadata.")
else:
    st.error("❌ Default metadata file not found.")
    st.stop()

# === Prepare and match sample IDs ===
expr_df.columns = expr_df.columns.astype(str)
expr_df = expr_df.rename(columns={expr_df.columns[0]: "GeneID", expr_df.columns[1]: "GeneSymbol"})
metadata.index = metadata.index.astype(str)

shared_samples = list(set(expr_df.columns[2:]) & set(metadata.index))
if not shared_samples:
    st.error("❌ No overlapping sample IDs between expression matrix and metadata.")
    st.stop()

expr_df = expr_df[["GeneID", "GeneSymbol"] + shared_samples]
metadata = metadata.loc[shared_samples]

# === Define custom orders for metadata columns ===
ordered_categories = {
    "Treatment": ["Postnatal", "Uninjured", "SCI"],
    "Time": ["P0", "P3", "P5", "P7", "P14", "P21", "P35", "P63", 
             "Uninjured", "2d", "5d", "14d", "28d", "42d", "70d"],
    "CellType": ["Astrocytes", "Flow"],
    "Genotype": ["73.12 GFAP Cre -RiboTag", "Aldh1l1-CreERT2 -RiboTag"],
    "Sex": ["F", "M"],
    "Replicate": [1, 2, 3, 4, 5]
}
# Convert to categorical dtype with the specified order
for col, order in ordered_categories.items():
    if col in metadata.columns:
        metadata[col] = pd.Categorical(metadata[col], categories=order, ordered=True)

# === Dynamic Filters: Create a dropdown for each metadata column ===
st.markdown("### 🎯 Filter samples by metadata (cascading, side-by-side)")

# Columns for side-by-side layout
col_objs = st.columns(6)

# Track filtered metadata and sample list
filtered_metadata = metadata.copy()
filtered_samples = metadata.index.tolist()
filter_values = {}

# List of fields in desired order
filter_order = ["Treatment", "Time", "CellType", "Genotype", "Sex", "Replicate"]

# Loop through each column and add multiselect inside the column
for i, col in enumerate(filter_order):
    col_data = filtered_metadata[col].dropna()
    # Use category order if set
    if pd.api.types.is_categorical_dtype(metadata[col]):
        options = [cat for cat in metadata[col].cat.categories if cat in col_data.values]
    else:
        options = sorted(col_data.unique().tolist())
    with col_objs[i]:
        selected = st.multiselect(f"{col}:", options, default=options, key=col)
        filter_values[col] = selected
    # Apply filtering if not selecting all
    if selected and set(selected) != set(options):
        filtered_metadata = filtered_metadata[filtered_metadata[col].isin(selected)]
        filtered_samples = filtered_metadata.index.tolist()

# Final filtered sample list
if not filtered_samples:
    st.warning("⚠️ No samples match the selected filters.")
    st.stop()

# Filter expression matrix to match
expr_df = expr_df[["GeneID", "GeneSymbol"] + filtered_samples]

# Sort by multiple ordered columns (if they exist)
sort_cols = [col for col in ordered_categories if col in filtered_metadata.columns]
filtered_metadata = filtered_metadata.sort_values(by=sort_cols)

# Update filtered_samples to match sorted metadata index
filtered_samples = filtered_metadata.index.tolist()

# Reorder expression matrix columns accordingly
expr_df = expr_df[["GeneID", "GeneSymbol"] + filtered_samples]

# === Checkpoint: Show filtered results ===
#st.markdown("### Filtered Sample Information")
#st.write(f"🧪 Number of samples after filtering: {len(filtered_samples)}")
#st.dataframe(filtered_metadata)

st.markdown("### Filtered Expression Matrix (preview)")
##st.dataframe(expr_df.set_index("GeneSymbol"))

# === 🧪 Summary: Mean ± SD Across Replicates ===
#st.markdown("## 🧮 Summary Table: Mean ± SD by Condition")

# 1. Melt expression matrix
df_long = expr_df.melt(id_vars=["GeneID", "GeneSymbol"], 
                       var_name="Sample", 
                       value_name="Expression")

# 2. Merge metadata
df_long = df_long.merge(filtered_metadata, left_on="Sample", right_index=True)

# 3. Define composite condition group
group_vars = ["Treatment", "Time", "CellType", "Genotype", "Sex"]
df_long["Group"] = df_long[group_vars].astype(str).agg("_".join, axis=1)

# 4. Group and compute mean ± SD
summary_df = df_long.groupby(["GeneSymbol", "Group"]).agg(
    Mean=("Expression", "mean"),
    SD=("Expression", "std"),
    Replicates=("Expression", "count")
).reset_index()

# 5. Display
#st.dataframe(summary_df)

# === Gene search with metadata ===
st.markdown("### 🔍 Search for a gene")

search_genes = st.text_area(
    "Enter one or more Gene Symbols or IDs (comma-separated or one per line):"
)
gene_list = [g.strip() for g in search_genes.replace(",", "\n").splitlines() if g.strip()]

for search_gene in gene_list:
    # Lookup gene symbol and ID
    lookup_row = expr_df[
        (expr_df["GeneSymbol"] == search_gene) | (expr_df["GeneID"] == search_gene)
    ]

    if lookup_row.empty:
        st.warning(f"⚠️ No matching gene found for `{search_gene}`.")
        continue

    gene_symbol = lookup_row["GeneSymbol"].values[0]
    gene_id = lookup_row["GeneID"].values[0]
    st.markdown(f"**Gene Symbol**: `{gene_symbol}`  &nbsp;&nbsp; | &nbsp;&nbsp; **Gene ID**: `{gene_id}`")

if search_gene:
    # Filter expression matrix for that gene
    matched_rows = expr_df[
        (expr_df["GeneSymbol"] == search_gene) | (expr_df["GeneID"] == search_gene)
    ]
    
    if matched_rows.empty:
        st.warning("⚠️ No matching gene found.")
    else:
        # Transpose expression values (samples as rows)
        expr_t = matched_rows.set_index("GeneSymbol").T
        expr_t.columns = ["Expression"]
        expr_t = expr_t.drop(index="GeneID", errors="ignore")  # just in case

        # Merge with metadata to show sample info
        expr_with_meta = expr_t.merge(metadata, left_index=True, right_index=True)

        expr_with_meta = expr_with_meta.sort_values(by=sort_cols)

        #st.markdown("### 📊 Expression of selected gene with sample metadata:")
        #st.dataframe(expr_with_meta)

# === Mean ± SD just for searched gene
# === Mean ± SD for searched gene — stratified by Sex + pooled
pooled_group_vars = ["Treatment", "Time", "CellType", "Genotype"]
sex_group_vars = pooled_group_vars + ["Sex"]

# 1. Stratified by Sex
sex_summary = (
    expr_with_meta
    .groupby(sex_group_vars)
    .agg(
        Mean=("Expression", "mean"),
        SD=("Expression", "std"),
        Replicates=("Expression", "count")
    )
    .reset_index()
)
sex_summary["SexGroup"] = sex_summary["Sex"]

# 2. Pooled across Sex (corrected: compute across all replicates)
expr_with_meta["Sex"] = "Both"  # overwrite Sex column
pooled_summary = (
    expr_with_meta
    .groupby(pooled_group_vars)
    .agg(
        Mean=("Expression", "mean"),
        SD=("Expression", "std"),
        Replicates=("Expression", "count")
    )
    .reset_index()
)
pooled_summary["Sex"] = "Both"
pooled_summary["SexGroup"] = "Both"

# 3. Combine
search_summary = pd.concat([sex_summary, pooled_summary], ignore_index=True)

# Remove groups with zero replicates or all-NaN values
search_summary = search_summary[search_summary["Replicates"] > 0]
search_summary = search_summary[search_summary["Mean"].notna()]

if gene_list:
    n_cols = 2
    rows = [gene_list[i:i + n_cols] for i in range(0, len(gene_list), n_cols)]

    for row_genes in rows:
        cols = st.columns(n_cols)
        for i, search_gene in enumerate(row_genes):
            with cols[i]:
                matched_rows = expr_df[
                    (expr_df["GeneSymbol"] == search_gene) | (expr_df["GeneID"] == search_gene)
                ]
                if matched_rows.empty:
                    st.warning(f"⚠️ No matching gene found for `{search_gene}`.")
                    continue

                expr_t = matched_rows.set_index("GeneSymbol").T
                expr_t.columns = ["Expression"]
                expr_t = expr_t.drop(index="GeneID", errors="ignore")

                expr_with_meta = expr_t.merge(metadata, left_index=True, right_index=True)
                expr_with_meta = expr_with_meta.sort_values(by=sort_cols)

                pooled_group_vars = ["Treatment", "Time", "CellType", "Genotype"]
                sex_group_vars = pooled_group_vars + ["Sex"]

                sex_summary = (
                    expr_with_meta
                    .groupby(sex_group_vars)
                    .agg(
                        Mean=("Expression", "mean"),
                        SD=("Expression", "std"),
                        Replicates=("Expression", "count")
                    )
                    .reset_index()
                )
                sex_summary["SexGroup"] = sex_summary["Sex"]

                expr_with_meta["Sex"] = "Both"
                pooled_summary = (
                    expr_with_meta
                    .groupby(pooled_group_vars)
                    .agg(
                        Mean=("Expression", "mean"),
                        SD=("Expression", "std"),
                        Replicates=("Expression", "count")
                    )
                    .reset_index()
                )
                pooled_summary["Sex"] = "Both"
                pooled_summary["SexGroup"] = "Both"

                search_summary = pd.concat([sex_summary, pooled_summary], ignore_index=True)
                search_summary = search_summary[search_summary["Replicates"] > 0]
                search_summary = search_summary[search_summary["Mean"].notna()]

                search_summary["Group"] = (
                    search_summary["SexGroup"].astype(str) + " | " +
                    search_summary["Treatment"].astype(str) + " | " +
                    search_summary["Time"].astype(str) + " | " +
                    search_summary["Genotype"].astype(str)
                )
                search_summary["LineGroup"] = (
                    search_summary["Treatment"].astype(str) + " | " +
                    search_summary["CellType"].astype(str) + " | " +
                    search_summary["SexGroup"].astype(str)
                )

                color_discrete_map = {
                    "Postnatal | Astrocytes | M": "#1b55a0",
                    "Postnatal | Flow | M": "#3EA3B5",

                    "Postnatal | Astrocytes | F": "#9e4f9d",
                    "Postnatal | Flow | F": "#e24fb1",

                    "Uninjured | Astrocytes | F": "#e174b4",
                    "Uninjured | Flow | F": "#e24da4",

                    "SCI | Astrocytes | M": "#1b78a0",
                    "SCI | Flow | M": "#3EA3B5",

                    "SCI | Astrocytes | F": "#a02e9e",
                    "SCI | Flow | F": "#e24fb1",

                    "Postnatal | Astrocytes | Both": "#474101",
                    "Postnatal | Flow | Both": "#AE8944",

                    "Uninjured | Astrocytes | Both": "#000000",
                    "Uninjured | Flow | Both": "#B061BC",

                    "SCI | Astrocytes | Both": "#000000",
                    "SCI | Flow | Both": "#A954B6"
                }

                fig = px.line(
                    search_summary,
                    x="Group",
                    y="Mean",
                    error_y="SD",
                    color="LineGroup",
                    line_group="CellType",
                    markers=True,
                    symbol="SexGroup",
                    title=f"Mean ± SD Expression of {search_gene}",
                    labels={"Mean": "Mean Expression", "Group": "Condition"},
                    height=500,
                    color_discrete_map=color_discrete_map
                )

                fig.update_layout(xaxis_tickangle=-45)
                tickvals = search_summary["Group"].tolist()
                ticktext = search_summary["Time"].astype(str).tolist()
                fig.update_xaxes(tickvals=tickvals, ticktext=ticktext, tickangle=-45)

                st.plotly_chart(fig, use_container_width=True)

