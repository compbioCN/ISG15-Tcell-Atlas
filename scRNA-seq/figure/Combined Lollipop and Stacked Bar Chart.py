# ========================
# Cell-type composition across cancer types
# - Stacked barplot: All T-cell subtypes
# - Lollipop plot: CD4+ Tisg only
# ========================

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scanpy as sc

# Ensure required columns exist
if "Cancer_Type" not in adata.obs.columns or "cell_type" not in adata.obs.columns:
    raise ValueError("'Cancer_Type' or 'cell_type' column not found in adata.obs")

# ==== Stacked barplot: cell-type composition ====

# Count cells per (Cancer_Type × cell_type) and convert to percentage
df_counts = adata.obs.groupby(["Cancer_Type", "cell_type"]).size().unstack(fill_value=0)
df_counts = df_counts.div(df_counts.sum(axis=1), axis=0) * 100  # row-wise %

# Define colors for selected cell types
colors = {
    "CD4+ Tcm": "#C7AED5",
    "CD4+ Teff": "#ACD48A",
    "CD4+ Tem": "#726BAE",
    "CD4+ Treg": "#EA945A",
    "CD4+ Tisg": "#E05F48",
    "CD4+ Tn": "#F7DBF0",
    "Th1-like T cells": "#a9dce6",
    "CD4+ Trm": "#6488B9",
    "CD4+ Tstr": "#8AB6D6"
}

# Match colors to cell types present in df_counts
color_list = [colors[cell] for cell in df_counts.columns if cell in colors]
if not color_list:
    raise ValueError("Color list is empty! Check df_counts.columns and colors dictionary.")

# Plot
fig, ax = plt.subplots(figsize=(10, 5))
df_counts.plot(kind="bar", stacked=True, ax=ax, color=color_list, alpha=0.85)

plt.xlabel("Cancer Type")
plt.ylabel("Percentage (%)")
plt.title("Proportion of Cell Types in Each Cancer Type")
plt.xticks(rotation=45, ha="right")
plt.legend(title="Cell Type", bbox_to_anchor=(1.05, 1), loc="upper left")
plt.tight_layout()
plt.savefig("cancer_type_stacked_bar.png", dpi=300)
plt.savefig("cancer_type_stacked_bar.pdf")
plt.show()

# ==== Lollipop plot: CD4+ Tisg only ====

# Count total cells and CD4+ Tisg cells per cancer type
total_cells = adata.obs["Cancer_Type"].value_counts()
tisg_cells = adata.obs[adata.obs["cell_type"] == "CD4+ Tisg"]["Cancer_Type"].value_counts()

# Calculate percentage
tisg_pct = (tisg_cells / total_cells) * 100
df_tisg = pd.DataFrame({
    "Cancer_Type": tisg_pct.index,
    "Percentage": tisg_pct.values
}).set_index("Cancer_Type")

# Reorder to match df_counts order (optional for visual consistency)
df_tisg = df_tisg.reindex(df_counts.index).dropna().reset_index()

# Plot
fig, ax = plt.subplots(figsize=(10, 5))
ax.stem(
    df_tisg["Cancer_Type"],
    df_tisg["Percentage"],
    linefmt="#726BAE",     # Line color
    markerfmt="#726BAE",   # Dot color
    basefmt=" "            # No baseline
)

plt.xlabel("Cancer Type")
plt.ylabel("Percentage of CD4+ Tisg (%)")
plt.title("CD4+ Tisg Proportion per Cancer Type (Lollipop Plot)")
plt.xticks(rotation=45, ha="right")
plt.tight_layout()
plt.savefig("CD4Tisg_lollipop.png", dpi=300)
plt.savefig("CD4Tisg_lollipop.pdf")
plt.show()
