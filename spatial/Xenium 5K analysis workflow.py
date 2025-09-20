# Xenium Spatial Transcriptomics Full Analysis Pipeline

import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from sklearn.neighbors import NearestNeighbors
from scipy.sparse import issparse
from scipy.ndimage import gaussian_filter1d
from xb.formatting import format_xenium_adata_final
from xb.preprocessing import preprocess_adata, main_preprocessing

# ----------------------------
# Step 1: Load and Format Xenium Data
# ----------------------------

input_path = "xenium/outs/output-XETG00365__0029785__Region_1__20241017__062442"  # Xenium output folder
output_path = "xenium/29785-1"
tag = "29785_1"

adata = format_xenium_adata_final(
    path=input_path,
    tag=tag,
    output_path=output_path,
    use_parquet=True,
    save=True
)

adata.write_h5ad(f"{output_path}/adataformat.h5ad")

# ----------------------------
# Step 2: Preprocessing
# ----------------------------

adata = main_preprocessing(
    adata=adata,
    target_sum=10000,
    mincounts=50,
    mingenes=20,
    neigh=30,
    npc=50,
    scale=False,
    hvg=False,
    total_clusters=20,
    norm=True,
    lg=True
)

clustering_params = {
    "min_counts_x_cell": 50,
    "min_genes_x_cell": 20,
    "normalization_target_sum": 10000,
    "scale": False,
    "hvg": False,
    "norm": True,
    "log1p": True,
    "n_neighbors": 30,
    "n_pcs": 50,
    "clustering_alg": "leiden",
    "resolutions": [1.4],
    "umap_min_dist": 0.6
}

adata = preprocess_adata(
    adata=adata,
    save=True,
    clustering_params=clustering_params,
    output_path=output_path
)

# ----------------------------
# Step 3: Manual Cell Type Annotation
# ----------------------------

cluster_annotation = {
    "0": "T cells", "1": "Vascular endothelial cells", "2": "DPT⁺/CXCL12⁺ CAFs",
    "3": "Muscle cells", "4": "MRC1⁺/CD163⁺ macrophages", "5": "Muscle cells",
    "6": "POSTN⁺/MMP14⁺ CAFs", "7": "Vascular smooth muscle cells",
    "8": "Basal-like malignant cells", "9": "Plasma cells",
    "10": "Basal-like malignant cells", "11": "Plasma cells",
    "12": "Well-differentiated malignant cells", "13": "POSTN⁺/MMP14⁺ CAFs",
    "14": "Cancer-testis-like malignant cells",
    "15": "Low-TP63⁺ pre-malignant squamous cells",
    "16": "POSTN⁺/MMP14⁺ CAFs", "17": "Well-differentiated malignant cells",
    "18": "High-TP63⁺ pre-malignant squamous cells", "19": "Mast cells",
    "20": "Vascular endothelial cells", "21": "Vascular endothelial cells"
}

adata.obs['celltype'] = adata.obs['leiden_1.4'].map(cluster_annotation)

# Color map
celltype_colors = {
    "T cells": "#E05F48", "Plasma cells": "#A5AA99",
    "POSTN⁺/MMP14⁺ CAFs": "#C7AED5", "MRC1⁺/CD163⁺ macrophages": "#ACD48A",
    "Basal-like malignant cells": "#F7DBF0", "Vascular endothelial cells": "#F5BC6E",
    "Well-differentiated malignant cells": "#EA945A",
    "Muscle cells": "#726BAE", "Cancer-testis-like malignant cells": "#86A667",
    "DPT⁺/CXCL12⁺ CAFs": "#276D9F", "Mast cells": "#8AB6D6",
    "Vascular smooth muscle cells": "#E5948E",
    "Low-TP63⁺ pre-malignant squamous cells": "#5D88BF",
    "High-TP63⁺ pre-malignant squamous cells": "#78C4D4"
}

adata.uns['celltype_colors'] = [
    celltype_colors.get(ct, "#000000")
    for ct in adata.obs['celltype'].astype('category').cat.categories
]

# ----------------------------
# Step 4: Spatial Celltype Plot
# ----------------------------

fig, ax = plt.subplots(figsize=(10, 10))
sc.pl.spatial(
    adata, color="celltype", spot_size=25,
    ax=ax, show=False, frameon=False, legend_loc=None
)

x0, x1 = ax.get_xlim()
y0, y1 = ax.get_ylim()
rect = Rectangle((x0, y0), x1 - x0, y1 - y0, fill=False, edgecolor='black', linewidth=2)
ax.add_patch(rect)

plt.tight_layout()
plt.savefig("spatial_celltype_bordered_highres.pdf", dpi=3000, bbox_inches='tight')
plt.close()

# ----------------------------
# Step 5: Hotspot Scoring
# ----------------------------

def compute_scores_hotspot(adata, genes, n_neighbors=30, neighborhood_factor=3, use_rep="spatial"):
    genes = [gene for gene in genes if gene in adata.var_names]
    coordinates = adata.obsm[use_rep]
    nbrs = NearestNeighbors(n_neighbors=n_neighbors).fit(coordinates)
    distances, indices = nbrs.kneighbors(coordinates)
    weights = np.exp(-distances / neighborhood_factor)
    weights /= weights.sum(axis=1, keepdims=True)
    counts_dense = adata[:, genes].X.toarray() if issparse(adata[:, genes].X) else adata[:, genes].X
    scores = [np.average(counts_dense[indices[i]], weights=weights[i], axis=0).mean() for i in range(len(coordinates))]
    return np.array(scores)

CD4T_ISG15_Signature = ["ISG15", "IFIT3", "IFIT1", "OASL", "MX1", "IFI44L", "IFI6", "CD4", "IFIT2",
    "EPSTI1", "CD3D", "RSAD2", "CD3E", "LTB", "OAS1", "STAT1", "LY6E", "OAS3"]

adata.obs["CD4T_ISG15_Score"] = compute_scores_hotspot(adata, CD4T_ISG15_Signature)

# ----------------------------
# Step 6: Hotspot Spatial Plot (T cells only)
# ----------------------------

adata.obs["CD4T_ISG15_Score_Tcells_only"] = adata.obs.apply(
    lambda row: row["CD4T_ISG15_Score"] if row["celltype"] == "T cells" else np.nan,
    axis=1
)

sc.settings.figdir = "."
sc.pl.spatial(
    adata,
    color="CD4T_ISG15_Score_Tcells_only",
    spot_size=20,
    vmax=1,
    cmap="viridis",
    na_color="lightgrey",
    show=False
)

fig = plt.gcf()
fig.patch.set_facecolor('black')
for ax in fig.axes:
    ax.set_facecolor('black')
    ax.tick_params(colors='white')
    for spine in ax.spines.values():
        spine.set_edgecolor('white')

plt.tight_layout()
plt.savefig("CD4T_ISG15_Tcells_only_blackbg.png", dpi=600, bbox_inches="tight", facecolor='black')
plt.close()

# ----------------------------
# Step 7: Hotspot Vertical Density Plot
# ----------------------------

obs = adata.obs.copy()
obs["y"] = adata.obsm["spatial"][:, 1]
tcell_obs = obs[~obs["CD4T_ISG15_Score_Tcells_only"].isna()]

n_bins = 200
y_min, y_max = tcell_obs["y"].min(), tcell_obs["y"].max()
bin_edges = np.linspace(y_min, y_max, n_bins + 1)
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

scores = tcell_obs["CD4T_ISG15_Score_Tcells_only"].values
y_coords = tcell_obs["y"].values
score_sum_per_bin = np.zeros(n_bins)
bin_indices = np.digitize(y_coords, bins=bin_edges) - 1
for i, score in zip(bin_indices, scores):
    if 0 <= i < n_bins:
        score_sum_per_bin[i] += score

smoothed_scores = gaussian_filter1d(score_sum_per_bin, sigma=2)

plt.figure(figsize=(4, 10))
plt.plot(smoothed_scores, bin_centers, color="lime", linewidth=2.5)
plt.fill_betweenx(bin_centers, smoothed_scores, color="lime", alpha=0.2)
plt.ylabel("Y Position (µm)", fontsize=12)
plt.xlabel("CD4T_ISG15 Score Density", fontsize=12)
plt.title("CD4T_ISG15 Vertical Density", fontsize=14)
plt.gca().invert_yaxis()
plt.tight_layout()
plt.savefig("CD4T_ISG15_VerticalDensityProfile_flipped.png", format='png', dpi=1200, bbox_inches='tight')
plt.show()

# Save final annotated object
adata.write_h5ad(f"{output_path}/adataanno.h5ad")
