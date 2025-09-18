# =========================
# Pan-cancer integration pipeline (Scanpy + BBKNN)
# =========================

import os
import numpy as np
import pandas as pd
import scanpy as sc
import bbknn
import matplotlib.pyplot as plt  # optional

# ---- Input/Output paths ----
in_h5ad  = "/home/guile/CvI/pancancer/adatapancer.h5ad"       # raw merged AnnData
out_h5ad = "/home/guile/CvI/pancancer/adataumap_leiden15.h5ad" # processed output

# ---- Load data ----
adata = sc.read_h5ad(in_h5ad)

# ---- Build batch key for BBKNN (Sample_Dataset) ----
adata.obs["batch"] = adata.obs["Sample"].astype(str) + "_" + adata.obs["Dataset"].astype(str)

# ---- Quality control at dataset level ----
# Keep only datasets with at least 300 T cells
# (assumes 'celltype' column exists and T cells labeled as "T cell")
t_counts = adata.obs.query("celltype == 'T cell'").groupby("Dataset").size()
keep_datasets = set(t_counts[t_counts >= 300].index)
adata = adata[adata.obs["Dataset"].isin(keep_datasets)].copy()

# ---- Normalize and log-transform ----
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# ---- Highly variable genes ----
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var["highly_variable"]].copy()

# ---- Batch correction with BBKNN ----
bbknn.bbknn(adata, batch_key="batch")

# ---- UMAP embedding ----
sc.tl.umap(adata)

# ---- Leiden clustering ----
sc.tl.leiden(adata, resolution=1.5, key_added="leiden")

# ---- Save processed object ----
adata.write_h5ad(out_h5ad)

