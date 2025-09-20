import os
import scanpy as sc
import anndata as ad

# -----------------------------
# Step 1: Load Pre-annotated Samples
# -----------------------------
adata_files = {
    "136341T": "PATH_TO/136341T.h5ad",
    "136992T": "PATH_TO/136992T.h5ad",
    "137056T": "PATH_TO/137056T.h5ad",
    "139540D": "PATH_TO/139540D.h5ad",
    "139907A": "PATH_TO/139907A.h5ad",
    "137486A": "PATH_TO/137486A.h5ad"
}

adatas = []
for sample_id, path in adata_files.items():
    adata = sc.read_h5ad(path)
    adata.obs["sample"] = sample_id
    adatas.append(adata)

# -----------------------------
# Step 2: Concatenate into One AnnData
# -----------------------------
adata = ad.concat(adatas, label="sample", keys=list(adata_files.keys()))

# -----------------------------
# Step 3: Preprocessing
# -----------------------------
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.scale(adata, max_value=10)
sc.pp.pca(adata)

# -----------------------------
# Step 4: Harmony Batch Correction
# -----------------------------
import scarches as sce  # ensure harmony is available via `scarches`
sce.pp.harmony_integrate(adata, key="sample")

# -----------------------------
# Step 5: Embedding and Neighbors
# -----------------------------
sc.pp.neighbors(adata, n_neighbors=30, n_pcs=50, use_rep='X_pca_harmony')
sc.tl.umap(adata)
adata.write_h5ad("merged_harmony_umap.h5ad")
