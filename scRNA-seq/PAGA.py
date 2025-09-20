import scanpy as sc
import anndata
from scipy import io
import numpy as np
import pandas as pd


def seurat_to_adata(counts_path,     # path to counts.mtx (R exported)
                    meta_path,       # path to metadata.csv
                    gene_name_path,  # path to gene_names.csv
                    pca_path,        # path to pca.csv
                    obsm_key,        # e.g. 'X_umap' or 'X_tsne'
                    reduction1,      # e.g. 'UMAP_1' or 'TSNE_1'
                    reduction2):     # e.g. 'UMAP_2' or 'TSNE_2'
    """
    Convert Seurat-exported data to AnnData object (Scanpy compatible).
    """
    # Load sparse expression matrix
    X = io.mmread(counts_path)
    adata = anndata.AnnData(X=X.transpose().tocsr())

    # Load metadata and assign to adata.obs
    cell_meta = pd.read_csv(meta_path)
    adata.obs = cell_meta
    adata.obs.index = adata.obs['barcode']

    # Load gene names
    with open(gene_name_path, 'r') as f:
        gene_names = f.read().splitlines()
    adata.var.index = gene_names

    # Load PCA coordinates
    pca = pd.read_csv(pca_path, index_col=0)
    pca = pca.loc[adata.obs.index]
    adata.obsm['X_pca'] = pca.to_numpy()

    # Construct 2D layout (UMAP or TSNE)
    adata.obsm[obsm_key] = np.vstack((
        adata.obs[reduction1].to_numpy(),
        adata.obs[reduction2].to_numpy()
    )).T

    return adata


# =======================================
# Example usage
# =======================================

sce_data = seurat_to_adata(
    counts_path='./counts.mtx',
    meta_path='./metadata.csv',
    gene_name_path='./gene_names.csv',
    pca_path='./pca.csv',
    obsm_key='X_umap',
    reduction1='UMAP_1',
    reduction2='UMAP_2'
)

# ---------------------------------------
# Downstream PAGA analysis
# ---------------------------------------

sc.pp.neighbors(sce_data, n_neighbors=20, n_pcs=30)
sc.tl.draw_graph(sce_data, layout='fa')  # optional layout
sc.tl.paga(sce_data, groups='celltype')

# PAGA plot with default and threshold filtering
sc.pl.paga(sce_data, node_size_scale=4.5, edge_width_scale=1.5, node_size_power=1)
sc.pl.paga(sce_data, node_size_scale=4.5, edge_width_scale=1.5, node_size_power=1, threshold=0.05)

# Optional: save PAGA plot to PDF
# sc.pl.paga(sce_data, node_size_scale=4.5, edge_width_scale=1.5, node_size_power=1,
#            threshold=0.05, save='paga_plot.pdf')
