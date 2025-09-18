# UMAP visualization and export (high-resolution PDF)
import os
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt

# Create UMAP plot
fig, ax = plt.subplots(figsize=(6, 6), dpi=300)  # High resolution

sc.pl.umap(
    adata,
    color='cell_type',     # Column in adata.obs to color by
    size=0.1,              # Point size
    alpha=0.8,             # Transparency
    frameon=False,         # Remove outer frame
    ax=ax,
    show=False             # Prevent immediate display
)

# Remove axis ticks and labels
ax.set_xticks([])
ax.set_yticks([])
ax.set_xticklabels([])
ax.set_yticklabels([])

# Save as PDF
plt.savefig("umap_plot.pdf", format="pdf")

# Display plot in notebook or interactive session
plt.show()

