# ==================================
# Violin plots of functional scores by cell type
# - Plot Exhaustion and Cytotoxicity scores
# - Box-style violin, per score panel
# ==================================

import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns

# Set plotting style
sns.set(style="whitegrid")

# Create violin plots per score
sc.pl.violin(
    adata,
    keys=["Exhaustion_score", "Cytotoxicity_score"],  # Scores to visualize
    groupby="cell_type",        # Group by cell type
    stripplot=False,            # No dots
    jitter=False,               # No jitter
    inner="box",                # Add boxplot inside violins
    scale="width",              # Uniform width
    rotation=45,                # Rotate x-axis labels
    linewidth=1,                # Border width
    multi_panel=True,           # Separate panels per score
    show=False                  # Delay rendering for possible post-editing
)

# Optional: adjust spacing
plt.tight_layout()
plt.show()
