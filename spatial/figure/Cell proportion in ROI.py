# -----------------------------------------
# Extract and Visualize ROI in Spatial Data
# -----------------------------------------

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import scanpy as sc

# ---------- ROI Configuration ----------
roi_center_x = 800
roi_center_y = 2596
roi_radius = 200
sample_id = list(adata.uns['spatial'].keys())[0]
scalefactor = adata.uns['spatial'][sample_id]['scalefactors']['tissue_hires_scalef']

# ---------- Compute Spatial ROI ----------
roi_center_x_spatial = roi_center_x / scalefactor
roi_center_y_spatial = roi_center_y / scalefactor
roi_radius_spatial = roi_radius / scalefactor

coordinates = adata.obsm['spatial']
distances = np.sqrt((coordinates[:, 0] - roi_center_x_spatial) ** 2 +
                    (coordinates[:, 1] - roi_center_y_spatial) ** 2)
in_roi = distances <= roi_radius_spatial

# ---------- Subset and Save ----------
cells_in_roi = np.sum(in_roi)
if cells_in_roi == 0:
    print("⚠️ No cells found in ROI.")
    exit()

roi_subset = adata[in_roi].copy()
subset_filename = f"ROI_subset_center_{roi_center_x}_{roi_center_y}_radius_{roi_radius}.h5ad"
roi_subset.write(subset_filename)
print(f"✅ ROI subset saved: {subset_filename} ({roi_subset.n_obs} cells)")

# ---------- Visualize ROI on Tissue Image ----------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 8))

# Plot full image
sc.pl.spatial(adata, img_key="hires", alpha_img=0.7, size=2, show=False, ax=ax1)
circle1 = Circle((roi_center_x, roi_center_y), roi_radius, edgecolor='orange',
                 facecolor='none', linewidth=3, linestyle='--', alpha=0.9)
ax1.add_patch(circle1)
ax1.set_title('Full Image with ROI', fontsize=14, fontweight='bold')

# Plot cropped ROI
sc.pl.spatial(roi_subset, img_key="hires", alpha_img=0.7, size=5, show=False, ax=ax2)
ax2.set_title(f'ROI Area ({cells_in_roi} cells)', fontsize=14, fontweight='bold')
margin = roi_radius * 0.5
ax2.set_xlim(roi_center_x - roi_radius - margin, roi_center_x + roi_radius + margin)
ax2.set_ylim(roi_center_y + roi_radius + margin, roi_center_y - roi_radius - margin)

plt.tight_layout()
plt.savefig(f"ROI_visualization_center_{roi_center_x}_{roi_center_y}_radius_{roi_radius}.pdf", format='pdf', dpi=300)
plt.close()
