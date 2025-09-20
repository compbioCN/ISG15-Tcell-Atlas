# ROI Extraction and Visualization for Spatial Transcriptomics (Scanpy)

import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Polygon, Point
from matplotlib.patches import Polygon as MplPolygon
import scanpy as sc

# === Load spatial coordinates and scale factor ===
coordinates = CD4Tcell.obsm['spatial']
scalefactor = CD4Tcell.uns['spatial']['137486A_outs']['scalefactors']['tissue_hires_scalef']

# === Define polygon ROI (in image coordinates, then scaled) ===
polygon_points = np.array([
    [1900, 4000],
    [1900, 8700],
    [6900, 8700],
    [6900, 4000],
    [1900, 4000]
])
scaled_polygon_points = polygon_points / scalefactor
polygon = Polygon(scaled_polygon_points)

# === Select cells inside the polygon ===
in_region = [polygon.contains(Point(coord)) for coord in coordinates]
CD4Tcell_subset = CD4Tcell[in_region].copy()

# === Spatial plot of CD4T cells with highlighted ROI ===
sc.pl.spatial(CD4Tcell, color='Tcelltype', size=3, show=False, legend_loc=None, alpha=0.8, alpha_img=0.5)
ax = plt.gca()
polygon_points_image = polygon_points  # original image coordinates
mpl_polygon = MplPolygon(polygon_points_image, closed=True, edgecolor='red', fill=False, linewidth=2, label='Region of Interest')
ax.add_patch(mpl_polygon)
ax.legend(handles=[mpl_polygon], loc='upper right')
plt.title('CD4Tcell Spatial Plot with Highlighted ROI')
plt.xlabel('X Coordinate')
plt.ylabel('Y Coordinate')
plt.grid(False)
plt.savefig("CD4Tcell_ROI_highlighted_no_annotation.pdf", format='pdf', bbox_inches='tight')
plt.show()

# === Plot ISG15 expression in CD8T ROI subset ===
plt.figure(figsize=(12, 12))
sc.pl.spatial(
    CD8Tcell_subset,
    color="ISG15",
    size=3,
    alpha=0.8,
    alpha_img=0.5,
    frameon=False,
    cmap="Reds",
    show=False,
    vmin=0,
    vmax=5
)
plt.savefig("ISG15_CD8T_cell_type.pdf", bbox_inches='tight', dpi=900)
plt.close()

# === Plot T cell types in CD8T ROI subset ===
plt.figure(figsize=(12, 12))
sc.pl.spatial(
    CD8Tcell_subset,
    color="Tcelltype",
    size=5,
    alpha=0.8,
    alpha_img=0.5,
    frameon=False,
    cmap="Spectral",
    show=False
)
plt.savefig("CD8T_cell_type.pdf", bbox_inches='tight', dpi=900)
plt.close()
