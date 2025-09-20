import os
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import bin2cell as b2c

# --------------------------
# Paths & Parameters
# --------------------------
data_root = "YOUR_PROJECT_FOLDER"
bin_path = f"{data_root}/binned_outputs/square_002um"
he_image_path = f"{data_root}/HE_image.tiff"
gex_image_path = "stardist/gex.tiff"
spaceranger_image_path = f"{data_root}/spatial"
os.makedirs("stardist", exist_ok=True)
mpp = 0.5  # microns per pixel

# --------------------------
# Load & Filter Raw Data
# --------------------------
adata = b2c.read_visium(
    bin_path,
    source_image_path=he_image_path,
    spaceranger_image_path=spaceranger_image_path
)
adata.var_names_make_unique()
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.filter_cells(adata, min_counts=1)

# --------------------------
# HE segmentation (nuclei)
# --------------------------
b2c.destripe(adata)
b2c.scaled_he_image(adata, mpp=mpp, save_path="stardist/he.tiff")

b2c.stardist(
    image_path="stardist/he.tiff",
    labels_npz_path="stardist/he.npz",
    stardist_model="2D_versatile_he",
    prob_thresh=0.01
)

b2c.insert_labels(
    adata,
    labels_npz_path="stardist/he.npz",
    basis="spatial",
    spatial_key="spatial",
    mpp=mpp,
    labels_key="labels_he"
)

# Filter and visualize
bdata = adata[adata.obs["labels_he"] > 0].copy()
bdata.obs["labels_he"] = bdata.obs["labels_he"].astype(str)
crop = b2c.get_crop(bdata, basis="spatial", spatial_key="spatial", mpp=mpp)

rendered = b2c.view_stardist_labels(
    image_path="stardist/he.tiff",
    labels_npz_path="stardist/he.npz",
    crop=crop
)
plt.imshow(rendered)
plt.savefig("he_nuclei_mask2.png", dpi=1200, bbox_inches="tight")
plt.close()

# --------------------------
# Label Expansion
# --------------------------
b2c.expand_labels(
    adata,
    labels_key="labels_he",
    expanded_labels_key="labels_he_expanded"
)

# --------------------------
# GEX segmentation (adjusted expression image)
# --------------------------
b2c.grid_image(adata, "n_counts_adjusted", mpp=mpp, sigma=5, save_path="stardist/gex.tiff")

b2c.stardist(
    image_path="stardist/gex.tiff",
    labels_npz_path="stardist/gex.npz",
    stardist_model="2D_versatile_fluo",
    prob_thresh=0.05,
    nms_thresh=0.5
)

b2c.insert_labels(
    adata,
    labels_npz_path="stardist/gex.npz",
    basis="array",
    mpp=mpp,
    labels_key="labels_gex"
)

# Filter and visualize
bdata = adata[adata.obs["labels_gex"] > 0].copy()
bdata.obs["labels_gex"] = bdata.obs["labels_gex"].astype(str)
crop = b2c.get_crop(bdata, basis="array", mpp=mpp)

rendered = b2c.view_stardist_labels(
    image_path="stardist/gex.tiff",
    labels_npz_path="stardist/gex.npz",
    crop=crop
)
plt.imshow(rendered)
plt.savefig("gex_cell_mask2.png", dpi=1200, bbox_inches="tight")
plt.close()

# --------------------------
# Combine HE + GEX labels
# --------------------------
b2c.salvage_secondary_labels(
    adata,
    primary_label="labels_he_expanded",
    secondary_label="labels_gex",
    labels_key="labels_joint"
)

# --------------------------
# Generate Cell-Level AnnData
# --------------------------
bdata = adata[adata.obs["labels_joint"] > 0].copy()
bdata.obs["labels_joint"] = bdata.obs["labels_joint"].astype(str)

cdata = b2c.bin_to_cell(
    adata,
    labels_key="labels_joint",
    spatial_keys=["spatial", "spatial_cropped_150_buffer"]
)

# Save output
cdata.write_h5ad("b2c.h5ad")
