# ===============================
# Slingshot Trajectory Inference
# ===============================

library(Seurat)
library(SingleCellExperiment)
library(slingshot)
library(ggplot2)
library(dplyr)
library(patchwork)

# Convert Seurat to SingleCellExperiment
sce <- as.SingleCellExperiment(CD4T, assay = "RNA")

# Run slingshot (adjust start.clus if needed)
sce <- slingshot(
  sce,
  reducedDim = "UMAP",
  clusterLabels = sce$celltype,
  start.clus = "CD4Tcm",  # Optional
  approx_points = 150
)

# Extract UMAP coordinates
umap_coords <- as.data.frame(reducedDims(sce)$UMAP)
colnames(umap_coords) <- c("UMAP1", "UMAP2")
umap_coords$celltype <- sce$celltype

# Extract lineage curves
sds <- SlingshotDataSet(sce)
curves <- slingCurves(sds)

# Combine all curve points into one dataframe
trajectory_data <- do.call(rbind, lapply(seq_along(curves), function(i) {
  curve_points <- curves[[i]]$s
  data.frame(
    UMAP1 = curve_points[, 1],
    UMAP2 = curve_points[, 2],
    Lineage = paste0("Lineage", i),
    Order = seq_len(nrow(curve_points))
  )
}))

# =====================================
# Visualize pseudotime over UMAP
# =====================================

plot_lineage <- function(sce, umap_df, lineage_id, curve_df) {
  ggplot(umap_df, aes(x = UMAP1, y = UMAP2)) +
    geom_point(aes(color = sce[[paste0("slingPseudotime_", lineage_id)]]), size = 0.8) +
    geom_path(data = subset(curve_df, Lineage == paste0("Lineage", lineage_id)),
              aes(x = UMAP1, y = UMAP2), color = "black", size = 1.5) +
    scale_color_viridis_c(name = "Pseudotime", na.value = "grey90") +
    theme_bw() +
    ggtitle(paste("Lineage", lineage_id)) +
    coord_fixed()
}

# Plot each lineage
p1 <- plot_lineage(sce, umap_coords, 1, trajectory_data)
p2 <- plot_lineage(sce, umap_coords, 2, trajectory_data)
p3 <- plot_lineage(sce, umap_coords, 3, trajectory_data)
p4 <- plot_lineage(sce, umap_coords, 4, trajectory_data)

# Combine (optional)
(p1 | p2) / (p3 | p4)

# =====================================
# Add pseudotime to Seurat metadata
# =====================================

pseudotime_df <- slingPseudotime(sce) %>% as.data.frame()
for (lineage in colnames(pseudotime_df)) {
  CD4T <- AddMetaData(CD4T, metadata = pseudotime_df[[lineage]], col.name = lineage)
}

# =====================================
# Pseudotime Density Plot per Lineage
# =====================================

# Define custom color palette
custom_colors <- c("#ACD48A", "#C7AED5", "#EA945A", "#8AB6D6", "#E05F48")

# Function to plot density of pseudotime
pseudotime_density <- function(seurat_obj, lineage_col, cluster_col, fill_colors) {
  df <- data.frame(
    celltype = seurat_obj[[cluster_col]][, 1],
    pseudotime = seurat_obj[[lineage_col]][, 1]
  ) %>% na.omit()

  ggplot(df, aes(x = pseudotime, fill = celltype)) +
    geom_density(alpha = 0.5) +
    scale_fill_manual(values = fill_colors) +
    theme_bw() +
    labs(title = paste("Pseudotime:", lineage_col), x = "Pseudotime", y = "Density")
}

# Example usage: Lineage 1
p_density1 <- pseudotime_density(CD4T, "Lineage1", "celltype", custom_colors)
print(p_density1)
