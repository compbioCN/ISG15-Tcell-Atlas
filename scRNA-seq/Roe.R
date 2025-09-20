# === Load Required Packages ===
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(Startrac)
library(ComplexHeatmap)
library(circlize)

# === Step 1: Prepare Metadata ===
# Extract required metadata from Seurat object
BigPDT <- data  # Assume your Seurat object is called `data`
meta_df <- BigPDT@meta.data[, c("Sample", "TimePoint", "celltype")]
colnames(meta_df) <- c("sample", "tissue", "celltype")

# Ensure tissue order is consistent for Ro/e calculation
meta_df$tissue <- factor(meta_df$tissue, levels = c("Non-aPD1", "aPD1"))

# === Step 2: Calculate Ro/e Matrix ===
Roe <- calTissueDist(
  meta_df,
  byPatient = FALSE,
  colname.cluster = "celltype",
  colname.patient = "sample",
  colname.tissue = "tissue",
  method = "chisq"
)

# Replace NA with 0 and cap values at 2
Roe[is.na(Roe)] <- 0
Roe_capped <- pmin(Roe, 2)

# === Step 3: Create Ro/e Heatmap ===
# Define custom color scale
col_fun <- colorRamp2(c(0, 1, 2), c("#8AB6D6", "white", "#E05F48"))

# Generate heatmap object
p_roe <- Heatmap(
  as.matrix(Roe_capped),
  name = "Ro/e",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_title_gp = gpar(fontsize = 10),
  column_title_gp = gpar(fontsize = 10),
  row_names_side = "left",
  column_names_rot = 45,
  gap = unit(3, "mm"),
  rect_gp = gpar(col = "white", lwd = 0.8),
  border = TRUE,
  border_gp = gpar(col = "black", lwd = 1),
  heatmap_legend_param = list(
    at = c(0, 1, 2),
    labels = c("0", "1", "2+"),
    title_gp = gpar(fontsize = 10),
    labels_gp = gpar(fontsize = 8)
  )
)

# === Step 4: Draw Heatmap ===
draw(p_roe)
