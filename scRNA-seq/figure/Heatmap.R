# =====================================================
# Heatmap of functional gene set scores (CD8+ & CD4+)
# =====================================================

library(Seurat)
library(dplyr)
library(tibble)
library(circlize)
library(ComplexHeatmap)
library(patchwork)

# -------------------------
# Function: Heatmap Plot
# -------------------------

plot_signature_heatmap <- function(seu, score_names, celltype_col = "celltype", row_split_vec, row_annot_vec, title = "Signature Score") {
  
  # Extract data from meta.data
  data <- seu@meta.data[, c(celltype_col, score_names)]
  
  # Average scores per cell type
  avg_data <- data %>%
    group_by(across(all_of(celltype_col))) %>%
    summarise(across(all_of(score_names), mean, na.rm = TRUE)) %>%
    column_to_rownames(celltype_col)
  
  # Transpose and rescale to [-1, 1]
  score_matrix <- t(as.matrix(avg_data))
  score_matrix <- t(apply(score_matrix, 1, rescale, to = c(-1, 1)))
  
  # Color gradient
  col_fun <- colorRamp2(c(-1, 0, 1), c("#ADD8E6", "white", "#C7AED5"))
  
  # Row annotations by functional module
  row_annot <- rowAnnotation(
    Signature = factor(row_annot_vec, levels = c("Differentiation", "Function", "Metabolism", "Apoptosis")),
    col = list(Signature = c(
      Differentiation = "#F7DBF0",
      Function = "#F5BC6E",
      Metabolism = "#E05F48",
      Apoptosis = "#8AB6D6"
    )),
    show_annotation_name = FALSE
  )
  
  # Draw heatmap
  Heatmap(
    score_matrix,
    name = title,
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_split = row_split_vec,
    row_title_gp = gpar(fontsize = 10),
    column_title_gp = gpar(fontsize = 10),
    row_gap = unit(2, "mm"),
    column_gap = unit(2, "mm"),
    top_annotation = NULL,
    left_annotation = row_annot,
    border = TRUE,
    heatmap_legend_param = list(
      title = title,
      at = c(-1, 0, 1),
      labels = c("-1", "0", "1"),
      title_gp = gpar(fontsize = 10),
      labels_gp = gpar(fontsize = 8)
    )
  )
}

# -------------------------
# CD8 module score heatmap
# -------------------------

# gene set names should match AddModuleScore `name` values
cd8_geneset_names <- names(genesetCD8)

cd8_row_split <- c(
  rep("Differentiation", 3),
  rep("Function", 11),
  rep("Metabolism", 3),
  rep("Apoptosis", 2)
)

p1 <- plot_signature_heatmap(
  seu = CD8T,
  score_names = cd8_geneset_names,
  row_split_vec = cd8_row_split,
  row_annot_vec = cd8_row_split,
  title = "CD8⁺ T Cell\nSignature Score"
)

# -------------------------
# CD4 module score heatmap
# -------------------------

cd4_geneset_names <- names(genesetCD4)

cd4_row_split <- c(
  rep("Differentiation", 3),
  rep("Function", 9),
  rep("Metabolism", 3),
  rep("Apoptosis", 2)
)

p2 <- plot_signature_heatmap(
  seu = CD4T,
  score_names = cd4_geneset_names,
  row_split_vec = cd4_row_split,
  row_annot_vec = cd4_row_split,
  title = "CD4⁺ T Cell\nSignature Score"
)

# -------------------------
# Combine & plot
# -------------------------

draw(p1 / p2)
