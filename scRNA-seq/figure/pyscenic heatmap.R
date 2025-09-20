# =========================================
# SCENIC RSS Heatmap: Fibroblast subtypes
# =========================================

library(SCENIC)
library(SCopeLoomR)
library(AUCell)
library(dplyr)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(grid)

# -----------------------------
# Step 1: Load AUC and metadata
# -----------------------------

sce_SCENIC <- open_loom("sample_SCENIC.loom")  # path to pySCENIC output
regulonAUC <- get_regulons_AUC(sce_SCENIC, column.attr.name = "RegulonsAUC")

# Seurat metadata
cellinfo <- CD4T@meta.data[, c("celltype", "Tissue", "nFeature_RNA", "nCount_RNA")]
colnames(cellinfo) <- c("celltype", "group", "nGene", "nUMI")
cellTypes <- cellinfo["celltype"]
selectedResolution <- "celltype"

# -----------------------------
# Step 2: Calculate RSS and AUC avg
# -----------------------------

rss <- calcRSS(
  AUC = getAUC(regulonAUC),
  cellAnnotation = cellTypes[colnames(regulonAUC), selectedResolution]
)
rss <- na.omit(rss)

cellsPerGroup <- split(rownames(cellTypes), cellTypes[[selectedResolution]])
regulonActivity_byGroup <- sapply(cellsPerGroup, function(cells) {
  rowMeans(getAUC(regulonAUC)[, cells])
})

regulonActivity_byGroup_Scaled <- t(scale(t(regulonActivity_byGroup)))
regulonActivity_byGroup_Scaled <- na.omit(regulonActivity_byGroup_Scaled)

# -----------------------------
# Step 3: Extract top N regulons per group
# -----------------------------

top_n <- 5
group_names <- colnames(regulonActivity_byGroup_Scaled)

top_regulons_list <- lapply(group_names, function(group) {
  top_genes <- head(order(regulonActivity_byGroup_Scaled[, group], decreasing = TRUE), top_n)
  regulonActivity_byGroup_Scaled[top_genes, , drop = FALSE]
})
top_rows_unique <- unique(do.call(rbind, top_regulons_list))

# -----------------------------
# Step 4: Annotations
# -----------------------------

row_group_assignment <- apply(top_rows_unique, 1, function(x) {
  group_names[which.max(x)]
})
gene_anno_df <- data.frame(RegulonGroup = row_group_assignment)
rownames(gene_anno_df) <- rownames(top_rows_unique)

group_palette <- setNames(
  RColorBrewer::brewer.pal(n = max(3, length(group_names)), name = "Set2")[1:length(group_names)],
  group_names
)

row_anno <- rowAnnotation(
  Group = gene_anno_df$RegulonGroup,
  col = list(Group = group_palette),
  show_annotation_name = FALSE
)

cell_anno <- data.frame(cell = colnames(top_rows_unique))
col_anno <- HeatmapAnnotation(
  df = cell_anno,
  col = list(cell = group_palette),
  show_annotation_name = FALSE,
  gp = gpar(col = "white", lwd = 2)
)

# -----------------------------
# Step 5: Plot heatmap
# -----------------------------

col_fun <- colorRamp2(c(-1.5, 0, 1.5), c("#0da9ce", "white", "#e74a32"))

Heatmap(
  as.matrix(top_rows_unique),
  name = "RSS-Z Score",
  col = col_fun,
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  column_names_side = "top",
  column_names_rot = 60,
  row_names_gp = gpar(fontsize = 12, fontface = "italic"),
  rect_gp = gpar(col = "white", lwd = 1.5),
  top_annotation = col_anno
) + row_anno
