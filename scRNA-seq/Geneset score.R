# ========================================================
# Module scoring of gene sets in CD8+ and CD4+ T cells
# ========================================================

library(Seurat)
library(dplyr)
library(tibble)

# -------- CD8 T cells module scoring --------

# Remove NA genes from each CD8 gene set
geneset_list_cd8 <- lapply(genesetCD8, na.omit)

# Add module score to Seurat object for each gene set
for (set_name in names(geneset_list_cd8)) {
  gene_set <- geneset_list_cd8[[set_name]]
  CD8T <- AddModuleScore(
    object = CD8T,
    features = list(gene_set),
    name = set_name,
    assay = "RNA"
  )
}

# -------- CD4 T cells module scoring --------

geneset_list_cd4 <- lapply(genesetCD4, na.omit)

for (set_name in names(geneset_list_cd4)) {
  gene_set <- geneset_list_cd4[[set_name]]
  CD4T <- AddModuleScore(
    object = CD4T,
    features = list(gene_set),
    name = set_name,
    assay = "RNA"
  )
}

# -------- Summary statistics for CD8 T cells --------

# Use known gene set score column names directly
score_cols <- names(geneset_list_cd8)
data <- CD8T@meta.data[, c("celltype", score_cols)]

# Compute average module score for each cell type
avg_data <- data %>%
  group_by(celltype) %>%
  summarise(across(everything(), mean, na.rm = TRUE))

# Convert to matrix-style dataframe (celltype as rownames)
avg_data_long <- avg_data %>%
  column_to_rownames("celltype")
