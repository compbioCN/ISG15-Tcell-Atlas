# ======================================================
# Step 1: Export count matrix for scFEA input (CD8/CD4)
# ======================================================

library(Seurat)

# Export count matrix
count_matrix <- CD8T@assays$RNA@counts
write.csv(count_matrix, file = "scFEA/input/count_matrix.csv", row.names = TRUE)

# Then run scFEA (Python command):
# ---------------------------------
# python scFEA.py \
#   --data_dir data \
#   --input_dir input \
#   --moduleGene_file module_gene_m168.csv \
#   --test_file count_matrix.csv \
#   --cName_file cName_c70_m168.csv \
#   --sc_imputation True \
#   --stoichiometry_matrix cmMat_c70_m168.csv \
#   --output_flux_file output/adj_flux.csv \
#   --output_balance_file output/adj_balance.csv

# ======================================================
# Step 2: Read flux output and calculate group averages
# ======================================================

library(dplyr)
library(ComplexHeatmap)
library(pheatmap)

# Load flux data from scFEA
adj_flux <- read.csv("scFEA/CD4/output/adj_flux.csv", row.names = 1)
rownames(adj_flux) <- gsub("\\.", "-", rownames(adj_flux))  # standardize cell IDs

# Load Seurat object
adj_scRNA <- CD4T  # or CD8T if needed

# Create combined group label
adj_scRNA$group_cells <- paste0(adj_scRNA$source, "_", adj_scRNA$celltype)

# Build cell annotation mapping
cell_anno <- data.frame(
  cellid = rownames(adj_scRNA@meta.data),
  group_cells = adj_scRNA$group_cells
)
cell_anno$cellid <- gsub("@", "-", cell_anno$cellid)  # ensure ID match
cell_anno <- cell_anno[order(cell_anno$group_cells), ]

# Align flux data to annotation
adj_flux <- adj_flux[cell_anno$cellid, ]

# Compute mean flux for each group
df_averages <- adj_flux %>%
  group_by(group = cell_anno$group_cells) %>%
  summarise_all(mean, na.rm = TRUE) %>%
  select(-group)

rownames(df_averages) <- unique(cell_anno$group_cells)
df_averages <- t(df_averages) %>% as.data.frame()

# ======================================================
# Step 3: Subset CD4 T subtypes of interest for heatmap
# ======================================================

selected_groups <- c(
  "Precancerous_CD4T-ISG15", "Tumor_CD4T-ISG15",
  "Precancerous_CD4Tcm",     "Tumor_CD4Tcm",
  "Precancerous_CD4Treg",    "Tumor_CD4Treg",
  "Precancerous_CD4Tstr",    "Tumor_CD4Tstr",
  "Precancerous_Th1",        "Tumor_Th1"
)

df_flux <- df_averages[, selected_groups]

# Remove rows with zero standard deviation
df_flux <- df_flux[apply(df_flux, 1, function(x) sd(x) != 0), ]

# Replace NA with 0
df_flux[is.na(df_flux)] <- 0

# ======================================================
# Step 4: Draw heatmap
# ======================================================

pheatmap(
  df_flux,
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  show_rownames = FALSE,
  scale = "row",
  color = colorRampPalette(c(
    "#2166AC", "#478ABF", "#90C0DC",
    "white",
    "#EF8C65", "#CF4F45", "#B2182B"
  ))(100),
  border = FALSE,
  heatmap_legend_param = list(title = "Flux"),
  fontsize_col = 8,
  treeheight_row = 20
)
