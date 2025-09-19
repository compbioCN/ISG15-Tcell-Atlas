# ===============================
# CD4+ T cell subtype flow plot across tissue sources (Alluvial plot)
# ===============================

# Load required libraries
library(Seurat)       # For Idents() if CD4T is a Seurat object
library(ggplot2)      # For plotting
library(ggalluvial)   # For geom_flow and stratum-based plotting

# Define color mapping for cell types
colors <- c(
  "CD4Tstr"              = "#EA945A", 
  "Endothelial cells"    = "#86A667", 
  "CD4T-ISG15"           = "#E05F48",
  "Fibroblasts"          = "#60A897",
  "CD8Trm"               = "#6488B9",
  "CD4Tcm"               = "#ACD48A",
  "Th1"                  = "#8AB6D6", 
  "CD8Tcm"               = "#726BAE",
  "Plasma cells"         = "#F5BC6E",
  "Proliferative T cells"= "#E5948E",
  "CD4Treg"              = "#C7AED5"
)

# Calculate proportions of each CD4T identity per source group
cell_ratio <- prop.table(table(Idents(CD4T), CD4T$source), margin = 2)
cell_ratio <- as.data.frame(cell_ratio)
colnames(cell_ratio) <- c("celltype", "group", "ratio")

# Create alluvial (flow) plot
p <- ggplot(cell_ratio, aes(
  x = group,
  y = ratio,
  fill = celltype,
  stratum = celltype,
  alluvium = celltype
)) +
  geom_col(width = 0.4, color = NA) +                                # stacked bar
  geom_flow(width = 0.4, alpha = 0.4, knot.pos = 0) +               # smooth flow
  scale_fill_manual(values = colors) +                              # cell type colors
  theme_classic() +
  theme(
    legend.position = "none",                                       # hide legend
    panel.border = element_rect(
      fill = NA, color = "black", size = 0.5, linetype = "solid"
    )
  )

# Show plot
print(p)

# Optional: Save to file
ggsave("CD4T_source_alluvial.pdf", p, width = 6, height = 4, dpi = 300)
