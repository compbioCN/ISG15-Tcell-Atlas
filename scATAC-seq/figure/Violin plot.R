# ------------------------------
# 1. Extract gene expression matrix from ArchR project
# ------------------------------
gene_matrix <- getMatrixFromProject(
  ArchRProj = projT, 
  useMatrix = "GeneIntegrationMatrix", 
  binarize = FALSE
)

# Convert matrix to data frame
expr_data <- as.data.frame(assay(gene_matrix))

# Set rownames to genes (must match gene_matrix rownames)
rownames(expr_data) <- rownames(gene_matrix)

# Extract ISG15 expression vector and transpose to match cell metadata
expr <- as.data.frame(t(expr_data["ISG15", , drop = FALSE]))

# ------------------------------
# 2. Combine with cell metadata
# ------------------------------
cell_metadata <- as.data.frame(projT@cellColData)
expr <- expr[rownames(cell_metadata), , drop = FALSE]
expr_df <- cbind(expr, celltype = cell_metadata$celltype)
colnames(expr_df)[1] <- "ISG15"  # Rename expression column for clarity

# ------------------------------
# 3. Define color palette for cell types
# ------------------------------
library(ggplot2)

celltype_colors <- c(
  "CD4Tstr"     = "#F5BC6E", 
  "CD4T-ISG15"  = "#E05F48", 
  "CD4Tcm"      = "#ACD48A",
  "Th1"         = "#8AB6D6", 
  "CD4Treg"     = "#C7AED5"
)

# ------------------------------
# 4. Create violin + boxplot
# ------------------------------
p <- ggplot(expr_df, aes(x = celltype, y = ISG15, fill = celltype)) +
  geom_violin(trim = FALSE, scale = "width", alpha = 0.9, color = NA) +              # Violin plot
  geom_boxplot(width = 0.15, alpha = 0.85, outlier.shape = NA,                       # Boxplot
               color = "gray30", fill = "white", size = 0.4) +
  scale_fill_manual(values = celltype_colors) +                                      # Cell type colors
  scale_y_continuous(expand = c(0.02, 0)) +                                          # Control Y-axis padding

  # ----------------------------
  # Theme customizations
  # ----------------------------
  theme_minimal(base_size = 12) +
  theme(
    panel.background  = element_rect(fill = "white", color = "black", size = 0.6),
    plot.background   = element_rect(fill = "white", color = "black", size = 0.6),
    panel.grid.major  = element_line(color = "gray90", size = 0.3),
    panel.grid.minor  = element_blank(),
    axis.text.x       = element_text(angle = 45, hjust = 1, size = 11, color = "gray20"),
    axis.text.y       = element_text(size = 10, color = "gray20"),
    axis.title        = element_text(size = 12, face = "bold", color = "gray10"),
    axis.line         = element_line(color = "gray30", size = 0.5),
    legend.position   = "right",  # Change to "none" if you want to hide
    plot.title        = element_text(hjust = 0.5, size = 14, face = "bold", color = "gray10"),
    plot.subtitle     = element_text(hjust = 0.5, size = 11, color = "gray40"),
    plot.margin       = margin(20, 20, 20, 20)
  ) +

  # ----------------------------
  # Axis and label settings
  # ----------------------------
  labs(
    x = "Cell Type",
    y = "ISG15 Gene Integration Score"
  )

# Display plot
print(p)
ggsave("ISG15_violin_plot.pdf", plot = p, width = 6, height = 5, dpi = 600)
