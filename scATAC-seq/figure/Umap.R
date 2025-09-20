# Load required libraries
library(ggplot2)
library(grid)

# Define custom color palette for cell types
colors <- c(
  "B cells" = "#EA945A", 
  "Endothelial cells" = "#86A667", 
  "Epithelial cells" = "#E05F48",
  "Fibroblasts" =  "#60A897",
  "Mast cells" =  "#6488B9",
  "Mural cells" = "#ACD48A",
  "Myeloid cells" = "#8AB6D6", 
  "NK" =  "#726BAE",
  "Plasma cells" =  "#F5BC6E",
  "Proliferative T cells" = "#E5948E",
  "T" = "#C7AED5"
)

# Extract metadata and UMAP coordinates
data <- as.data.frame(sc@cellColData)
A <- sc@embeddings$UMAP$df
stopifnot(identical(rownames(data), rownames(A)))  # Ensure consistency

# Construct UMAP data frame
umap_data <- data.frame(
  UMAP1 = A[, 1],
  UMAP2 = A[, 2],
  CellType = data[["predictedGroup_ArchR"]]
)

# Base UMAP scatter plot
p <- ggplot(umap_data, aes(x = UMAP1, y = UMAP2, color = CellType)) +
  geom_point(size = 0.01, alpha = 0.15) +
  scale_color_manual(values = colors) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_blank()
  )

# Add cell type labels at cluster centers
label_positions <- aggregate(cbind(UMAP1, UMAP2) ~ CellType, data = umap_data, FUN = mean)
p <- p + geom_text(
  data = label_positions, 
  aes(x = UMAP1, y = UMAP2, label = CellType), 
  color = "black", size = 4, fontface = "bold"
)

# Add total cell count to top-left corner
cell_count <- nrow(umap_data)
p <- p + annotate("text", 
                  x = min(umap_data$UMAP1), 
                  y = max(umap_data$UMAP2), 
                  label = paste0("k=", format(cell_count, big.mark = ","), " cells"),
                  hjust = -0.1, vjust = 1.5, size = 4, fontface = "bold")

# Add UMAP axis arrows in bottom-left corner
arrow_base_x <- min(umap_data$UMAP1) + 1
arrow_base_y <- min(umap_data$UMAP2) + 1
p <- p + 
  geom_segment(aes(x = arrow_base_x, y = arrow_base_y, 
                   xend = arrow_base_x + 2, yend = arrow_base_y), 
               arrow = arrow(length = unit(0.3, "cm")), 
               size = 0.6, color = "black") +
  geom_segment(aes(x = arrow_base_x, y = arrow_base_y, 
                   xend = arrow_base_x, yend = arrow_base_y + 2), 
               arrow = arrow(length = unit(0.3, "cm")), 
               size = 0.6, color = "black") +
  annotate("text", x = arrow_base_x + 2.5, y = arrow_base_y, 
           label = "UMAP1", hjust = 0, size = 4, fontface = "bold") +
  annotate("text", x = arrow_base_x, y = arrow_base_y + 2.5, 
           label = "UMAP2", vjust = 0, size = 4, fontface = "bold")

# Display the plot
print(p)

# Save the plot (optional)
ggsave("UMAP_labeled_arrow.pdf", plot = p, width = 6, height = 6, dpi = 600)
