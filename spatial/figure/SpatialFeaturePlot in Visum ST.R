# === Load Required Libraries ===
library(Seurat)
library(ggplot2)
library(patchwork)

# === Define Unified Color Gradient ===
score_colors <- c(
  "#333333", "#1B0C42FF", "#4B0C6BFF", "#781C6DFF", "#A52C60FF",
  "#CF4446FF", "#ED6925FF", "#FB9A06FF", "#F7D03CFF", "#FCFFA4FF"
)

# === Loop Through Each Spatial Seurat Object ===
for (i in seq_along(spatial_list)) {
  sample_name <- names(spatial_list)[i]
  sample_seurat <- spatial_list[[i]]

  # --- Plot 1: Spatial Cell Type Distribution ---
  p1 <- SpatialDimPlot(
    object = sample_seurat,
    label = FALSE,
    pt.size.factor = 8000,
    image.alpha = 0.6,
    crop = FALSE,
    stroke = 0.1,
    shape = 21,
    label.box = TRUE,
    label.size = 6,
    label.color = "white"
  ) &
    scale_fill_manual(values = colors) &
    theme_bw() &
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
    )

  # --- Plot 2: CD4 ISG15 Score Heatmap ---
  p2 <- SpatialFeaturePlot(
    sample_seurat,
    features = "CD4ISG15_Score",
    pt.size.factor = 8000,
    image.alpha = 0.0
  ) +
    scale_fill_gradientn(
      colours = score_colors,
      limits = c(0.0, 1.2),
      breaks = c(0.0, 0.6, 1.2),
      na.value = "black",
      name = "CD4ISG15_Score"
    ) +
    ggtitle("CD4ISG15 Score") +
    theme_void(base_size = 14) +
    Seurat::DarkTheme() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold", color = "white"),
      panel.border = element_rect(color = "white", fill = NA, size = 2),
      legend.title = element_text(color = "white"),
      legend.text = element_text(color = "white")
    )

  # --- Plot 3: CD8 ISG15 Score Heatmap ---
  p3 <- SpatialFeaturePlot(
    sample_seurat,
    features = "CD8ISG15_Score",
    pt.size.factor = 8000,
    image.alpha = 0.0
  ) +
    scale_fill_gradientn(
      colours = score_colors,
      limits = c(0.0, 1.2),
      breaks = c(0.0, 0.6, 1.2),
      na.value = "black",
      name = "CD8ISG15_Score"
    ) +
    ggtitle("CD8ISG15 Score") +
    theme_void(base_size = 14) +
    Seurat::DarkTheme() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold", color = "white"),
      panel.border = element_rect(color = "white", fill = NA, size = 2),
      legend.title = element_text(color = "white"),
      legend.text = element_text(color = "white")
    )

  # --- Combine and Save ---
  combined_plot <- p1 | p2 | p3

  ggsave(
    filename = paste0("Spatial_ISG15_", sample_name, ".png"),
    plot = combined_plot,
    dpi = 1200,
    width = 16,
    height = 6
  )
}
