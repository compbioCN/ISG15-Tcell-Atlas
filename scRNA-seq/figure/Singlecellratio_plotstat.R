library(Seurat)
library(ggplot2)
library(reshape2)
library(ggpubr)
Singlecellratio_plotstat <- function(seu,
                                     by = "cell.type",
                                     meta.include = NULL,
                                     group_by = NULL,
                                     shape_by = NULL,
                                     custom_fill_colors = c("#8AB6D6", "#F5BC6E", "#E05F48"),
                                     group_by.point = NULL,
                                     color_by = NULL,
                                     pb = FALSE,
                                     comparisons = my_comparisons,
                                     ncol = NULL,
                                     label = c("p.format", "p.signif"),
                                     label.x = NA,
                                     pt.size = 4) {
  
  by <- match.arg(by)
  if (is.null(group_by)) group_by <- "null.group"
  
  shapes <- if (!is.null(shape_by)) c(16, 15, 3, 7, 8, 18, 5, 6, 2, 4, 1, 17) else NULL
  
  # Compute % composition of each cell type per sample
  fq <- prop.table(table(seu@meta.data$celltype, seu@meta.data$orig.ident), margin = 2) * 100
  df <- melt(fq, value.name = "freq", varnames = c("cell.type", "orig.ident"))
  
  # Get unique group/sample metadata
  uniques <- apply(seu@meta.data, 2, function(x) length(unique(x)))
  ei <- unique(seu@meta.data[, names(uniques[uniques <= 100])])
  ei <- unique(ei[, colnames(ei) %in% meta.include])
  df <- merge(df, ei, by = "orig.ident")
  df$null.group <- "1"
  df$orig.ident <- as.factor(df$orig.ident)
  
  if (is.null(ncol)) {
    n_types <- length(unique(df$cell.type))
    ncol <- ifelse(n_types > 20, 5, ifelse(n_types > 9, 4, 3))
  }

  # Basic ggplot setup
  p <- ggplot(df, aes_string(y = "freq", x = group_by)) +
    labs(x = NULL, y = "Proportion (%)") +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 12),
      axis.ticks.x = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black", size = 12),
      axis.text.y = element_text(color = "black", size = 12),
      axis.title.y = element_text(size = 14, face = "bold")
    )
  
  # Plot options
  if (by == "cell.type" && color_by == "cell.type") {
    p + 
      facet_wrap(group_by, scales = "free_x") +
      geom_bar(aes_string(x = "orig.ident", fill = "factor(cell.type)"),
               position = "fill", stat = "identity") +
      scale_fill_manual("Cell Type", values = custom_fill_colors) +
      scale_y_continuous(expand = c(0, 0), labels = seq(0, 100, 25)) +
      theme(panel.border = element_blank())
    
  } else {
    p + 
      facet_wrap("cell.type", scales = "free_y", ncol = ncol) +
      guides(fill = FALSE) +
      geom_boxplot(aes_string(x = group_by), alpha = 0.25, outlier.color = NA) +
      geom_point(
        size = pt.size,
        position = position_jitter(width = 0.25),
        aes_string(x = group_by, y = "freq", color = color_by, shape = shape_by)
      ) +
      scale_shape_manual(values = shapes) +
      theme(panel.grid.major = element_line(color = "grey", size = 0.25)) +
      scale_color_manual(values = custom_fill_colors) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
      ggpubr::stat_compare_means(
        mapping = aes_string(group_by),
        comparisons = comparisons,
        label = label,
        method = "t.test"
      )
  }
}
