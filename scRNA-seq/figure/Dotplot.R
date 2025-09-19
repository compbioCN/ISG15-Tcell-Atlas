library(Seurat)
library(ggplot2)

# Define marker gene sets per cell type
features <- list(
  "B cells"               = c("MS4A1", "BANK1", "TNFRSF13C"),
  "Endothelial cells"     = c("FLT1", "SELE", "TM4SF1"),
  "Epithelial cells"      = c("KRT14", "KRT13", "KRT5"),
  "Fibroblasts"           = c("DCN", "COL1A1", "LUM"),
  "Mast cells"            = c("TPSB2", "CPA3", "TPSAB1"),
  "Mural cells"           = c("TAGLN", "ACTA2", "RGS5"),
  "Myeloid cells"         = c("LYZ", "IL1B", "C1QA"),
  "NK"                    = c("GNLY", "NKG7", "GZMB"),
  "Plasma cells"          = c("IGLC2", "IGKC", "IGHG1"),
  "Proliferative T cells" = c("MKI67", "STMN1", "TOP2A"),
  "T"                     = c("CD8A", "CD8B", "GZMK")
)

# Generate DotPlot with full customization
DotPlot(object = OSCC, features = features) &
  theme_bw(base_size = 12) &
  geom_point(shape = 21, aes(size = pct.exp), stroke = 1.2) &
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(
      color = 'black', angle = 45, hjust = 1, vjust = 1,
      face = "bold", size = 10
    ),
    axis.text.y = element_text(color = 'black', face = "bold", size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 10),
    plot.margin = unit(c(1, 1, 1, 1), 'cm'),
    panel.border = element_rect(color = "black", size = 1.2),
    panel.spacing = unit(0.2, "cm"),
    legend.box.background = element_rect(colour = "black", size = 0.5),
    legend.key.width = unit(0.4, "cm"),
    legend.key.height = unit(0.6, "cm"),
    legend.title = element_text(color = 'black', face = "bold", size = 11),
    legend.text = element_text(size = 10),
    plot.tag.position = c(0.3, 1.1),
    plot.tag = element_text(size = 14, face = "bold")
  ) &
  scale_color_gradientn(colours = colorRampPalette(c("white", "#276D9F"))(100)) &
  labs(tag = "Marker Genes") &
  guides(
    size = guide_legend(
      title = "Proportion of\nExpressing Cells",
      title.position = "top",
      title.hjust = 0.5
    ),
    colour = guide_colorbar(
      title = "Average\nExpression",
      title.position = "top",
      title.hjust = 0.5
    )
  )
ggsave("celltype_marker_dotplot.pdf", width = 10, height = 6, dpi = 300)
