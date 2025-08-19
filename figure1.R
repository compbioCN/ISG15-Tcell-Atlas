library(Seurat)
library(ggplot2)
library(grid)
colors <- c(
  "#726BAE", "#C7AED5", "#F7DBF0", "#EA945A", "#60A897", "#F5BC6E", "#86A667", 
  "#276D9F", "#6488B9", "#ACD48A", "#A5AA99", "#E5948E", "#E05F48", "#5D88BF", 
  "#78C4D4", "#8AB6D6"
)
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

plot_custom_umap <- function(seurat_obj, meta_col = "celltype", reduction = "umap") {

  umap_coords <- Embeddings(seurat_obj, reduction = reduction)
  meta_data <- seurat_obj@meta.data
  umap_data <- data.frame(
    UMAP1 = umap_coords[, 1],
    UMAP2 = umap_coords[, 2],
    CellType = meta_data[[meta_col]]
  )
  
  p <- ggplot(umap_data, aes(x = UMAP1, y = UMAP2, color = CellType)) +
    geom_point(size = 0.0001, alpha = 0.04) +  
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
  label_positions <- aggregate(cbind(UMAP1, UMAP2) ~ CellType, data = umap_data, FUN = mean)
  p <- p + geom_text(
    data = label_positions, 
    aes(x = UMAP1, y = UMAP2, label = CellType), 
    color = "black", size = 4, fontface = "bold"
  )
  
  cell_count <- nrow(umap_data)
  p <- p + annotate("text", x = min(umap_data$UMAP1), y = max(umap_data$UMAP2), 
                    label = paste0("k=", format(cell_count, big.mark = ","), " cells"),
                    hjust = -0.1, vjust = 1.5, size = 4, fontface = "bold")
  
  # 添加 UMAP 箭头图标（左下角）
  arrow_base_x <- min(umap_data$UMAP1) + 1
  arrow_base_y <- min(umap_data$UMAP2) + 1
  p <- p + 
    geom_segment(aes(x = arrow_base_x, y = arrow_base_y, xend = arrow_base_x + 2, yend = arrow_base_y), 
                 arrow = arrow(length = unit(0.3, "cm")), size = 0.6, color = "black") +  # 横向箭头
    geom_segment(aes(x = arrow_base_x, y = arrow_base_y, xend = arrow_base_x, yend = arrow_base_y + 2), 
                 arrow = arrow(length = unit(0.3, "cm")), size = 0.6, color = "black") +  # 纵向箭头
    annotate("text", x = arrow_base_x + 2.5, y = arrow_base_y, label = "UMAP1", 
             hjust = 0, size = 4, fontface = "bold") +  # 横向箭头文本
    annotate("text", x = arrow_base_x, y = arrow_base_y + 2.5, label = "UMAP2", 
             vjust = 0, size = 4, fontface = "bold")  # 纵向箭头文本
  
  return(p)
}


umap_plot <- plot_custom_umap(OSCC, meta_col = "celltype")
print(umap_plot)


colors <- c(
  "Normal" = "#8AB6D6", 
  "Precancerous" = "#F5BC6E", 
  "Tumor" = "#E05F48"
)
umap_plot <- plot_custom_umap(OSCC, meta_col = "source")
print(umap_plot)


colorlist <- c("#ea5c6f", "#f7905a", "#e187cb", "#fb948d", "#fe9e37", 
               "#e2b159", "#ebed6f", "#b2db87", "#7ee7bb", "#64cccf", 
               "#a9dce6", "#a48cbe", "#e4b7d6", "#6a3d9a")

plot_umap <- function(seurat_obj, meta_col, reduction = "umap", colors, legend = FALSE) {
  umap_coords <- Embeddings(seurat_obj, reduction = reduction)
  meta_data <- seurat_obj@meta.data
  
  umap_data <- data.frame(
    UMAP1 = umap_coords[, 1],
    UMAP2 = umap_coords[, 2],
    Annotation = meta_data[[meta_col]]
  )
  
  p <- ggplot(umap_data, aes(x = UMAP1, y = UMAP2, color = Annotation)) +
    geom_point(size = 0.001, alpha = 0.04) +  # 增大点的大小
    scale_color_manual(values = colors) +
    theme_minimal() +
    theme(
      legend.position = ifelse(legend, "bottom", "none"),  # 图例显示在底部或隐藏
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      panel.background = element_blank(),
      legend.text = element_text(size = 12), 
      legend.title = element_text(size = 14)  
    )
  
  return(p)
}

# 假设 meta_col 是分组注释列
p1 <- plot_umap(OSCC, meta_col = "source", colors = colors, legend = FALSE)

library(ggplot2)

# 定义颜色映射
colors <- c(
  "Normal" = "#8AB6D6", 
  "Precancerous" = "#F5BC6E", 
  "Tumor" = "#E05F48"
)

# 创建空数据框，仅用于生成图例
legend_data <- data.frame(
  Category = names(colors),
  Value = 1:length(colors)  # 假数据，仅用于绘制
)

# 绘制图例
legend_plot <- ggplot(legend_data, aes(x = Value, y = Value, color = Category)) +
  geom_point(size = 4) +  # 图例中的点大小
  scale_color_manual(values = colors) +
  theme_void() +  # 去除坐标轴和网格
  theme(
    legend.position = "bottom",  # 图例放在底部
    legend.text = element_text(size = 12),  # 图例文本大小
    legend.title = element_blank()  # 移除图例标题
  )

# 显示图例
print(legend_plot)

library(ggplot2)

# 自定义颜色映射
colorlist <- c("#ea5c6f", "#f7905a", "#e187cb", "#fb948d", "#fe9e37", 
               "#e2b159", "#ebed6f", "#b2db87", "#7ee7bb", "#64cccf", 
               "#a9dce6", "#a48cbe", "#e4b7d6", "#6a3d9a")

# 映射到分组 P01 到 P13
group_names <- paste0("P", sprintf("%02d", 1:13))  # 确保分组名称为 P01, P02, ..., P13
colors <- setNames(colorlist[1:13], group_names)

# 创建空数据框，仅用于生成图例
legend_data <- data.frame(
  Category = names(colors),
  Value = 1:length(colors)  # 假数据，仅用于绘制
)

# 按 P01 到 P13 顺序绘制图例
legend_plot <- ggplot(legend_data, aes(x = Value, y = Value, color = Category)) +
  geom_point(size = 4) +  # 图例中的点大小
  scale_color_manual(values = colors, breaks = names(colors)) +  # 确保按顺序显示
  guides(color = guide_legend(
    ncol = 4, byrow = TRUE, override.aes = list(size = 5)  # 每列 4 个，点大小调整为 5
  )) +
  theme_void() +  # 去除坐标轴和网格
  theme(
    legend.position = "bottom",  # 图例放在底部
    legend.text = element_text(size = 10),  # 图例文本大小
    legend.title = element_blank()  # 移除图例标题
  )

# 显示图例
print(legend_plot)



#marker点图
A <- FindAllMarkersMAESTRO(OSCC,min.pct = 0.25,logfc.threshold = 0.25)
top10 <-A %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)
features <- list(
  "B cells" = c("MS4A1", "BANK1", "TNFRSF13C"),
  "Endothelial cells" = c("FLT1", "SELE", "TM4SF1"),
  "Epithelial cells" = c("KRT14", "KRT13", "KRT5"),
  "Fibroblasts" = c("DCN", "COL1A1", "LUM"),
  "Mast cells" = c("TPSB2", "CPA3", "TPSAB1"),
  "Mural cells" = c("TAGLN", "ACTA2", "RGS5"),
  "Myeloid cells" = c("LYZ", "IL1B", "C1QA"),
  "NK" = c("GNLY", "NKG7", "GZMB"),
  "Plasma cells" = c("IGLC2", "IGKC", "IGHG1"),
  "Proliferative T cells" = c("MKI67", "STMN1", "TOP2A"),
  "T" = c("CD8A", "CD8B","GZMK")
)

DotPlot(object = OSCC, features = features) &
  theme_bw(base_size = 12) &  # 设置白色背景，基础字体大小
  geom_point(shape = 21, aes(size = pct.exp), stroke = 1.2) &  # 设置点的样式
  theme(
    axis.title = element_blank(),  # 移除坐标轴标题
    axis.text.x = element_text(color = 'black', angle = 45, hjust = 1, vjust = 1, face = "bold", size = 10),  # X轴标签倾斜45度
    axis.text.y = element_text(color = 'black', face = "bold", size = 10),  # Y轴标签样式
    panel.grid.major = element_blank(),  # 移除主网格线
    panel.grid.minor = element_blank(),  # 移除次网格线
    strip.background = element_blank(),  # 移除分面背景
    strip.text = element_text(face = "bold", size = 10),  # 设置分面标签样式
    plot.margin = unit(c(1, 1, 1, 1), 'cm'),  # 设置图形边距
    panel.border = element_rect(color = "black", size = 1.2, linetype = "solid"),  # 设置边框样式
    panel.spacing = unit(0.2, "cm"),  # 设置分面之间的间距
    legend.box.background = element_rect(colour = "black", size = 0.5),  # 图例背景边框
    legend.key.width = unit(0.4, "cm"),  # 图例键宽度
    legend.key.height = unit(0.6, "cm"),  # 图例键高度
    legend.title = element_text(color = 'black', face = "bold", size = 11),  # 图例标题样式
    legend.text = element_text(size = 10)  # 图例文字大小
  ) &
  scale_color_gradientn(colours = colorRampPalette(c("white", "#276D9F"))(100)) &  # 设置颜色梯度
  labs(tag = "Marker Genes") &  # 添加标签
  theme(
    plot.tag.position = c(0.3, 1.1),  # 设置标签位置
    plot.tag = element_text(size = 14, face = "bold")  # 标签样式
  ) &
  guides(
    size = guide_legend(title = "Proportion of\nExpressing Cells", title.position = "top", title.hjust = 0.5),  # 调整比例图例标题
    colour = guide_colorbar(title = "Average\nExpression", title.position = "top", title.hjust = 0.5)  # 调整颜色图例标题
  )


#细胞比例柱状图
source("/home/guile/Code/Singlecellratio_plotstat.R")
my_comparisons <- list(c("Normal", "Precancerous"),c("Normal","Tumor"),c("Precancerous","Tumor"))
Singlecellratio_plotstat <- function(seu, by = "cell.type", meta.include = NULL, 
                                     group_by = NULL, shape_by = NULL,
                                     custom_fill_colors = c("#8AB6D6", "#F5BC6E", "#E05F48"), 
                                     group_by.point = NULL, color_by = NULL, 
                                     pb = FALSE, comparisons = my_comparisons, 
                                     ncol = NULL, label = c("p.format", "p.signif"), 
                                     label.x = NA, pt.size = 4) {
  
  by <- match.arg(by)
  if (is.null(group_by)) {
    group_by <- "null.group"
  }
  
  shapes <- if (!is.null(shape_by)) c(16, 15, 3, 7, 8, 18, 5, 6, 2, 4, 1, 17) else NULL
  
  fq <- prop.table(table(seu@meta.data$celltype, seu@meta.data[,"orig.ident"]), margin = 2) * 100
  df <- reshape2::melt(fq, value.name = "freq", varnames = c("cell.type", "orig.ident"))
  
  uniques <- apply(seu@meta.data, 2, function(x) length(unique(x)))
  ei <- unique(seu@meta.data[, names(uniques[uniques <= 100])])
  ei <- unique(ei[, colnames(ei) %in% meta.include])
  df <- merge(df, ei, by = "orig.ident")
  df <- cbind(df, null.group = paste("1"))
  df$orig.ident <- as.factor(df$orig.ident)
  
  if (is.null(x = ncol)) {
    ncol <- ifelse(length(unique(df$celltype)) > 20, 5, ifelse(length(unique(df$celltype)) > 9, 4, 3))
  }
  
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
  
  if (by == "cell.type" && color_by == "cell.type") {
    p + 
      facet_wrap(group_by, scales = "free_x") +
      geom_bar(aes_string(x = "orig.ident", fill = "factor(cell.type)"), position = "fill", stat = "identity") +
      scale_fill_manual("Cell Type", values = custom_fill_colors) +
      scale_y_continuous(expand = c(0, 0), labels = seq(0, 100, 25)) +
      theme(panel.border = element_blank())
  } else {
    switch(
      by,
      cell.type = p + 
        facet_wrap("cell.type", scales = "free_y", ncol = ncol) + 
        guides(fill = FALSE) +
        geom_boxplot(aes_string(x = group_by), alpha = 0.25, outlier.color = NA) +
        geom_point(size = pt.size, position = position_jitter(width = 0.25), 
                   aes_string(x = group_by, y = "freq", color = color_by, shape = shape_by)) +
        scale_shape_manual(values = shapes) +
        theme(panel.grid.major = element_line(color = "grey", size = 0.25)) +
        scale_color_manual(values = custom_fill_colors) +
        scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
        ggpubr::stat_compare_means(mapping = aes_string(group_by), comparisons = comparisons, label = label, method = "t.test")
    )
  }
}
Singlecellratio_plotstat(
  OSCC,
  group_by = "source",
  meta.include = c("source", "orig.ident"),
  comparisons = my_comparisons,
  color_by = 'source',
  group_by.point = "orig.ident",
  label.x = 1,
  pt.size = 3,
  label = 'p.format',
  ncol = 4,
  custom_fill_colors = c("#8AB6D6", "#F5BC6E", "#E05F48")
)

###CD8umap

colors <- c(
  "CD8Tex" = "#EA945A", 
  "Endothelial cells" = "#86A667", 
  "CD8T-ISG15" = "#E05F48",
  "Fibroblasts" =  "#60A897",
  "CD8Trm" =  "#6488B9",
  "CD8Teff" = "#ACD48A",
  "MAIT" = "#8AB6D6", 
  "CD8Tcm" =  "#726BAE",
  "Plasma cells" =  "#F5BC6E",
  "Proliferative T cells" = "#E5948E",
  "CD8Tem" = "#C7AED5"
)

# 定义 UMAP 图函数
plot_custom_umap <- function(seurat_obj, meta_col = "celltype", reduction = "umap") {
  # 提取 UMAP 坐标和细胞类型
  umap_coords <- Embeddings(seurat_obj, reduction = reduction)
  meta_data <- seurat_obj@meta.data
  
  # 创建绘图数据框
  umap_data <- data.frame(
    UMAP1 = umap_coords[, 1],
    UMAP2 = umap_coords[, 2],
    CellType = meta_data[[meta_col]]
  )
  
  # 绘制 UMAP 图
  p <- ggplot(umap_data, aes(x = UMAP1, y = UMAP2, color = CellType)) +
    geom_point(size = 0.015, alpha = 0.35) +  # 调整点的大小和透明度
    scale_color_manual(values = colors) +  # 应用自定义颜色
    theme_minimal() +  # 去除多余背景
    theme(
      legend.position = "none",  # 隐藏图例
      axis.title = element_blank(),  # 去除坐标轴标题
      axis.text = element_blank(),  # 去除坐标轴刻度
      axis.ticks = element_blank(),  # 去除坐标轴刻度线
      panel.grid = element_blank(),  # 去除网格线
      panel.background = element_blank()  # 设置背景为空白
    )
  
  # 添加标签
  label_positions <- aggregate(cbind(UMAP1, UMAP2) ~ CellType, data = umap_data, FUN = mean)
  p <- p + geom_text(
    data = label_positions, 
    aes(x = UMAP1, y = UMAP2, label = CellType), 
    color = "black", size = 4, fontface = "bold"  # 设置标签样式
  )
  
  # 添加左上角细胞总数目信息
  cell_count <- nrow(umap_data)
  p <- p + annotate("text", x = min(umap_data$UMAP1), y = max(umap_data$UMAP2), 
                    label = paste0("k=", format(cell_count, big.mark = ","), " cells"),
                    hjust = -0.1, vjust = 1.5, size = 4, fontface = "bold")
  
  # 添加 UMAP 箭头图标（左下角）
  arrow_base_x <- min(umap_data$UMAP1) 
  arrow_base_y <- min(umap_data$UMAP2) 
  p <- p + 
    geom_segment(aes(x = arrow_base_x, y = arrow_base_y, xend = arrow_base_x + 2, yend = arrow_base_y), 
                 arrow = arrow(length = unit(0.3, "cm")), size = 0.6, color = "black") +  # 横向箭头
    geom_segment(aes(x = arrow_base_x, y = arrow_base_y, xend = arrow_base_x, yend = arrow_base_y + 2), 
                 arrow = arrow(length = unit(0.3, "cm")), size = 0.6, color = "black") +  # 纵向箭头
    annotate("text", x = arrow_base_x + 2.5, y = arrow_base_y, label = "UMAP1", 
             hjust = 0, size = 4, fontface = "bold") +  # 横向箭头文本
    annotate("text", x = arrow_base_x, y = arrow_base_y + 2.5, label = "UMAP2", 
             vjust = 0, size = 4, fontface = "bold")  # 纵向箭头文本
  
  return(p)
}

# 使用示例
# 假设 Seurat 对象为 OSCC，注释列为 "celltype"
umap_plot <- plot_custom_umap(CD8T, meta_col = "celltype")
print(umap_plot)

###CD4umap
colors <- c(
  "CD4Tstr" = "#EA945A", 
  "Endothelial cells" = "#86A667", 
  "CD4T-ISG15" = "#E05F48",
  "Fibroblasts" =  "#60A897",
  "CD8Trm" =  "#6488B9",
  "CD4Tcm" = "#ACD48A",
  "Th1" = "#8AB6D6", 
  "CD8Tcm" =  "#726BAE",
  "Plasma cells" =  "#F5BC6E",
  "Proliferative T cells" = "#E5948E",
  "CD4Treg" = "#C7AED5"
)
umap_plot <- plot_custom_umap(CD4T, meta_col = "celltype")
print(umap_plot)

#CD8细胞比例柱状图
colors <- c(
  "CD8Tex" = "#EA945A", 
  "Endothelial cells" = "#86A667", 
  "CD8T-ISG15" = "#E05F48",
  "Fibroblasts" = "#60A897",
  "CD8Trm" = "#6488B9",
  "CD8Teff" = "#ACD48A",
  "MAIT" = "#8AB6D6", 
  "CD8Tcm" = "#726BAE",
  "Plasma cells" = "#F5BC6E",
  "Proliferative T cells" = "#E5948E",
  "CD8Tem" = "#C7AED5"
)

Cellratio <- prop.table(table(Idents(CD8T), CD8T$source), margin = 2)
Cellratio <- as.data.frame(Cellratio)
colnames(Cellratio) <- c("celltype", "group", "ratio")
ggplot(Cellratio,aes(x=group,y=ratio,fill=celltype,stratum=celltype,alluvium=celltype))+
geom_col(width = 0.4,color=NA)+
  geom_flow(width=0.4,alpha=0.4,knot.pos=0)+ # knot.pos参数可以使连线变直
  scale_fill_manual(values=colors)+
  theme_classic()+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))+
  theme(legend.position = 'none')
print(p)

#CD4细胞比例柱状图
colors <- c(
  "CD4Tstr" = "#EA945A", 
  "Endothelial cells" = "#86A667", 
  "CD4T-ISG15" = "#E05F48",
  "Fibroblasts" =  "#60A897",
  "CD8Trm" =  "#6488B9",
  "CD4Tcm" = "#ACD48A",
  "Th1" = "#8AB6D6", 
  "CD8Tcm" =  "#726BAE",
  "Plasma cells" =  "#F5BC6E",
  "Proliferative T cells" = "#E5948E",
  "CD4Treg" = "#C7AED5"
)
Cellratio <- prop.table(table(Idents(CD4T), CD4T$source), margin = 2)
Cellratio <- as.data.frame(Cellratio)
colnames(Cellratio) <- c("celltype", "group", "ratio")
ggplot(Cellratio,aes(x=group,y=ratio,fill=celltype,stratum=celltype,alluvium=celltype))+
  geom_col(width = 0.4,color=NA)+
  geom_flow(width=0.4,alpha=0.4,knot.pos=0)+ # knot.pos参数可以使连线变直
  scale_fill_manual(values=colors)+
  theme_classic()+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))+
  theme(legend.position = 'none')
print(p)

###ISG15分子表达小提琴图
library(ggpubr)
library(ggimage)
library(ggplot2)
library(Seurat)
ks_VlnExp <- function(object,
                      group,
                      group_order,
                      features, #仅支持单个基因
                      comparisons,
                      cols = NULL,
                      pieSize = NULL) {
  
  # Seurat 对象的元数据
  meta <- object@meta.data
  if(group %in% colnames(meta)){
    Idents(object) <- group
    Idents(object) <- factor(Idents(object), levels = group_order)
  } else {
    stop("Your group name does not exist. Please provide a correct name in the colnames of metadata.")
  }
  
  # 获取基因表达数据
  exp <- FetchData(object, vars = features)
  colnames(exp) <- 'gene'
  exp$Group <- object@meta.data[,group]
  
  features_exp <- list()
  for (i in seq_along(group_order)) {
    exp1 <- subset(exp, Group == group_order[i])
    colnames(exp1)[1] <- "gene"
    features_exp[[i]] <- exp1
    names(features_exp)[i] <- group_order[i]
  }
  
  # 计算每组中基因表达的百分比
  pct_list <- lapply(features_exp, function(df) {
    pct <- sum(df$gene > 0) / nrow(df)
    return(as.data.frame(pct))
  })
  exp_pct <- do.call(rbind, pct_list)
  colnames(exp_pct) <- "exp_gene"
  exp_pct$non_exp <- 1 - exp_pct$exp_gene
  
  # 定义饼图的颜色
  if(is.null(cols)) {
    colors_map <- c("#766DA7","#C2ADCE","#F98400", "#5BBCD6",'#7F3C8D' ,'#11A579', '#3969AC','#E73F74')
    cols <- colors_map[1:length(group_order)]
  }
  
  # 绘制饼图
  plot_pie <- function(i) {
    df1 <- pivot_longer(exp_pct[i,], cols = everything(), names_to = "type", values_to = "value")
    df1$labels <- scales::percent(df1$value, accuracy = 0.1)
    df1$labels[2] <- ''
    
    ggplot(df1, aes(x= '', y = value, fill = type)) +
      geom_col(color = 'black') + # 饼图边框设置为黑色
      coord_polar(theta = 'y') +
      theme_void() + 
      theme(legend.position = "none") +
      geom_text(aes(label = labels), position = position_stack(vjust = 0.5), size = 4) +
      scale_fill_manual(values = c(cols[i], "grey80")) # 修改饼图填充颜色
  }
  
  exp_pct$pie <- lapply(seq_along(group_order), plot_pie)
  
  exp_pct$x <- seq_along(group_order)
  exp_pct$y <- max(exp$gene, na.rm = TRUE) + 2
  
  if(is.null(pieSize)) {
    pieSize <- 1.5
  }
  
  exp_pct$width <- pieSize
  exp_pct$height <- pieSize
  
  # 绘制小提琴图并添加箱线图
  label.y.pos <- seq(max(exp$gene, na.rm = TRUE), max(exp$gene, na.rm = TRUE) + 100, by = 0.3)
  
  p <- VlnPlot(object, features = features, pt.size = 0) &  # 隐藏散点图
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) &  # 添加箱线图
    theme_bw() &
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_text(color = 'black', face = "bold", size = 12),
      axis.text.y = element_text(color = 'black', face = "bold"),
      axis.title.y = element_text(color = 'black', face = "bold", size = 15),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", size = 1.2, linetype = "solid"),
      panel.spacing = unit(0.12, "cm"),
      plot.title = element_text(hjust = 0.5, face = "bold.italic"),
      legend.position = 'none'
    ) &
    stat_compare_means(
      method = "t.test",
      hide.ns = FALSE,
      comparisons = comparisons,
      label = "p.signif",
      bracket.size = 0.8,
      tip.length = 0,
      size = 5,
      vjust = 0.6,
      label.y = label.y.pos
    ) &
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) &
    scale_fill_manual(values = cols) &
    geom_subview(aes(x = x, y = y, subview = pie, width = width, height = height), data = exp_pct) &
    ylim(0, max(exp$gene, na.rm = TRUE) + 3)
  
  return(p)
}
ks_VlnExp(object = CD4T, group="source",group_order=c("Normal","Precancerous","Tumor"),
          features="ISG15",comparisons=list(c("Normal","Precancerous"),c("Normal","Tumor"),c("Precancerous","Tumor")),cols = c(c("#8AB6D6", "#F5BC6E", "#E05F48")))

###基因评分热图
geneset_list <- lapply(genesetCD8, function(x) na.omit(x))
for (set_name in names(geneset_list)) {
  gene_set <- geneset_list[[set_name]]
  CD8T <- AddModuleScore(
    object = CD8T,
    features = list(gene_set),
    name = set_name,assay = "RNA"
  )
}
meta_columns <- colnames(CD8T@meta.data)
num_last <- 19
old_names <- meta_columns[(length(meta_columns) - num_last + 1):length(meta_columns)]  # 取最后19列
new_names <- gsub("1$", "", old_names)  # 去掉后缀 "1"
colnames(CD8T@meta.data)[(length(meta_columns) - num_last + 1):length(meta_columns)] <- new_names
tail(colnames(CD8T@meta.data), num_last)

geneset_list <- lapply(genesetCD4, function(x) na.omit(x))
for (set_name in names(geneset_list)) {
  gene_set <- geneset_list[[set_name]]
  CD4T <- AddModuleScore(
    object = CD4T,
    features = list(gene_set),
    name = set_name,assay = "RNA"
  )
}
meta_columns <- colnames(CD4T@meta.data)
num_last <- 17
old_names <- meta_columns[(length(meta_columns) - num_last + 1):length(meta_columns)]  # 取最后19列
new_names <- gsub("1$", "", old_names)  # 去掉后缀 "1"
colnames(CD4T@meta.data)[(length(meta_columns) - num_last + 1):length(meta_columns)] <- new_names
tail(colnames(CD4T@meta.data), num_last)


data <- CD8T@meta.data[,c(18,26:44)]
library(dplyr)
avg_data <- data %>%
  group_by(celltype) %>%
  summarise(across(everything(), mean, na.rm = TRUE))
# 除 celltype 列外计算每列的均值
avg_data_long <- avg_data %>%
  column_to_rownames("celltype")
FunctionScoreMatrix <- as.data.frame(t(avg_data_long))
FunctionScoreMatrix <- t(apply(FunctionScoreMatrix, 1, rescale, to = c(-1, 1)))
# 定义颜色
col_fun <- colorRamp2(c(-1, 0, 1), c("#ADD8E6", "white", "#C7AED5"))

# 行注释
row_annot <- rowAnnotation(Signature = rep(c("Differentiation", "Function", "Metabolism", "Apoptosis"), 
                                           times = c(3, 11, 3, 2)),
                           col = list(Signature = c(
                             Differentiation = "#F7DBF0",  # 黄色
                             Function = "#F5BC6E",         # 浅蓝色
                             Metabolism = "#E05F48",# 棕褐色
                             Apoptosis = "#8AB6D6"         # 番茄红
                           )),
                           show_annotation_name = FALSE)

"Normal" = "#8AB6D6", 
"Precancerous" = "#F5BC6E", 
"Tumor" = "#E05F48"


# 绘制热图
p1 <- Heatmap(FunctionScoreMatrix,
        name = "Signature Score",
        col = col_fun,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        row_split = c(rep("Differentiation", 3),
                      rep("Function", 11),
                      rep("Metabolism", 3),
                      rep("Apoptosis", 2)),  # 按模块分块
        row_title_gp = gpar(fontsize = 10),
        column_title_gp = gpar(fontsize = 10),
        row_gap = unit(2, "mm"),
        column_gap = unit(2, "mm"),
        top_annotation = NULL,
        left_annotation = row_annot,
        border = TRUE,
        heatmap_legend_param = list(
          title = "Signature Score",
          at = c(-1, 0, 1),
          labels = c("-1", "0", "1"),
          title_gp = gpar(fontsize = 10),
          labels_gp = gpar(fontsize = 8)
        ))


data <- CD4T@meta.data[,c(18,26:42)]
library(dplyr)
avg_data <- data %>%
  group_by(celltype) %>%
  summarise(across(everything(), mean, na.rm = TRUE))  # 除 celltype 列外计算每列的均值
avg_data_long <- avg_data %>%
  column_to_rownames("celltype")
FunctionScoreMatrix <- as.data.frame(t(avg_data_long))
FunctionScoreMatrix <- t(apply(FunctionScoreMatrix, 1, rescale, to = c(-1, 1)))
# 定义颜色
col_fun <- colorRamp2(c(-1, 0, 1), c("#ADD8E6", "white", "#C7AED5"))

# 行注释
row_annot <- rowAnnotation(Signature = rep(c("Differentiation", "Function", "Metabolism", "Apoptosis"), 
                                           times = c(3, 9, 3, 2)),
                           col = list(Signature = c(
                             Differentiation = "#F7DBF0",  # 黄色
                             Function = "#F5BC6E",         # 浅蓝色
                             Metabolism = "#E05F48",# 棕褐色
                             Apoptosis = "#8AB6D6"         # 番茄红
                           )),
                           show_annotation_name = FALSE)

"Normal" = "#8AB6D6", 
"Precancerous" = "#F5BC6E", 
"Tumor" = "#E05F48"


# 绘制热图
p2 <- Heatmap(FunctionScoreMatrix,
              name = "Signature Score",
              col = col_fun,
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              row_split = c(rep("Differentiation", 3),
                            rep("Function", 9),
                            rep("Metabolism", 3),
                            rep("Apoptosis", 2)),  # 按模块分块
              row_title_gp = gpar(fontsize = 10),
              column_title_gp = gpar(fontsize = 10),
              row_gap = unit(2, "mm"),
              column_gap = unit(2, "mm"),
              top_annotation = NULL,
              left_annotation = row_annot,
              border = TRUE,
              heatmap_legend_param = list(
                title = "Signature Score",
                at = c(-1, 0, 1),
                labels = c("-1", "0", "1"),
                title_gp = gpar(fontsize = 10),
                labels_gp = gpar(fontsize = 8)
              ))
library(patchwork)
library(ComplexHeatmap)
p1 / p2
# 包装为 ggplot 对象
p1_gg <- as.ggplot(p1)
p2_gg <- as.ggplot(p2)

# 使用 patchwork 进行上下排列
p1_gg / p2_gg


#GSE181919umap+柱状图+小提琴图
load("/home/guile/CvI/important dataset/GSE181919注释好T.RData")
colors <- c(
  "#726BAE", "#C7AED5", "#F7DBF0", "#EA945A", "#60A897", "#F5BC6E", "#86A667", 
  "#276D9F", "#6488B9", "#ACD48A", "#A5AA99", "#E5948E", "#E05F48", "#5D88BF", 
  "#78C4D4", "#8AB6D6"
)
# 自定义颜色映射
colors <- c(
  "Tprolif" = "#E5948E",       # Proliferative T cells
  "CD4Treg" = "#C7AED5",       # Regulatory T cells (Treg)
  "CD4Tfh" = "#726BAE",        # Follicular helper T cells (Tfh)
  "GSE181919rm" = "#60A897",        # Tissue-resident memory CD8+ T cells (Trm)
  "CD4Tcm" = "#ACD48A",        # Central memory CD4+ T cells (Tcm)
  "GSE181919ex" = "#8AB6D6",        # Exhausted CD8+ T cells (Tex)
  "Th17" = "#86A667",          # Th17 cells
  "GSE181919em" = "#F5BC6E",        # Effector memory CD8+ T cells (Tem)
  "GSE181919eff" = "#5D88BF",       # Effector CD8+ T cells (Teff)
  "CD4Tstr" = "#6488B9",       # Stressed CD4+ T cells
  "GSE181919str" = "#78C4D4",       # Stressed CD8+ T cells
  "GSE181919-ISG15" = "#E05F48",    # CD8+ T cells with ISG15 expression
  "CD4T-ISG15" = "#EA945A"     # CD4+ T cells with ISG15 expression
)


# 定义 UMAP 图函数
plot_custom_umap <- function(seurat_obj, meta_col = "celltype", reduction = "umap") {
  # 提取 UMAP 坐标和细胞类型
  umap_coords <- Embeddings(seurat_obj, reduction = reduction)
  meta_data <- seurat_obj@meta.data
  
  # 创建绘图数据框
  umap_data <- data.frame(
    UMAP1 = umap_coords[, 1],
    UMAP2 = umap_coords[, 2],
    CellType = meta_data[[meta_col]]
  )
  
  # 绘制 UMAP 图
  p <- ggplot(umap_data, aes(x = UMAP1, y = UMAP2, color = CellType)) +
    geom_point(size = 0.4, alpha = 0.4) +  # 调整点的大小和透明度
    scale_color_manual(values = colors) +  # 应用自定义颜色
    theme_minimal() +  # 去除多余背景
    theme(
      legend.position = "none",  # 隐藏图例
      axis.title = element_blank(),  # 去除坐标轴标题
      axis.text = element_blank(),  # 去除坐标轴刻度
      axis.ticks = element_blank(),  # 去除坐标轴刻度线
      panel.grid = element_blank(),  # 去除网格线
      panel.background = element_blank()  # 设置背景为空白
    )
  
  # 添加标签
  label_positions <- aggregate(cbind(UMAP1, UMAP2) ~ CellType, data = umap_data, FUN = mean)
  p <- p + geom_text(
    data = label_positions, 
    aes(x = UMAP1, y = UMAP2, label = CellType), 
    color = "black", size = 4, fontface = "bold"  # 设置标签样式
  )
  
  # 添加左上角细胞总数目信息
  cell_count <- nrow(umap_data)
  p <- p + annotate("text", x = min(umap_data$UMAP1), y = max(umap_data$UMAP2), 
                    label = paste0("k=", format(cell_count, big.mark = ","), " cells"),
                    hjust = -0.1, vjust = 1.5, size = 4, fontface = "bold")
  
  # 添加 UMAP 箭头图标（左下角）
  arrow_base_x <- min(umap_data$UMAP1) + 1
  arrow_base_y <- min(umap_data$UMAP2) + 1
  p <- p + 
    geom_segment(aes(x = arrow_base_x, y = arrow_base_y, xend = arrow_base_x + 2, yend = arrow_base_y), 
                 arrow = arrow(length = unit(0.3, "cm")), size = 0.6, color = "black") +  # 横向箭头
    geom_segment(aes(x = arrow_base_x, y = arrow_base_y, xend = arrow_base_x, yend = arrow_base_y + 2), 
                 arrow = arrow(length = unit(0.3, "cm")), size = 0.6, color = "black") +  # 纵向箭头
    annotate("text", x = arrow_base_x + 2.5, y = arrow_base_y, label = "UMAP1", 
             hjust = 0, size = 4, fontface = "bold") +  # 横向箭头文本
    annotate("text", x = arrow_base_x, y = arrow_base_y + 2.5, label = "UMAP2", 
             vjust = 0, size = 4, fontface = "bold")  # 纵向箭头文本
  
  return(p)
}

# 使用示例
# 假设 Seurat 对象为 OSCC，注释列为 "celltype"
p1 <- plot_custom_umap(GSE181919, meta_col = "celltype")
print(umap_plot)
Cellratio <- prop.table(table(Idents(GSE181919), GSE181919$Tissue), margin = 2)
Cellratio <- as.data.frame(Cellratio)
colnames(Cellratio) <- c("celltype", "group", "ratio")
p2 <- ggplot(Cellratio,aes(x=group,y=ratio,fill=celltype,stratum=celltype,alluvium=celltype))+
  geom_col(width = 0.4,color=NA)+
  geom_flow(width=0.4,alpha=0.4,knot.pos=0)+ # knot.pos参数可以使连线变直
  scale_fill_manual(values=colors)+
  theme_classic()+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))+
  theme(legend.position = 'none')

unique(Cellratio$group)

p3 <- ks_VlnExp(object = GSE181919, group="Tissue",group_order=c("Normal","OLK","HNSCC","LN"),
          features="ISG15",comparisons=list(c("Normal","OLK"),c("OLK","HNSCC"),c("HNSCC","LN")),cols = c(c("#8AB6D6", "#F5BC6E", "#E05F48","#726BAE")))



#HRA001006umap+柱状图+小提琴图
load("/home/guile/CvI/important dataset/HRA001006注释好T.RData")
colors <- c(
  "#726BAE", "#C7AED5", "#F7DBF0", "#EA945A", "#60A897", "#F5BC6E", "#86A667", 
  "#276D9F", "#6488B9", "#ACD48A", "#A5AA99", "#E5948E", "#E05F48", "#5D88BF", 
  "#78C4D4", "#8AB6D6"
)

colors <- c(
  "CD8Tcm" = "#ACD48A",        # Central memory CD8+ T cells
  "CD8Tem" = "#F5BC6E",        # Effector memory CD8+ T cells
  "CD8Trm" = "#60A897",        # Tissue-resident memory CD8+ T cells
  "CD4Tcm" = "#C7AED5",        # Central memory CD4+ T cells
  "Th17" = "#86A667",          # Th17 cells
  "CD4Treg" = "#726BAE",       # Regulatory T cells (Treg)
  "CD8Tex" = "#8AB6D6",        # Exhausted CD8+ T cells
  "CD8Teff" = "#6488B9",       # Effector CD8+ T cells
  "CD4-ISG15" = "#EA945A",     # CD4+ T cells with ISG15 expression
  "CD8Tstr" = "#78C4D4",       # Stressed CD8+ T cells
  "Tfh" = "#5D88BF",           # Follicular helper T cells
  "CD8-ISG15" = "#E05F48",     # CD8+ T cells with ISG15 expression
  "CD8Tn" = "#E5948E",         # Naive CD8+ T cells
  "MAIT" = "#276D9F"           # Mucosal-associated invariant T cells
)

# 定义 UMAP 图函数
plot_custom_umap <- function(seurat_obj, meta_col = "celltype", reduction = "umap") {
  # 提取 UMAP 坐标和细胞类型
  umap_coords <- Embeddings(seurat_obj, reduction = reduction)
  meta_data <- seurat_obj@meta.data
  
  # 创建绘图数据框
  umap_data <- data.frame(
    UMAP1 = umap_coords[, 1],
    UMAP2 = umap_coords[, 2],
    CellType = meta_data[[meta_col]]
  )
  
  # 绘制 UMAP 图
  p <- ggplot(umap_data, aes(x = UMAP1, y = UMAP2, color = CellType)) +
    geom_point(size = 0.1, alpha = 0.25) +  # 调整点的大小和透明度
    scale_color_manual(values = colors) +  # 应用自定义颜色
    theme_minimal() +  # 去除多余背景
    theme(
      legend.position = "none",  # 隐藏图例
      axis.title = element_blank(),  # 去除坐标轴标题
      axis.text = element_blank(),  # 去除坐标轴刻度
      axis.ticks = element_blank(),  # 去除坐标轴刻度线
      panel.grid = element_blank(),  # 去除网格线
      panel.background = element_blank()  # 设置背景为空白
    )
  
  # 添加标签
  label_positions <- aggregate(cbind(UMAP1, UMAP2) ~ CellType, data = umap_data, FUN = mean)
  p <- p + geom_text(
    data = label_positions, 
    aes(x = UMAP1, y = UMAP2, label = CellType), 
    color = "black", size = 4, fontface = "bold"  # 设置标签样式
  )
  
  # 添加左上角细胞总数目信息
  cell_count <- nrow(umap_data)
  p <- p + annotate("text", x = min(umap_data$UMAP1), y = max(umap_data$UMAP2), 
                    label = paste0("k=", format(cell_count, big.mark = ","), " cells"),
                    hjust = -0.1, vjust = 1.5, size = 4, fontface = "bold")
  
  # 添加 UMAP 箭头图标（左下角）
  arrow_base_x <- min(umap_data$UMAP1) + 1
  arrow_base_y <- min(umap_data$UMAP2) + 1
  p <- p + 
    geom_segment(aes(x = arrow_base_x, y = arrow_base_y, xend = arrow_base_x + 2, yend = arrow_base_y), 
                 arrow = arrow(length = unit(0.3, "cm")), size = 0.6, color = "black") +  # 横向箭头
    geom_segment(aes(x = arrow_base_x, y = arrow_base_y, xend = arrow_base_x, yend = arrow_base_y + 2), 
                 arrow = arrow(length = unit(0.3, "cm")), size = 0.6, color = "black") +  # 纵向箭头
    annotate("text", x = arrow_base_x + 2.5, y = arrow_base_y, label = "UMAP1", 
             hjust = 0, size = 4, fontface = "bold") +  # 横向箭头文本
    annotate("text", x = arrow_base_x, y = arrow_base_y + 2.5, label = "UMAP2", 
             vjust = 0, size = 4, fontface = "bold")  # 纵向箭头文本
  
  return(p)
}

# 使用示例
# 假设 Seurat 对象为 OSCC，注释列为 "celltype"
p4 <- plot_custom_umap(NKT, meta_col = "celltype")
print(p4)
Cellratio <- prop.table(table(Idents(NKT), NKT$Tissue), margin = 2)
Cellratio <- as.data.frame(Cellratio)
colnames(Cellratio) <- c("celltype", "group", "ratio")
p5 <- ggplot(Cellratio,aes(x=group,y=ratio,fill=celltype,stratum=celltype,alluvium=celltype))+
  geom_col(width = 0.4,color=NA)+
  geom_flow(width=0.4,alpha=0.4,knot.pos=0)+ # knot.pos参数可以使连线变直
  scale_fill_manual(values=colors)+
  theme_classic()+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))+
  theme(legend.position = 'none')

unique(Cellratio$group)

p6 <- ks_VlnExp(object = NKT, group="Tissue",group_order=c("Normal","OLK","OSCC"),
                features="ISG15",comparisons=list(c("Normal","OLK"),c("OLK","OSCC"),c("Normal","OSCC")),cols = c(c("#8AB6D6", "#F5BC6E", "#E05F48","#726BAE")))

p1|p4

#亚群映射大群
plot_custom_umap <- function(seurat_obj, meta_col = "celltype", reduction = "umap") {
  # 提取 UMAP 坐标和细胞类型
  umap_coords <- Embeddings(seurat_obj, reduction = reduction)
  meta_data <- seurat_obj@meta.data
  
  # 创建绘图数据框
  umap_data <- data.frame(
    UMAP1 = umap_coords[, 1],
    UMAP2 = umap_coords[, 2],
    CellType = meta_data[[meta_col]]
  )
  
  # 绘制 UMAP 图
  p <- ggplot(umap_data, aes(x = UMAP1, y = UMAP2, color = CellType)) +
    geom_point(size = 0.01, alpha = 0.1) +  # 调整点的大小和透明度
    scale_color_manual(values = colors) +  # 应用自定义颜色
    theme_minimal() +  # 去除多余背景
    theme(
      legend.position = "none",  # 隐藏图例
      axis.title = element_blank(),  # 去除坐标轴标题
      axis.text = element_blank(),  # 去除坐标轴刻度
      axis.ticks = element_blank(),  # 去除坐标轴刻度线
      panel.grid = element_blank(),  # 去除网格线
      panel.background = element_blank()  # 设置背景为空白
    )
  
  # 添加标签
  label_positions <- aggregate(cbind(UMAP1, UMAP2) ~ CellType, data = umap_data, FUN = mean)
  p <- p + geom_text(
    data = label_positions, 
    aes(x = UMAP1, y = UMAP2, label = CellType), 
    color = "black", size = 4, fontface = "bold"  # 设置标签样式
  )
  
  # 添加左上角细胞总数目信息
  cell_count <- nrow(umap_data)
  p <- p + annotate("text", x = min(umap_data$UMAP1), y = max(umap_data$UMAP2), 
                    label = paste0("k=", format(cell_count, big.mark = ","), " cells"),
                    hjust = -0.1, vjust = 1.5, size = 4, fontface = "bold")
  
  # 添加 UMAP 箭头图标（左下角）
  arrow_base_x <- min(umap_data$UMAP1) + 1
  arrow_base_y <- min(umap_data$UMAP2) + 1
  p <- p + 
    geom_segment(aes(x = arrow_base_x, y = arrow_base_y, xend = arrow_base_x + 2, yend = arrow_base_y), 
                 arrow = arrow(length = unit(0.3, "cm")), size = 0.6, color = "black") +  # 横向箭头
    geom_segment(aes(x = arrow_base_x, y = arrow_base_y, xend = arrow_base_x, yend = arrow_base_y + 2), 
                 arrow = arrow(length = unit(0.3, "cm")), size = 0.6, color = "black") +  # 纵向箭头
    annotate("text", x = arrow_base_x + 2.5, y = arrow_base_y, label = "UMAP1", 
             hjust = 0, size = 4, fontface = "bold") +  # 横向箭头文本
    annotate("text", x = arrow_base_x, y = arrow_base_y + 2.5, label = "UMAP2", 
             vjust = 0, size = 4, fontface = "bold")  # 纵向箭头文本
  
  return(p)
}
colors <- c(
  "no" = "#5D88BF",        # Central memory CD8+ T cells
  "CD4ISG" = "#EA945A",        # Effector memory CD8+ T cells
  "CD8ISG" = "#E05F48"       # Mucosal-associated invariant T cells
)
p2 <- plot_custom_umap(GSE181919all, meta_col = "CD8ISG")
P3 <- p1|p2
GSE181919all <- seurat_metadata

load("/home/guile/CvI/important dataset/GSE181919注释好T.RData")
CD8ISG15 <- GSE181919[,GSE181919$celltype%in%c("CD8T-ISG15")]
GSE181919all$CD8ISG <- ifelse(colnames(GSE181919all)%in%colnames(CD8ISG15),"CD8ISG","no")

load("/home/guile/CvI/data/HRA001006.Rdata")
HRA001 <- NKT

HRA001006 <- seurat_metadata
CD8ISG15 <- NKT[,NKT$celltype%in%c("CD8-ISG15")]
HRA001006$CD8ISG <- ifelse(colnames(HRA001006)%in%colnames(CD8ISG15),"CD8ISG","no")
p2 <- plot_custom_umap(HRA001006, meta_col = "CD8ISG")
P3 <- p1|p2

P1/P2/P3
P1
P2
P3

#韦恩图
library(ggplot2)
library(ggVennDiagram)
venn_list <- list(
  GSE = top100_CD8GSE,
  HRA = top100_CD8HRA,
  MY  = top100_CD8my
)

p1 <- ggVennDiagram(venn_list, label_alpha = 0, set_size = 5) +
  scale_fill_gradientn(colors = c(
    rgb(1, 1, 1, alpha = 0),      # 白色，完全透明
    rgb(0.36, 0.53, 0.75, alpha = 0.5) # #5D88BF 带透明度
  )) +
  scale_color_manual(values = c("black")) +  # 设置边框颜色
  theme_void()

p2 <- ggVennDiagram(venn_list, label_alpha = 0, set_size = 5) +
  scale_fill_gradientn(colors = c(
    rgb(1, 1, 1, alpha = 0),          # 完全透明的白色
    rgb(0.92, 0.58, 0.35, alpha = 0.5) # #EA945A 带透明度
  )) +
  scale_color_manual(values = c("black")) +  # 设置边框颜色为黑色
  theme_void() +
  theme(
    plot.title = element_text(size = 14, hjust = 0.5, face = "bold")
  )
p1|p2
common_genesCD4 <- Reduce(intersect, list(top100_CD4GSE, top100_CD4HRA, top100_CD4my))
common_genesCD8 <- Reduce(intersect, list(top100_CD8GSE, top100_CD8HRA, top100_CD8my))


#基因列表图
# 加载必要的包
library(ggplot2)

# 定义基因分类
isg_genes <- c("ISG15", "IFIT3", "IFIT2", "IFIT1", "IFI6", "IFI44L", "EPSTI1", "MX1", "OAS1", "OAS3", "OASL", "RSAD2")
t_cell_markers <- c("CD4", "CD3E", "CD3D", "LTB")
signaling_genes <- c("STAT1", "LY6E")

# 创建数据框，包含基因和分类信息
genes_data <- data.frame(
  Gene = c(isg_genes, t_cell_markers, signaling_genes),
  Category = factor(c(
    rep("ISGs", length(isg_genes)),
    rep("T Cell Markers", length(t_cell_markers)),
    rep("Signaling Genes", length(signaling_genes))
  ), levels = c("ISGs", "T Cell Markers", "Signaling Genes"))
)


# 定义分类颜色与透明度
category_colors <- c(
  "ISGs" =  alpha("#E05F48", 0.5),         # 高级蓝色带透明度
  "T Cell Markers" = alpha("#EA945A", 0.5), # 柔和橙色带透明度
  "Signaling Genes" =  alpha("#5D88BF", 0.5) # 质感红色带透明度
)

# 绘制方框图：美化适用于发表
p1 <- ggplot(genes_data, aes(x = 1, y = rev(seq_along(Gene)), fill = Category)) +
  # 添加带阴影的方框（增强质感）
  geom_tile(color = "grey90", size = 0.3, width = 0.92, height = 0.92) +  
  # 居中显示基因名称
  geom_text(aes(label = Gene), size = 4.5, color = "black") +  # 去掉字体加粗
  # 设置颜色和图例
  scale_fill_manual(values = category_colors, name = "Gene Categories") +
  # 设置高级主题
  theme_minimal(base_size = 12) +
  theme(
    panel.background = element_rect(fill = "white", color = NA), # 背景纯白
    plot.background = element_rect(fill = "white", color = NA),  # 绘图区域背景
    axis.title = element_blank(),      # 移除坐标轴标题
    axis.text = element_blank(),       # 移除坐标轴文本
    axis.ticks = element_blank(),      # 移除坐标轴刻度
    panel.grid = element_blank(),      # 移除网格线
    legend.position = "right",         # 图例位置
    legend.title = element_text(size = 12),    # 图例标题样式（不加粗）
    legend.text = element_text(size = 10, color = "grey30"),  # 图例文本颜色与大小
    plot.title = element_text(hjust = 0.5, size = 16, color = "#333333"), # 标题美化（不加粗）
    plot.subtitle = element_text(hjust = 0.5, size = 12, face = "italic", color = "grey40"), # 副标题
    plot.caption = element_text(hjust = 1, size = 10, color = "grey50") # 数据来源说明
  ) 


# 加载必要的包
library(ggplot2)
library(scales)  # 控制颜色和透明度

# 定义基因分类和列表
categories <- c(
  "Interferon-Stimulated Genes (ISGs)",
  "Cytotoxicity-Related Genes",
  "Exhaustion-Associated Genes",
  "T Cell Markers"
)

gene_lists <- list(
  ISGs = c("OASL", "ISG15", "IFIT3", "IFIT1", "IFIT2", "IFI6", "IFI44L", "MX1", "RSAD2", "CMPK2", "OAS1", "EPSTI1", "IRF7", "STAT1", "LY6E", "IFNG"),
  TCellMarkers = c("CD8A", "CD8B", "CD3D", "CD3E"),Cytotoxic = c("GZMK", "GZMA", "PRF1", "CCL5", "NKG7"),
  Exhaustion = c("LAG3", "CD38")
)

# 创建数据框
genes_data <- data.frame(
  Gene = unlist(gene_lists),
  Category = rep(names(gene_lists), times = sapply(gene_lists, length))
)

# 将分类名映射到英文
category_labels <- c(
  ISGs = "Interferon-Stimulated Genes (ISGs)",
  Cytotoxic = "Cytotoxicity-Related Genes",
  Exhaustion = "Exhaustion-Associated Genes",
  TCellMarkers = "T Cell Markers"
)
genes_data$Category <- factor(genes_data$Category, levels = names(category_labels), labels = category_labels)

# 定义分类颜色与透明度
category_colors <- c(
  "Interferon-Stimulated Genes (ISGs)" = alpha("#E05F48", 0.5),
  "Cytotoxicity-Related Genes" = alpha("#726BAE", 0.5),
  "Exhaustion-Associated Genes" = alpha("#ACD48A", 0.5),
  "T Cell Markers" = alpha("#EA945A", 0.5)
)



p2 <- ggplot(genes_data, aes(x = 1, y = rev(seq_along(Gene)), fill = Category)) +
  # 添加带阴影的方框（增强质感）
  geom_tile(color = "grey90", size = 0.3, width = 0.92, height = 0.92) +  
  # 居中显示基因名称
  geom_text(aes(label = Gene), size = 4.5, color = "black") +  # 去掉字体加粗
  # 设置颜色和图例
  scale_fill_manual(values = category_colors, name = "Gene Categories") +
  # 设置高级主题
  theme_minimal(base_size = 12) +
  theme(
    panel.background = element_rect(fill = "white", color = NA), # 背景纯白
    plot.background = element_rect(fill = "white", color = NA),  # 绘图区域背景
    axis.title = element_blank(),      # 移除坐标轴标题
    axis.text = element_blank(),       # 移除坐标轴文本
    axis.ticks = element_blank(),      # 移除坐标轴刻度
    panel.grid = element_blank(),      # 移除网格线
    legend.position = "right",         # 图例位置
    legend.title = element_text(size = 12),    # 图例标题样式（不加粗）
    legend.text = element_text(size = 10, color = "grey30"),  # 图例文本颜色与大小
    plot.title = element_text(hjust = 0.5, size = 16, color = "#333333"), # 标题美化（不加粗）
    plot.subtitle = element_text(hjust = 0.5, size = 12, face = "italic", color = "grey40"), # 副标题
    plot.caption = element_text(hjust = 1, size = 10, color = "grey50") # 数据来源说明
  ) 
# 打印图表
print(p)
p1|p2


### bulk映射
genelist <- list(common_genesCD4,common_genesCD8)
names(genelist) <- c("CD4-ISG15 Signature", "CD8-ISG15 Signature")
TCGAscore <- gsva(expr=as.matrix(logTPM), 
                      gset.idx.list=genelist, 
                      mx.diff=F,
                      kcdf="Gaussian", #CPM, RPKM, TPM数据用"Gaussian"， read count用"Poisson"，
                      parallel.sz=16)
TCGAscore <- as.data.frame(t(TCGAscore))
TCGAscore$group <- ifelse(substr(rownames(TCGAscore),14,16)%in%c("01"),"Tumor","no")
TCGAscore$group <- ifelse(substr(rownames(TCGAscore),14,16)%in%c("11"),"Normol",TCGAscore$group)
TCGAscore <- TCGAscore[TCGAscore$group%in%c("Normol","Tumor"),]

library(ggplot2)
library(ggsignif)
library(gghalves)
library(ggdist)

# 添加 Signature 列，将 CD4 和 CD8 的数据合并
TCGAscore_long <- TCGAscore %>%
  tidyr::pivot_longer(
    cols = c(`CD4-ISG15 Signature`, `CD8-ISG15 Signature`),
    names_to = "Signature",
    values_to = "Score"
  )

# 自定义颜色列表
Custom.color <- c(alpha("#8AB6D6", 0.5),alpha("#F5BC6E", 0.5),alpha("#E05F48", 0.5) )

# 提取唯一的 group 列表
groups <- unique(TCGAscore_long$group)
comb_list <- combn(groups, 2, simplify = FALSE)

# 过滤掉数据点不足的组
valid_combinations <- lapply(comb_list, function(comb) {
  group1 <- TCGAscore_long$Score[TCGAscore_long$group == comb[1]]
  group2 <- TCGAscore_long$Score[TCGAscore_long$group == comb[2]]
  if (length(group1) > 1 && length(group2) > 1) {
    return(comb)
  } else {
    return(NULL)
  }
})
valid_combinations <- Filter(Negate(is.null), valid_combinations)

# 计算每组的样本例数
sample_counts <- TCGAscore_long %>%
  group_by(group, Signature) %>%
  summarise(n = n(), .groups = "drop")

# 绘制分面图并添加样本例数
P <- ggplot(TCGAscore_long, aes(x = group, y = Score, fill = group)) +
  geom_jitter(mapping = aes(color = group), width = 0.2, alpha = 0.6, size = 1) +
  geom_boxplot(position = position_nudge(x = 0.2), width = 0.15, outlier.size = 0, outlier.alpha = 0) +
  stat_halfeye(mapping = aes(fill = group), width = 0.4, .width = 0, justification = -1.2, point_colour = NA, alpha = 0.5) +
  scale_fill_manual(values = Custom.color) +
  scale_color_manual(values = Custom.color) +
  facet_wrap(~Signature, scales = "free_y") +  # 分面显示 CD4 和 CD8
  xlab("Group") +
  ylab("Signature Score") +
  theme_minimal() +
  theme(
    axis.ticks.x = element_line(size = 0.5, color = "black"),
    panel.grid.major = element_line(color = "gray90", size = 0.5),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    legend.position = "none",
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12),
    strip.text = element_text(size = 14),  # 调整分面标签字体大小
    plot.title = element_text(hjust = 0.5, size = 16)
  ) +
  geom_signif(comparisons = valid_combinations, step_increase = 0.1, map_signif_level = FALSE, vjust = 0.5, hjust = 0) +
  # 在横坐标标签下方添加样本例数
  geom_text(data = sample_counts, aes(x = group, y = min(TCGAscore_long$Score) - 0.1 * abs(min(TCGAscore_long$Score)), 
                                      label = paste0("n=", n)),
            inherit.aes = FALSE, size = 4, color = "black", vjust = 1.5)

# 输出图形
print(P)

# 输出图形
print(P)
Ptcga <- P

7:4
GSE39366score <- gsva(expr=as.matrix(expr), 
                  gset.idx.list=genelist, 
                  mx.diff=F,
                  kcdf="Gaussian", #CPM, RPKM, TPM数据用"Gaussian"， read count用"Poisson"，
                  parallel.sz=64)
GSE39366score <- as.data.frame(t(GSE39366score))
GSE39366score <- cbind(GSE39366score,clinical)
library(GEOquery)
gse_id <- "GSE117973"
gse_data <- getGEO(gse_id, GSEMatrix = TRUE)
gse <- gse_data[[1]]
if (length(gset) > 1) idx <- grep("GPL9053", attr(gset, "names")) else idx <- 1
gset <- gse_data[[idx]]

ex <- exprs(gset)
clinical_info <- pData(gse)
clinical_info <- clinical_info[,c(2,45:57)]
GSE26549score <- cbind(GSE26549score,clinical_info)
GSE39366score_long <- GSE39366score%>%
  tidyr::pivot_longer(
    cols = c(`CD4-ISG15 Signature`, `CD8-ISG15 Signature`),
    names_to = "Signature",
    values_to = "Score"
  )

Custom.color <- c(alpha("#8AB6D6", 0.5),alpha("#E05F48", 0.5), "#726BAE" )
GSE39366score_long$group <- GSE39366score_long$condition
# 提取唯一的 group 列表
groups <- unique(GSE39366score_long$group)
comb_list <- combn(groups, 2, simplify = FALSE)

# 过滤掉数据点不足的组
valid_combinations <- lapply(comb_list, function(comb) {
  group1 <- GSE39366score_long$Score[GSE39366score_long$group == comb[1]]
  group2 <- GSE39366score_long$Score[GSE39366score_long$group == comb[2]]
  if (length(group1) > 1 && length(group2) > 1) {
    return(comb)
  } else {
    return(NULL)
  }
})
valid_combinations <- Filter(Negate(is.null), valid_combinations)

# 计算每组的样本例数
sample_counts <- GSE39366score_long %>%
  group_by(group, Signature) %>%
  summarise(n = n(), .groups = "drop")

# 绘制分面图并添加样本例数
P <- ggplot(GSE39366score_long, aes(x = group, y = Score, fill = group)) +
  geom_jitter(mapping = aes(color = group), width = 0.2, alpha = 0.6, size = 1) +
  geom_boxplot(position = position_nudge(x = 0.2), width = 0.15, outlier.size = 0, outlier.alpha = 0) +
  stat_halfeye(mapping = aes(fill = group), width = 0.4, .width = 0, justification = -1.2, point_colour = NA, alpha = 0.5) +
  scale_fill_manual(values = Custom.color) +
  scale_color_manual(values = Custom.color) +
  facet_wrap(~Signature, scales = "free_y") +  # 分面显示 CD4 和 CD8
  xlab("Group") +
  ylab("Signature Score") +
  theme_minimal() +
  theme(
    axis.ticks.x = element_line(size = 0.5, color = "black"),
    panel.grid.major = element_line(color = "gray90", size = 0.5),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    legend.position = "none",
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12),
    strip.text = element_text(size = 14),  # 调整分面标签字体大小
    plot.title = element_text(hjust = 0.5, size = 16)
  ) +
  geom_signif(comparisons = valid_combinations, step_increase = 0.1, map_signif_level = FALSE, vjust = 0.5, hjust = 0) +
  # 在横坐标标签下方添加样本例数
  geom_text(data = sample_counts, aes(x = group, y = min(GSE39366score_long$Score) - 0.1 * abs(min(TCGAscore_long$Score)), 
                                      label = paste0("n=", n)),
            inherit.aes = FALSE, size = 4, color = "black", vjust = 1.5)
final_plot <- plot_grid(
  plot_grid(Ptcga, PGSE42743,PGSE30784, ncol = 3),        # 上排两张
  plot_grid(PGSE13601, PGSE78060, ncol = 3), # 下排三张
  ncol = 1,                                    # 总体分两行
  rel_heights = c(1, 1)                        # 两行高度相等
)
print(final_plot)




