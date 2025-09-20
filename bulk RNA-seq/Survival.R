# === Load Required Libraries ===
library(survival)
library(survminer)
library(tibble)
library(ggplot2)
library(ggpmisc)  # for geom_table()

# === Step 1: Create Grouping Variables ===
# Split by ISG15 signature and Tex signature using median cutoff
data$ISG15_group <- ifelse(data$CD8T_ISG15_Signature > median(data$CD8T_ISG15_Signature), "high", "low")
data$Tex_group   <- ifelse(data$CD8Tex_Signature > median(data$CD8Tex_Signature), "high", "low")

# Combine into a 4-group classification
data$combined_group <- paste0("CD8T-ISG15(", data$ISG15_group, ")_CD8Tex(", data$Tex_group, ")")

# Set group level order
group_levels <- c(
  "CD8T-ISG15(low)_CD8Tex(low)",
  "CD8T-ISG15(low)_CD8Tex(high)",
  "CD8T-ISG15(high)_CD8Tex(low)",
  "CD8T-ISG15(high)_CD8Tex(high)"
)
group_short <- c("low_low", "low_high", "high_low", "high_high")
data$combined_group <- factor(data$combined_group, levels = group_levels)

# Legend labels for plotting
legend_labels <- c(
  "CD8T-ISG15(low) & CD8Tex(low)",
  "CD8T-ISG15(low) & CD8Tex(high)",
  "CD8T-ISG15(high) & CD8Tex(low)",
  "CD8T-ISG15(high) & CD8Tex(high)"
)

# === Step 2: Fit Kaplan-Meier Survival Model ===
fit <- survfit(Surv(OS.time, OS) ~ combined_group, data = data)

# Log-rank test for global difference
fitd <- survdiff(Surv(OS.time, OS) ~ combined_group, data = data, na.action = na.exclude)
p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
p.label <- paste0("log-rank test P", ifelse(p.val < 0.001, " < 0.001", paste0(" = ", round(p.val, 3))))

# === Step 3: Pairwise Log-rank Test ===
ps <- pairwise_survdiff(Surv(OS.time, OS) ~ combined_group, data = data, p.adjust.method = "none")

# Rename rows/cols with short names for compact display
rownames(ps$p.value) <- group_short[-1]
colnames(ps$p.value) <- group_short[-length(group_short)]

# Format p-value table
p_table <- as.data.frame(ifelse(round(ps$p.value, 3) < 0.001, "<0.001", round(ps$p.value, 3)))
p_table[is.na(p_table)] <- "-"
df_star <- tibble(x = 0, y = 0, tb = list(p_table))  # position for table placement

# === Step 4: Custom Theme and Color Palette ===
theme_custom <- theme(
  text = element_text(family = "sans", size = 12),
  axis.title = element_text(size = 14, face = "bold"),
  axis.text = element_text(size = 12),
  legend.position = "top",
  legend.text = element_text(size = 12),
  legend.title = element_text(size = 13, face = "bold"),
  panel.grid = element_blank(),
  panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
  plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
  plot.margin = margin(10, 10, 10, 10)
)

palette_custom <- c("#8AB6D6", "#C7AED5", "#ACD48A", "#F5BC6E")  # 4-group colors

# === Step 5: Plot KM Curve ===
p <- ggsurvplot(
  fit = fit,
  data = data,
  conf.int = FALSE,
  risk.table = TRUE,
  risk.table.col = "strata",
  palette = palette_custom,
  xlim = c(0, 120),
  break.time.by = 12,
  legend.title = "",
  legend.labs = legend_labels,
  xlab = "Time (months)",
  ylab = "Overall survival",
  risk.table.y.text = FALSE,
  tables.height = 0.3
)

# Add p-value and pairwise p-value table to the plot
p$plot <- p$plot + theme_custom +
  annotate("text", x = 0, y = 0.55, hjust = 0, fontface = 4, label = p.label) +
  geom_table(data = df_star, aes(x = x, y = y, label = tb), table.rownames = TRUE)

# === Step 6: Print Final Plot ===
print(p)
