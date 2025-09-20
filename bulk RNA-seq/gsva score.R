# =========================================
# Bulk Signature Scoring & Visualization
# =========================================

library(GSVA)
library(ggplot2)
library(dplyr)
library(ggsignif)
library(gghalves)
library(ggdist)
library(tidyr)
library(cowplot)

# -----------------------------
# Define gene sets
# -----------------------------
genelist <- list(
  "CD4-ISG15 Signature" = common_genesCD4,
  "CD8-ISG15 Signature" = common_genesCD8
)

# -----------------------------
# GSVA: TCGA TPM Matrix
# -----------------------------
TCGAscore <- gsva(expr = as.matrix(logTPM),
                  gset.idx.list = genelist,
                  mx.diff = FALSE,
                  kcdf = "Gaussian",
                  parallel.sz = 16) %>%
  t() %>%
  as.data.frame()

# Annotate tumor vs normal
TCGAscore$group <- "Other"
TCGAscore$group[substr(rownames(TCGAscore), 14, 16) == "01"] <- "Tumor"
TCGAscore$group[substr(rownames(TCGAscore), 14, 16) == "11"] <- "Normal"
TCGAscore <- TCGAscore[TCGAscore$group %in% c("Tumor", "Normal"), ]

# Long format
TCGAscore_long <- TCGAscore %>%
  pivot_longer(cols = starts_with("CD"), names_to = "Signature", values_to = "Score")

Custom.color <- c(alpha("#8AB6D6", 0.5), alpha("#F5BC6E", 0.5), alpha("#E05F48", 0.5))

# Statistical comparisons
groups <- unique(TCGAscore_long$group)
valid_combinations <- combn(groups, 2, simplify = FALSE) %>%
  Filter(function(x) {
    all(sapply(x, function(g) sum(TCGAscore_long$group == g) > 1))
  }, .)

sample_counts <- TCGAscore_long %>%
  group_by(group, Signature) %>%
  summarise(n = n(), .groups = "drop")

Ptcga <- ggplot(TCGAscore_long, aes(x = group, y = Score, fill = group)) +
  geom_jitter(aes(color = group), width = 0.2, alpha = 0.6, size = 1) +
  geom_boxplot(position = position_nudge(x = 0.2), width = 0.15, outlier.alpha = 0) +
  stat_halfeye(width = 0.4, .width = 0, justification = -1.2, point_colour = NA, alpha = 0.5) +
  scale_fill_manual(values = Custom.color) +
  scale_color_manual(values = Custom.color) +
  facet_wrap(~Signature, scales = "free_y") +
  labs(x = "Group", y = "Signature Score") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    panel.border = element_rect(fill = NA, color = "black")
  ) +
  geom_signif(comparisons = valid_combinations, step_increase = 0.1, vjust = 0.5) +
  geom_text(data = sample_counts,
            aes(x = group, y = min(TCGAscore_long$Score) - 0.1 * abs(min(Score)),
                label = paste0("n=", n)),
            inherit.aes = FALSE, size = 4, vjust = 1.5)

# -----------------------------
# GSVA: GSE39366
# -----------------------------
GSE39366score <- gsva(expr = as.matrix(expr),
                      gset.idx.list = genelist,
                      mx.diff = FALSE,
                      kcdf = "Gaussian",
                      parallel.sz = 64) %>%
  t() %>%
  as.data.frame()

GSE39366score$group <- clinical$condition
GSE39366score_long <- GSE39366score %>%
  pivot_longer(cols = starts_with("CD"), names_to = "Signature", values_to = "Score")

Custom.color <- c(alpha("#8AB6D6", 0.5), alpha("#E05F48", 0.5), "#726BAE")

groups <- unique(GSE39366score_long$group)
valid_combinations <- combn(groups, 2, simplify = FALSE) %>%
  Filter(function(x) {
    all(sapply(x, function(g) sum(GSE39366score_long$group == g) > 1))
  }, .)

sample_counts <- GSE39366score_long %>%
  group_by(group, Signature) %>%
  summarise(n = n(), .groups = "drop")

PGSE39366 <- ggplot(GSE39366score_long, aes(x = group, y = Score, fill = group)) +
  geom_jitter(aes(color = group), width = 0.2, alpha = 0.6, size = 1) +
  geom_boxplot(position = position_nudge(x = 0.2), width = 0.15, outlier.alpha = 0) +
  stat_halfeye(width = 0.4, .width = 0, justification = -1.2, point_colour = NA, alpha = 0.5) +
  scale_fill_manual(values = Custom.color) +
  scale_color_manual(values = Custom.color) +
  facet_wrap(~Signature, scales = "free_y") +
  labs(x = "Group", y = "Signature Score") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    panel.border = element_rect(fill = NA, color = "black")
  ) +
  geom_signif(comparisons = valid_combinations, step_increase = 0.1, vjust = 0.5) +
  geom_text(data = sample_counts,
            aes(x = group, y = min(GSE39366score_long$Score) - 0.1 * abs(min(Score)),
                label = paste0("n=", n)),
            inherit.aes = FALSE, size = 4, vjust = 1.5)

# -----------------------------
# Combine All Plots
# -----------------------------
# Add PGSE42743, PGSE30784, PGSE13601, PGSE78060 as defined elsewhere

final_plot <- plot_grid(
  plot_grid(Ptcga, PGSE42743, PGSE30784, ncol = 3),
  plot_grid(PGSE13601, PGSE78060, PGSE39366, ncol = 3),
  ncol = 1, rel_heights = c(1, 1)
)

print(final_plot)
