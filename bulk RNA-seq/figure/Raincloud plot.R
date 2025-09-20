# Load required libraries
library(ggplot2)     # Core plotting
library(ggsignif)    # Significance annotations
library(ggdist)      # For half-eye (raincloud) effect
library(gghalves)    # Optional half-violin for aesthetics

# Example (for structure)
# data <- data.frame(Group = c(...), Score = c(...))

# Define group order (optional)
group_order <- unique(as.character(data$Group))

# Generate all valid pairwise group comparisons
comparison_list <- list()
for (i in 1:(length(group_order) - 1)) {
  for (j in (i + 1):length(group_order)) {
    g1 <- group_order[i]
    g2 <- group_order[j]
    n1 <- sum(data$Group == g1)
    n2 <- sum(data$Group == g2)
    if (n1 > 1 && n2 > 1) {
      comparison_list[[length(comparison_list) + 1]] <- c(g1, g2)
    }
  }
}

# Define custom color palette (edit to match your groups)
custom_colors <- c("#B7E3ED", "#D87E8F", "#C3C2CD", "#D6BEA4", "#D2E0A0")

# Create raincloud plot
p <- ggplot(data, aes(x = Group, y = Score, fill = Group)) +

  # Add jittered points
  geom_jitter(aes(color = Group), width = 0.2, alpha = 0.6, size = 1) +

  # Add boxplot
  geom_boxplot(
    position = position_nudge(x = 0.2),
    width = 0.15,
    outlier.shape = NA
  ) +

  # Add half-eye (raincloud) density
  stat_halfeye(
    width = 0.4,
    .width = 0,
    justification = -1.2,
    point_colour = NA,
    alpha = 0.5
  ) +

  # Set custom colors
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +

  # Set axis limits
  expand_limits(x = c(0.5, length(unique(data$Group)) + 0.5)) +
  coord_cartesian(ylim = c(0, 0.008)) +

  # Axis labels
  xlab("Group") +
  ylab("Score") +

  # Styling
  theme_minimal() +
  theme(
    axis.ticks.x     = element_line(size = 0.5, color = "black"),
    panel.grid.major = element_line(color = "gray90", size = 0.5),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "black", fill = NA, size = 1),
    legend.position  = "none",
    axis.title.x     = element_text(size = 14),
    axis.title.y     = element_text(size = 14),
    axis.text.x      = element_text(size = 12, angle = 45, hjust = 1),
    axis.text.y      = element_text(size = 12),
    plot.title       = element_text(hjust = 0.5, size = 16)
  ) +

  # Add significance comparisons
  geom_signif(
    comparisons = comparison_list,
    step_increase = 0.1,
    map_signif_level = TRUE,
    vjust = 0.5,
    hjust = 0
  )

# Save plot
ggsave("Score_raincloud.pdf", p, width = 6, height = 3, dpi = 600)
