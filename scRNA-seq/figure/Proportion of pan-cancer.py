# ================================
# Bar plot of CD4⁺ Tisg proportions across tissue sources (no manual override)
# ================================

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick

# Extract relevant metadata
df = adata.obs[["Source_reclassified", "cell_type"]].copy()

# Count total cells and CD4+ Tisg cells per source
total_counts = df.groupby("Source_reclassified", observed=True).size().reset_index(name="total")
tisg_counts = df[df["cell_type"] == "CD4+ Tisg"].groupby("Source_reclassified", observed=True).size().reset_index(name="CD4_Tisg_count")

# Merge and compute proportion
merged = pd.merge(total_counts, tisg_counts, on="Source_reclassified", how="left")
merged["CD4_Tisg_count"] = merged["CD4_Tisg_count"].astype(float).fillna(0)
merged["CD4_Tisg_ratio"] = merged["CD4_Tisg_count"] / merged["total"]

# Sort by ratio for plotting
merged = merged.sort_values("CD4_Tisg_ratio", ascending=True).reset_index(drop=True)

# Plot settings
colors = ["#8AB6D6", "#C7AED5", "#F5BC6E", "#E05F48"]
sns.set_style("whitegrid")
sns.set_context("notebook", font_scale=1.2)

plt.figure(figsize=(5, 5))
bars = plt.bar(
    merged["Source_reclassified"],
    merged["CD4_Tisg_ratio"],
    color=colors,
    edgecolor="black",
    alpha=0.95
)

# Annotate each bar with percentage
for bar, ratio in zip(bars, merged["CD4_Tisg_ratio"]):
    xpos = bar.get_x() + bar.get_width() / 2
    height = bar.get_height()
    label = f"{height:.2%}"

    if height < 0.02:
        plt.text(xpos, height + 0.002, label, ha="center", va="bottom", fontsize=11, fontweight="semibold")
    else:
        plt.text(xpos, height / 2, label, ha="center", va="center", fontsize=11, fontweight="semibold", color="white")

# Format y-axis as percentage
plt.ylim(0, merged["CD4_Tisg_ratio"].max() * 1.25)
plt.gca().yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1.0))

# Labels and title
plt.ylabel("Proportion of CD4⁺ Tisg Cells", fontsize=13)
plt.xlabel("Tissue Type", fontsize=13)
plt.title("Proportion of CD4⁺ Tisg Cells in Tissues", fontsize=13, fontweight="semibold", pad=15)

# Save figure
plt.tight_layout()
plt.savefig("CD4_Tisg_barplot.png", dpi=900)
plt.show()
