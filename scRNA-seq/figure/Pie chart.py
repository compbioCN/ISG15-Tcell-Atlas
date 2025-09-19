import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def plot_cell_type_pie_chart(adata, 
                              font_size=16,
                              figsize=(18, 14),
                              label_distance=1.65,
                              output_pdf='cell_type_distribution.pdf'):
    """
    Draws a pie chart showing the distribution of cell types.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with 'cell_type' in `.obs` and 'cell_type_colors' in `.uns`.
    font_size : int, optional
        Font size for labels. Default is 16.
    figsize : tuple, optional
        Size of the figure. Default is (18, 14).
    label_distance : float, optional
        Distance of labels from the pie center. Default is 1.65.
    output_pdf : str, optional
        Output filename for the PDF. PNG will also be saved.
    """

    # Count cells per cell type
    cell_type_counts = adata.obs['cell_type'].value_counts()
    total_cells = cell_type_counts.sum()
    cell_type_percentages = (cell_type_counts / total_cells * 100).round(2)

    # Color mapping from adata.uns
    cell_types_list = adata.obs['cell_type'].cat.categories.tolist()
    cell_type_colors = adata.uns['cell_type_colors']
    color_dict = dict(zip(cell_types_list, cell_type_colors))
    colors = [color_dict.get(ct, '#808080') for ct in cell_type_counts.index]

    fig, ax = plt.subplots(figsize=figsize)

    # Create pie wedges
    wedges, _ = ax.pie(cell_type_counts.values,
                       colors=colors,
                       startangle=90,
                       counterclock=False,
                       wedgeprops=dict(width=1, edgecolor='white', linewidth=2))

    # Adjust label angles to avoid overlap
    def adjust_label_positions(angles, min_distance=25):
        n = len(angles)
        adjusted = angles.copy()
        sorted_indices = sorted(range(n), key=lambda i: angles[i])
        for i in range(1, n):
            idx = sorted_indices[i]
            prev_idx = sorted_indices[i-1]
            if adjusted[idx] - adjusted[prev_idx] < min_distance:
                adjusted[idx] = adjusted[prev_idx] + min_distance
        if 360 - adjusted[sorted_indices[-1]] + adjusted[sorted_indices[0]] < min_distance:
            total_adjustment = (n * min_distance - 360) / n
            for i, idx in enumerate(sorted_indices):
                adjusted[idx] -= total_adjustment * i
        return adjusted

    wedge_angles = [(w.theta2 - w.theta1) / 2 + w.theta1 for w in wedges]
    adjusted_angles = adjust_label_positions(wedge_angles, min_distance=30)

    for i, (wedge, ct, count, pct) in enumerate(zip(wedges, cell_type_counts.index,
                                                    cell_type_counts.values, cell_type_percentages.values)):
        # Original wedge angle
        wedge_center = (wedge.theta2 - wedge.theta1) / 2 + wedge.theta1
        wedge_center_rad = np.deg2rad(wedge_center)
        label_angle_rad = np.deg2rad(adjusted_angles[i])

        # Line start & elbow
        x_start = 1.0 * np.cos(wedge_center_rad)
        y_start = 1.0 * np.sin(wedge_center_rad)
        x_elbow = 1.15 * np.cos(wedge_center_rad)
        y_elbow = 1.15 * np.sin(wedge_center_rad)

        # Label position
        x_label = label_distance * np.cos(label_angle_rad)
        y_label = label_distance * np.sin(label_angle_rad)
        ha = 'left' if x_label > 0 else 'right'
        x_end = x_label - 0.05 if x_label > 0 else x_label + 0.05

        # Draw connector lines
        ax.plot([x_start, x_elbow], [y_start, y_elbow], 'gray', lw=1.5, alpha=0.6)
        ax.plot([x_elbow, x_end], [y_elbow, y_label], 'gray', lw=1.5, alpha=0.6)

        # Label text
        label = f"{ct}\n({count:,}; {pct}%)"
        bbox_props = dict(boxstyle="round,pad=0.6", facecolor="white",
                          edgecolor="lightgray", linewidth=0.8, alpha=0.9)

        ax.text(x_label, y_label, label,
                ha=ha, va='center',
                fontsize=font_size,
                bbox=bbox_props,
                fontfamily='sans-serif')

    ax.set_aspect('equal')
    ax.set_xlim(-label_distance - 0.8, label_distance + 0.8)
    ax.set_ylim(-label_distance - 0.8, label_distance + 0.8)
    ax.axis('off')
    plt.tight_layout()
    plt.subplots_adjust(left=0.02, right=0.98, top=0.98, bottom=0.02)

    # Save as PDF and PNG
    plt.savefig(output_pdf, format='pdf', dpi=300, bbox_inches='tight')
    png_out = output_pdf.replace('.pdf', '.png')
    plt.savefig(png_out, format='png', dpi=300, bbox_inches='tight')
    print(f"✓ Saved: {output_pdf}, {png_out}")
    plt.show()

    # Console summary
    print("\n" + "="*60)
    print("Cell Type Summary")
    print("="*60)
    for ct, count, pct in zip(cell_type_counts.index, cell_type_counts.values, cell_type_percentages.values):
        print(f"{ct:<25} {count:>10,} cells ({pct:>6.2f}%)")
    print("-"*60)
    print(f"{'Total':<25} {total_cells:>10,} cells (100.00%)")
    print("="*60)
# Call the function (basic)
plot_cell_type_pie_chart(adata)

# Optional: customize font and layout
plot_cell_type_pie_chart(
    adata,
    font_size=20,
    figsize=(22, 18),
    label_distance=1.8,
    output_pdf="my_celltype_pie.pdf"
)
