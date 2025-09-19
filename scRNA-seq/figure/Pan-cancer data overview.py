import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle

def plot_adata_overview(adata, cancer_type_col, sample_col, patient_col, source_col, figsize=(12, 8)):
    """
    Generate a multi-panel summary plot of an AnnData object, grouped by cancer type.

    Panels:
    1. Cell count per cancer type (log scale)
    2. Sample count per cancer type (log scale)
    3. Subject count per cancer type (log scale)
    4. Tissue source composition (stacked bar plot)

    Parameters
    ----------
    adata : AnnData
        Annotated AnnData object with .obs columns for grouping and metadata.
    cancer_type_col : str
        Column name in .obs for cancer type.
    sample_col : str
        Column name in .obs for sample ID.
    patient_col : str
        Column name in .obs for patient ID.
    source_col : str
        Column name in .obs for tissue source (e.g., Tumor, Normal).
    figsize : tuple
        Size of the entire figure.
    
    Returns
    -------
    fig, axes, stats_df, comp_df
    """

    obs_df = adata.obs.copy()

    # Summary stats per cancer type
    stats_list = []
    for cancer_type in obs_df[cancer_type_col].unique():
        subset = obs_df[obs_df[cancer_type_col] == cancer_type]
        stats_list.append({
            'cancer_type': cancer_type,
            'cell_count': len(subset),
            'sample_count': subset[sample_col].nunique(),
            'patient_count': subset[patient_col].nunique()
        })
    stats_df = pd.DataFrame(stats_list).sort_values("cell_count", ascending=True)

    # Tissue source composition
    composition_data = []
    for cancer_type in stats_df["cancer_type"]:
        subset = obs_df[obs_df[cancer_type_col] == cancer_type]
        source_counts = subset[source_col].value_counts()
        total = len(subset)
        comp_dict = {"cancer_type": cancer_type}
        for source in source_counts.index:
            comp_dict[source] = source_counts[source] / total
        composition_data.append(comp_dict)
    comp_df = pd.DataFrame(composition_data).fillna(0)
    comp_df = comp_df.set_index("cancer_type").reindex(stats_df["cancer_type"])

    # Setup figure and color palette
    fig, axes = plt.subplots(1, 4, figsize=figsize, gridspec_kw={'width_ratios': [1, 1, 1, 1.5]})
    bar_colors = ['#7FB3D380', '#C5A3C780', '#F4A46080']  # semi-transparent

    # Panel 1: Cell count
    ax1 = axes[0]
    ax1.barh(range(len(stats_df)), stats_df["cell_count"], color=bar_colors[0])
    ax1.set_yticks(range(len(stats_df)))
    ax1.set_yticklabels(stats_df["cancer_type"], fontsize=10)
    ax1.set_xlabel("Count")
    ax1.set_title("Cell")
    ax1.set_xscale("log")

    # Panel 2: Sample count
    ax2 = axes[1]
    ax2.barh(range(len(stats_df)), stats_df["sample_count"], color=bar_colors[1])
    ax2.set_yticks(range(len(stats_df)))
    ax2.set_yticklabels([])  # hide
    ax2.set_xlabel("Count")
    ax2.set_title("Sample")
    ax2.set_xscale("log")

    # Panel 3: Subject (patient) count
    ax3 = axes[2]
    ax3.barh(range(len(stats_df)), stats_df["patient_count"], color=bar_colors[2])
    ax3.set_yticks(range(len(stats_df)))
    ax3.set_yticklabels([])
    ax3.set_xlabel("Count")
    ax3.set_title("Subject")
    ax3.set_xscale("log")

    # Panel 4: Tissue source composition (stacked bar)
    ax4 = axes[3]
    sources = list(comp_df.columns)
    predefined_colors = {
        'Blood': '#7FB3D3', 'Normal': '#C5A3C7', 'Tumor': '#F4A460', 'Metastatic': '#E74C3C',
        'Primary': '#F4A460', 'Met': '#E74C3C', 'Norm': '#C5A3C7'
    }
    source_colors = {}
    default_colors = ['#95A5A6', '#34495E', '#2ECC71', '#9B59B6', '#F39C12']
    color_idx = 0
    for source in sources:
        matched = False
        for key in predefined_colors:
            if key.lower() in source.lower():
                source_colors[source] = predefined_colors[key]
                matched = True
                break
        if not matched:
            source_colors[source] = default_colors[color_idx % len(default_colors)]
            color_idx += 1

    # Define source drawing order
    draw_order = ['Blood', 'Normal', 'Tumor', 'Metastatic']
    ordered_sources = []
    for preferred in draw_order:
        for source in sources:
            if preferred.lower() in source.lower() and source not in ordered_sources:
                ordered_sources.append(source)
    for source in sources:
        if source not in ordered_sources:
            ordered_sources.append(source)

    # Draw stacked bars
    left = np.zeros(len(comp_df))
    for source in ordered_sources:
        if source in comp_df.columns:
            values = comp_df[source].values
            ax4.barh(range(len(comp_df)), values, left=left, color=source_colors[source], label=source)
            left += values
    ax4.set_yticks(range(len(stats_df)))
    ax4.set_yticklabels([])
    ax4.set_xlabel("Cell fraction")
    ax4.set_xlim(0, 1)
    ax4.set_title("Tissue")

    # Legend
    handles = [Rectangle((0, 0), 1, 1, color=source_colors[s]) for s in ordered_sources if s in comp_df.columns]
    labels = [s for s in ordered_sources if s in comp_df.columns]
    ax4.legend(handles, labels, bbox_to_anchor=(1.05, 1), loc='upper left', fontsize='small')

    plt.tight_layout()
    return fig, axes, stats_df, comp_df
fig, axes, stats_df, comp_df = plot_adata_overview(
    adata,
    cancer_type_col='Cancer_Type',
    sample_col='Sample',
    patient_col='Patient',
    source_col='Source_reclassified'
)
plt.savefig("adata_overview.pdf", bbox_inches='tight')
plt.show()
