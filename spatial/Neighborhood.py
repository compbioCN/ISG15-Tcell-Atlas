import scanpy as sc
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import networkx as nx
from mpl_chord_diagram import chord_diagram
import squidpy as sq

# Ensure Tcelltype exists and is categorical
if 'Tcelltype' not in cdata.obs:
    raise ValueError("Tcelltype column is missing in cdata.obs")
cdata.obs['Tcelltype'] = cdata.obs['Tcelltype'].astype('category')

# Spatial neighbor graph
sq.gr.spatial_neighbors(cdata, coord_type="generic", delaunay=False, n_neighs=10)
sq.pl.spatial_scatter(
    cdata,
    color="Tcelltype",
    connectivity_key="spatial_connectivities",
    size=1,
    figsize=(10, 10),
)

# Neighborhood enrichment
level_ = 'Tcelltype'
sq.gr.nhood_enrichment(cdata, cluster_key=level_)

# Visualize enrichment heatmap with custom color and settings
custom_cmap = mcolors.LinearSegmentedColormap.from_list("custom_cmap", ["#a9dce6", "#E05F48"])
plt.figure(figsize=(12, 12))
sq.pl.nhood_enrichment(
    cdata,
    cluster_key=level_,
    figsize=(6, 6),
    cmap=custom_cmap,
    fontsize=10,
    annot=True,
    vmin=-50,
    vmax=50,
)
plt.savefig("Epi-T_neighbor_enrichment.pdf", bbox_inches='tight', dpi=900)
plt.close()

# Network graph based on z-score > 1
zscore_matrix = np.nan_to_num(cdata.uns[f"{level_}_nhood_enrichment"]['zscore'])
categories = cdata.obs[level_].cat.categories
G = nx.Graph()
for i, source in enumerate(categories):
    for j, target in enumerate(categories):
        if i >= j:
            continue
        zscore = zscore_matrix[i, j]
        if zscore > 1:
            G.add_edge(source, target, weight=zscore)

pos = nx.spring_layout(G, seed=42)
node_size = cdata.obs[level_].value_counts(sort=False).loc[categories].values
node_colors = sns.color_palette("Spectral", len(categories))
plt.figure(figsize=(10, 10))
nx.draw_networkx_nodes(G, pos, node_size=node_size, node_color=node_colors, alpha=0.8)
nx.draw_networkx_edges(G, pos, alpha=0.5)
nx.draw_networkx_labels(G, pos, font_size=12, font_color="black")
plt.title("Tcelltype Interaction Network")
plt.axis("off")
plt.show()

# Chord diagram of top interactions
sq.gr.interaction_matrix(cdata, cluster_key=level_, normalized=False)
interaction_matrix = pd.DataFrame(
    cdata.uns[f"{level_}_interactions"],
    index=cdata.obs[level_].cat.categories,
    columns=cdata.obs[level_].cat.categories,
)
thresh = interaction_matrix.sum().quantile(0.6)
df_filt = interaction_matrix.loc[:, interaction_matrix.sum() > thresh]
df_filt = df_filt.loc[df_filt.index.intersection(df_filt.columns)]

with plt.rc_context({"figure.figsize": (10, 10), "figure.dpi": 100}):
    chord_diagram(
        df_filt,
        names=list(df_filt.columns),
        rotate_names=True,
        fontsize=10,
        alpha=0.9,
    )
    plt.title("Tcelltype Interaction Matrix")
    plt.savefig("Tcelltype_interaction_chord.svg", bbox_inches="tight")
    plt.show()
