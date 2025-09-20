import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.neighbors import NearestNeighbors

# Step 1: Extract spatial coordinates and cell types
spatial_coordinates = CD8Tcell.obsm['spatial']
tcelltype = CD8Tcell.obs['Tcelltype']

# Step 2: Compute nearest neighbors (6 neighbors including self)
nbrs = NearestNeighbors(n_neighbors=6, algorithm='auto').fit(spatial_coordinates)
distances, indices = nbrs.kneighbors(spatial_coordinates)

# Step 3: Extract neighbor types (ignore self by skipping index 0)
neighbor_types = pd.DataFrame(
    [[tcelltype[idx] for idx in neighbor_list] for neighbor_list in indices],
    columns=[f'Neighbor_{i}' for i in range(6)]
)
neighbor_types['Tcelltype'] = tcelltype.values

# Step 4: Reshape for distribution analysis
neighbor_types_long = neighbor_types.melt(
    id_vars='Tcelltype',
    value_vars=[f'Neighbor_{i}' for i in range(1, 6)],
    var_name='Neighbor',
    value_name='Neighbor_Tcelltype'
)

# Step 5: Count and normalize neighbor distribution
neighbor_distribution = neighbor_types_long.groupby(['Tcelltype', 'Neighbor_Tcelltype']).size().reset_index(name='Count')
neighbor_distribution['Proportion'] = neighbor_distribution.groupby('Tcelltype')['Count'].transform(lambda x: x / x.sum())

# Step 6: Barplot of neighbor composition
neighbor_pivot = neighbor_distribution.pivot(index='Tcelltype', columns='Neighbor_Tcelltype', values='Proportion').fillna(0)

plt.figure(figsize=(10, 6))
neighbor_pivot.plot(kind='bar', stacked=True, figsize=(10, 6), colormap='tab20')
plt.title('Neighbor Cell Type Distribution for Each Tcelltype')
plt.xlabel('Tcelltype')
plt.ylabel('Proportion of Neighbor Types')
plt.legend(title='Neighbor Tcelltype', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.show()

# Step 7: Assign dominant neighbor type to each cell
CD8Tcell.obs['major_neighbor_type'] = neighbor_types.iloc[:, 1:].mode(axis=1)[0]
