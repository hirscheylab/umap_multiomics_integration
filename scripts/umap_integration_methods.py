import pandas as pd
import numpy as np
from umap import UMAP
import matplotlib.pyplot as plt
import funs_umap_integration_methods as umap_funs
import os

# Create the output directory if it doesn't exist
output_dir = '../results_umap_integration_methods'
os.makedirs(output_dir, exist_ok=True)

# Load Multi-omics Data (Expression and miRNA)
expression_data = pd.read_csv('../data/DepMap/expression.txt', delimiter='\t').T
mirna_data = pd.read_csv('../data/DepMap/mirna.txt', delimiter='\t').T

omics_matrices = [expression_data, mirna_data]

# Load Metadata (lineages)
sample_info = pd.read_csv('../data/DepMap/sample_info.csv')
sample_info['lineage'] = sample_info.apply(
    lambda x: x['lineage_subtype'] if x['lineage'] in ["skin", "blood"] else x['lineage'], axis=1
)
sample_info = sample_info[sample_info['DepMap_ID'].isin(expression_data.index)]
lineages = sample_info['lineage']

# Apply UMAP Individually to Each Omics Dataset
umap_individual = UMAP(n_neighbors=15, n_components=10, min_dist=0.01, metric='euclidean')
umap_individual_embeddings = [umap_individual.fit_transform(omics_data) for omics_data in omics_matrices]

# Apply UMAP to the Individual UMAP Embeddings Concatenation
umap_integrated = UMAP(n_neighbors=10, n_components=2, min_dist=0.1, spread=2, metric='euclidean')

concatenation_mapper = umap_integrated.fit_transform(np.hstack(umap_individual_embeddings))

# Apply 2D UMAP to the Raw Omic Datasets 
expression_mapper = umap_integrated.fit_transform(expression_data)
mirna_mapper = umap_integrated.fit_transform(mirna_data)

# Apply Joint Matrix Factorization (JMF) to UMAP Embeddings
shared_matrix, unique_matrices = umap_funs.joint_matrix_factorization(umap_individual_embeddings, num_factors=10)

# Apply UMAP to the Shared Matrix from JMF
shared_mapper = umap_integrated.fit_transform(shared_matrix)

# Benchmark Methods to Integrate UMAPs
intersection_result = expression_mapper * mirna_mapper
union_result = expression_mapper + mirna_mapper
subtraction_result = expression_mapper - mirna_mapper
concatenate_raw_omics_result = umap_integrated.fit_transform(np.hstack(omics_matrices))

concatenation_mapper # UBMI
shared_mapper

subtraction_result_jmf = concatenation_mapper - shared_mapper
intersection_result_jmf = concatenation_mapper * shared_mapper
union_result_jmf = concatenation_mapper + shared_mapper
concatenation_result_jmf = umap_integrated.fit_transform(np.hstack([np.hstack(umap_individual_embeddings), shared_matrix]))

# Applying HDBSCAN and Compute Metrics for Each Result
intersection_metrics = umap_funs.apply_hdbscan(intersection_result, lineages)
union_metrics = umap_funs.apply_hdbscan(union_result, lineages)
subtraction_metrics = umap_funs.apply_hdbscan(subtraction_result, lineages)
concatenate_raw_omics_metrics = umap_funs.apply_hdbscan(concatenate_raw_omics_result, lineages)

concatenation_metrics = umap_funs.apply_hdbscan(concatenation_mapper, lineages)
shared_metrics = umap_funs.apply_hdbscan(shared_mapper, lineages)

subtraction_metrics_jmf = umap_funs.apply_hdbscan(subtraction_result_jmf, lineages)
intersection_metrics_jmf = umap_funs.apply_hdbscan(intersection_result_jmf, lineages)
union_metrics_jmf = umap_funs.apply_hdbscan(union_result_jmf, lineages)
concatenate_metrics_jmf = umap_funs.apply_hdbscan(concatenation_result_jmf, lineages)

# Plot Results and Metrics
results = [intersection_result, union_result, subtraction_result, concatenate_raw_omics_result,
           concatenation_mapper, shared_mapper,
           subtraction_result_jmf, intersection_result_jmf, union_result_jmf, concatenation_result_jmf]

results_metrics = [
    ("Intersection", intersection_metrics),
    ("Union", union_metrics),
    ("Subtraction", subtraction_metrics),
    ("Concatenate Raw Omics", concatenate_raw_omics_metrics),
    ("Concatenation", concatenation_metrics),
    ("Shared", shared_metrics),
    ("Subtraction JMF", subtraction_metrics_jmf),
    ("Intersection JMF", intersection_metrics_jmf),
    ("Union JMF", union_metrics_jmf),
    ("Concatenate JMF", concatenate_metrics_jmf)
]

fig, axs = plt.subplots(2, 5, figsize=(25, 12))
axs = axs.flatten()

for ax, (title, metrics), result in zip(axs, results_metrics, results):
    purity, ari, silhouette = metrics
    purity, ari, silhouette = map(umap_funs.truncate, [purity, ari, silhouette])
    ax.scatter(result[:, 0], result[:, 1], cmap='viridis', alpha=0.7)
    ax.set_title(f"{title}: Purity: {purity}, Silhouette: {silhouette}") # ARI: {ari}
    ax.set_xlabel('UMAP 1')
    ax.set_ylabel('UMAP 2')

plt.tight_layout()
plt.savefig(f'{output_dir}/results_grid.pdf')
plt.show()