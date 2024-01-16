import pandas as pd
import numpy as np
import hdbscan
from sklearn.metrics import silhouette_score
import math
from sklearn.metrics import adjusted_rand_score

def truncate(number, decimals=2):
    if number == "N/A":
        return number
    else:
        stepper = 10.0 ** decimals
        return math.trunc(stepper * number) / stepper

def joint_matrix_factorization(matrices, num_factors, learning_rate=0.01, iterations=1000, tolerance=1e-4):
    """
    A simple implementation of Joint Matrix Factorization with convergence check.

    :param matrices: List of matrices to factorize.
    :param num_factors: Number of factors to decompose into.
    :param learning_rate: Learning rate for gradient descent.
    :param iterations: Number of iterations for the optimization loop.
    :param tolerance: Threshold for convergence.
    :return: A tuple of the shared matrix and a list of unique matrices.
    """
        
    num_matrices = len(matrices)
    shared_matrix = np.random.rand(matrices[0].shape[0], num_factors) * 0.01
    unique_matrices = [np.random.rand(num_factors, matrix.shape[1]) * 0.01 for matrix in matrices]

    for it in range(iterations):
        total_error = 0

        for i in range(num_matrices):
            error = matrices[i] - np.dot(shared_matrix, unique_matrices[i])
            total_error += np.mean(error ** 2)

            shared_gradient = -2 * np.dot(error, unique_matrices[i].T)
            unique_gradient = -2 * np.dot(shared_matrix.T, error)

            # Gradient clipping
            shared_gradient = np.clip(shared_gradient, -1, 1)
            unique_gradient = np.clip(unique_gradient, -1, 1)

            shared_matrix -= learning_rate * shared_gradient
            unique_matrices[i] -= learning_rate * unique_gradient

        avg_error = total_error / num_matrices

        if avg_error < tolerance:
            print(f"Convergence reached at iteration {it+1}")
            break
        elif np.isnan(avg_error): # Check for NaN values
            print("NaN encountered in error calculation.")
            break

    return shared_matrix, unique_matrices

def calculate_cluster_purity(clusters, lineages):
    unique_clusters = np.unique(clusters)
    cluster_purity_list = []

    for cluster in unique_clusters:
        if cluster == -1: # Exclude noise points
            continue
        # Subset lineages in the current cluster
        lineages_in_cluster = lineages[clusters == cluster]
        # Most common lineage in the cluster
        most_common_lineage = lineages_in_cluster.value_counts().idxmax()
        # Purity calculation
        purity = lineages_in_cluster.value_counts().max() / len(lineages_in_cluster)
        cluster_purity_list.append(purity)

    # Overall purity (weighted by cluster size)
    cluster_sizes = pd.Series(clusters[clusters != -1]).value_counts()
    overall_purity = sum(np.array(cluster_purity_list) * cluster_sizes) / cluster_sizes.sum()

    return overall_purity

def safe_silhouette_score(embedding, labels):
    unique_labels = np.unique(labels)
    if len(unique_labels[unique_labels != -1]) > 1:
        return silhouette_score(embedding, labels)
    else:
        return None

def apply_hdbscan(embedding, lineages, min_size=10):
    clusterer = hdbscan.HDBSCAN(min_cluster_size=min_size, gen_min_span_tree=True)
    cluster_labels = clusterer.fit_predict(embedding)

    # Exclude noise points for ARI calculation
    valid_indices = cluster_labels != -1
    filtered_labels = cluster_labels[valid_indices]
    filtered_lineages = lineages[valid_indices]

    # Calculate purity
    purity = calculate_cluster_purity(cluster_labels, lineages)

    # Calculate Adjusted Rand Index
    ari = adjusted_rand_score(filtered_lineages, filtered_labels)

    # Calculate silhouette
    silhouette = safe_silhouette_score(embedding[valid_indices, :], filtered_labels)

    return purity, ari, silhouette