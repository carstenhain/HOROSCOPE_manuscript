"""Pairwise similarity calculations and cluster-node selection."""

import multiprocessing as mp

import numpy as np # type: ignore
import scipy.cluster.hierarchy as hierarchy # type: ignore
from tqdm import tqdm # type: ignore


_PRESENCE_MATRIX = None
_SAMPLE_TOTALS = None


def _similarity_for_pair(pair):
    first_index, second_index = pair
    shared_kmers = np.minimum(
        _PRESENCE_MATRIX[first_index], _PRESENCE_MATRIX[second_index] # type: ignore
    ).sum()
    denominator = _SAMPLE_TOTALS[first_index] + _SAMPLE_TOTALS[second_index] # type: ignore
    return first_index, second_index, 2 * shared_kmers / denominator


def calculate_pairwise_similarity(kmer_presence, workers):
    """Calculate the old pairwise Dice similarity score for every sample pair."""

    global _PRESENCE_MATRIX, _SAMPLE_TOTALS

    _PRESENCE_MATRIX = kmer_presence.to_numpy(copy=False)
    _SAMPLE_TOTALS = _PRESENCE_MATRIX.sum(axis=1)
    sample_count = _PRESENCE_MATRIX.shape[0]
    sample_pairs = [
        (first_index, second_index)
        for first_index in range(sample_count)
        for second_index in range(first_index, sample_count)
    ]
    similarity = np.zeros((sample_count, sample_count))

    if workers > 1 and mp.get_start_method() == "fork":
        with mp.Pool(workers) as pool:
            results = pool.imap_unordered(_similarity_for_pair, sample_pairs)
            for first_index, second_index, value in tqdm(results, total=len(sample_pairs)):
                similarity[first_index, second_index] = value
                similarity[second_index, first_index] = value
    else:
        for first_index, second_index in tqdm(sample_pairs):
            _, _, value = _similarity_for_pair((first_index, second_index))
            similarity[first_index, second_index] = value
            similarity[second_index, first_index] = value

    _PRESENCE_MATRIX = None
    _SAMPLE_TOTALS = None
    return similarity


def load_or_calculate_similarity(config, kmer_counts):
    """Load a cached similarity matrix or calculate and cache it."""

    path = config.paths.count_arrays / "{}_pairwise_similiarity.filtered_data.pairwise.npy".format(
        config.chromosome
    )
    if not config.recalculate_similarity:
        return np.load(path)

    similarity = calculate_pairwise_similarity(kmer_counts.clip(upper=1), config.workers)
    np.save(path, similarity)
    return similarity


def build_linkage(similarity):
    """Build the weighted linkage matrix used by the original workflow."""
    
    """
    # cleaner solution would be to use squareform on the distances --> reorders the haplotypes --> functionally similiar results though
    distance_matrix = 1 - similarity
    np.fill_diagonal(distance_matrix, 0)
    condensed_distances = squareform(distance_matrix, checks=False)
    return hierarchy.linkage(condensed_distances, method="weighted")
    """
    
    return hierarchy.linkage(similarity, method="weighted")


def select_tagging_nodes(linkage, config):
    """Select curated dendrogram nodes within the configured tree-height range."""

    root = hierarchy.to_tree(linkage)
    nodes_in_range = []
    stack = [root]

    while stack:
        node = stack.pop()
        inspect_children = False
        if config.tree_min <= node.dist <= config.tree_max:
            nodes_in_range.append(node)
            inspect_children = True
        elif node.dist > config.tree_max:
            inspect_children = True

        if inspect_children and not node.is_leaf():
            if node.left.dist < config.tree_min and not node.left.is_leaf():
                nodes_in_range.append(node.left)
            if node.right.dist < config.tree_min and not node.right.is_leaf():
                nodes_in_range.append(node.right)

        if node.right:
            stack.append(node.right)
        if node.left:
            stack.append(node.left)

    return [node for node in nodes_in_range if node.id in config.cluster_ids]


def add_cluster_membership_columns(kmer_counts, tagging_nodes):
    """Add one binary membership column per selected dendrogram node."""

    counts_with_clusters = kmer_counts.copy()
    for node in tagging_nodes:
        membership = np.zeros(counts_with_clusters.shape[0])
        membership[list(node.pre_order(lambda child: child.id))] = 1
        counts_with_clusters[node.id] = membership
    return counts_with_clusters