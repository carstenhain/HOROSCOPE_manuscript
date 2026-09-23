"""Figures produced by the centromere clustering workflow."""

import numpy as np # type: ignore
import matplotlib.pyplot as plt # type: ignore
from matplotlib import gridspec # type: ignore
from matplotlib import patches # type: ignore
import scipy.cluster.hierarchy as hierarchy # type: ignore


def plot_kmer_count_distribution(config, all_kmers, filtered_kmers):
    """Plot k-mer prevalence and maximum-copy-number distributions."""

    figure, axes = plt.subplots(2, 3, figsize=(12, 4))
    all_presence = all_kmers.shape[0] - (all_kmers == 0).sum()
    filtered_presence = filtered_kmers.shape[0] - (filtered_kmers == 0).sum()

    axes[0, 0].hist(all_presence, bins=all_kmers.shape[0])
    axes[0, 0].set_yscale("log")
    axes[0, 0].set(xlabel="Num samples", ylabel="Count", title="61mers present in x samples")
    axes[0, 1].hist(all_kmers.max(), stacked=True, bins=100)
    axes[0, 1].set_yscale("log")
    axes[0, 1].set(xlabel="Max. count", ylabel="Count", title="Maximum count of 61mer")
    axes[0, 2].scatter(all_presence, all_kmers.max(), s=0.2)
    axes[0, 2].set(xlabel="61mers present in x samples", ylabel="Maximum count of 61mer")

    axes[1, 0].hist(filtered_presence, bins=max(1, int(filtered_kmers.shape[0] * 0.8)))
    axes[1, 0].set_yscale("log")
    axes[1, 0].set(xlabel="Num samples", ylabel="Count", title="Filtered 61mers present in x samples")
    axes[1, 1].hist(filtered_kmers.max(), stacked=True, bins=100)
    axes[1, 1].set_yscale("log")
    axes[1, 1].set(xlabel="Max. count", ylabel="Count", title="Filtered maximum count of 61mer")
    axes[1, 2].scatter(filtered_presence, filtered_kmers.max(), s=0.2)
    axes[1, 2].set(xlabel="61mers present in x samples", ylabel="Maximum count of 61mer")

    figure.tight_layout()
    figure.savefig(
        config.paths.figures / "{}_kmer_count_distribution.png".format(config.chromosome)
    )
    plt.close(figure)


def plot_clustering(config, active_centromeres, metadata, similarity, linkage):
    """Save the dendrogram, similarity matrix, and metadata tracks."""

    figure = plt.figure(figsize=(80, 80))
    grid = gridspec.GridSpec(
        2, 5, width_ratios=[32, 4, 1, 1, 1], height_ratios=[1, 4]
    )
    dendrogram_axis = plt.subplot(grid[0, 0])
    dendrogram = hierarchy.dendrogram(linkage, ax=dendrogram_axis, orientation="top")
    dendrogram_axis.axis("off")
    dendrogram_axis.set_xticks([])
    order = dendrogram["leaves"]

    matrix_axis = plt.subplot(grid[1, 0])
    ordered_similarity = similarity[order, :][:, order]
    image = matrix_axis.imshow(ordered_similarity, aspect="auto", origin="lower", cmap="viridis")

    contig_axis = plt.subplot(grid[1, 1])
    contig_lengths = active_centromeres["CONTIG_LENGTH"].astype(int).to_numpy()[order]
    contig_axis.barh(range(active_centromeres.shape[0]), contig_lengths, height=1, label="Contig")
    contig_axis.set_ylim(-0.5, active_centromeres.shape[0] - 0.5)
    contig_axis.set(xlabel="LENGTH")
    contig_axis.set_yticks([])

    population_axis = plt.subplot(grid[1, 2])
    population_axis.barh(
        range(metadata.shape[0]),
        width=1,
        height=1,
        color=np.asarray(metadata["SUPERPOPCOL"].tolist())[order],
    )
    population_axis.set(xlabel="SUPERPOP", xlim=(0, 1), ylim=(-0.5, metadata.shape[0] - 0.5))
    population_axis.set_xticks([])
    population_axis.set_yticks([])

    color_axis = plt.subplot(grid[1, 4])
    colorbar = plt.colorbar(image, cax=color_axis)
    colorbar.ax.set_ylabel("Fraction shared kmers", rotation=-90, va="bottom")
    ordered_sample_ids = active_centromeres["ID"].to_numpy()[order]
    matrix_axis.set_xticks(np.arange(active_centromeres.shape[0]))
    matrix_axis.set_xticklabels(ordered_sample_ids, rotation=-90)
    matrix_axis.set_yticks(np.arange(active_centromeres.shape[0]))
    matrix_axis.set_yticklabels(ordered_sample_ids)

    figure.tight_layout()
    figure.savefig(config.paths.figures / "{}_cluster.png".format(config.chromosome))
    matrix_axis.set_rasterized(True)
    figure.savefig(config.paths.figures / "{}_cluster.svg".format(config.chromosome))
    plt.close(figure)
    return order


def _dendrogram_coordinates(dendrogram):
    return {
        distances[1]: np.mean([coordinates[1], coordinates[2]])
        for distances, coordinates in zip(dendrogram["dcoord"], dendrogram["icoord"])
    }


def plot_tagging_nodes(config, active_centromeres, linkage, tagging_nodes):
    """Show selected tagging nodes over the dendrogram."""

    figure, axis = plt.subplots(figsize=(25, 4))
    dendrogram = hierarchy.dendrogram(
        linkage, ax=axis, labels=active_centromeres["ID"].tolist(), orientation="top"
    )
    coordinates = _dendrogram_coordinates(dendrogram)
    axis.add_patch(
        patches.Rectangle(
            (np.amin(dendrogram["icoord"]), config.tree_min),
            np.amax(dendrogram["icoord"]) - np.amin(dendrogram["icoord"]),
            config.tree_max - config.tree_min,
            facecolor="lightblue",
            alpha=0.2,
        )
    )
    for node in tagging_nodes:
        axis.scatter(coordinates[node.dist], node.dist, c="r", s=20)
        axis.text(coordinates[node.dist], node.dist, str(node.id))

    figure.savefig(config.paths.figures / "{}_tree_with_nodes.png".format(config.chromosome))
    figure.savefig(config.paths.figures / "{}_tree_with_nodes.svg".format(config.chromosome))
    plt.close(figure)


def plot_tagging_kmers(config, kmer_counts, tagging_kmers, linkage, tagging_nodes):
    """Plot the best tagging k-mer for each selected dendrogram node."""

    nodes_with_kmers = [
        node
        for node in tagging_nodes
        if not tagging_kmers[tagging_kmers["CLUSTER"] == "{}_{}".format(config.chromosome, node.id)].empty
    ]
    figure = plt.figure(figsize=(25, 5 + len(nodes_with_kmers) / 2))
    grid = gridspec.GridSpec(2, 1, height_ratios=[1, max(0.1, len(nodes_with_kmers) / 10)])
    tree_axis = plt.subplot(grid[0, 0])
    heatmap_axis = plt.subplot(grid[1, 0])
    dendrogram = hierarchy.dendrogram(linkage, ax=tree_axis, orientation="top")
    coordinates = _dendrogram_coordinates(dendrogram)

    heatmap_data = []
    heatmap_labels = []
    colour_map = plt.cm.viridis # type: ignore
    normalization = plt.Normalize(vmin=0, vmax=30) # type: ignore
    for node in nodes_with_kmers:
        node_kmers = tagging_kmers[
            tagging_kmers["CLUSTER"] == "{}_{}".format(config.chromosome, node.id)
        ]
        for label in ("ALL", "UNIQUE", "RELATIVE"):
            label_kmers = node_kmers if label == "ALL" else node_kmers[node_kmers["TYPE"] == label]
            if label_kmers.empty:
                continue
            best_kmer = label_kmers.sort_values(
                by="MEAN_COUNT_IN", ascending=False
            )["KMER"].iloc[0]
            heatmap_data.append(kmer_counts[best_kmer].to_numpy()[dendrogram["leaves"]])
            heatmap_labels.append("{}_{}".format(node.id, label))

        size = 20
        if node_kmers.shape[0] > 100:
            size = 40
        if node_kmers.shape[0] > 1000:
            size = 100
        tree_axis.scatter(coordinates[node.dist], node.dist, color="k", s=size + 5, zorder=9)
        tree_axis.scatter(
            coordinates[node.dist],
            node.dist,
            color=colour_map(normalization(node_kmers["MEAN_COUNT_IN"].max())),
            s=size,
            zorder=10,
        )
        tree_axis.text(coordinates[node.dist], node.dist, str(node.id), zorder=11)

    if heatmap_data:
        heatmap_axis.imshow(
            np.clip(np.asarray(heatmap_data), a_min=0, a_max=30),
            aspect="auto",
            cmap="viridis",
        )
        heatmap_axis.set_yticks(range(len(heatmap_labels)), heatmap_labels)
    else:
        heatmap_axis.set_axis_off()

    figure.savefig(
        config.paths.figures / "{}_tagging_nodes_kmers.png".format(config.chromosome)
    )
    figure.savefig(
        config.paths.figures / "{}_tagging_nodes_kmers.svg".format(config.chromosome)
    )
    plt.close(figure)