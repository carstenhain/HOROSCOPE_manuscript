"""End-to-end centromere k-mer clustering workflow."""

import argparse
import warnings

import os
import pandas as pd # type: ignore
from pathlib import Path

from clustering import (
    add_cluster_membership_columns,
    build_linkage,
    load_or_calculate_similarity,
    select_tagging_nodes,
)
from config import make_config # type: ignore
from data import (
    build_cluster_metadata,
    filter_clustering_kmers,
    load_active_centromeres,
    load_kmer_counts,
    load_permissive_kmers,
    load_population_colors,
    restrict_to_kmers,
    write_sample_order,
    filter_kmers
)
from length import (
    add_hor_measurements,
    find_normalization_kmers,
    correlate_length_kmers,
)
from output import write_cluster_assignments, write_tagging_kmers
from plots import (
    plot_clustering,
    plot_kmer_count_distribution,
    plot_tagging_kmers,
    plot_tagging_nodes,
)
from tagging import empty_kmer_results, find_tagging_kmers


def _find_total_length_kmers(config, all_kmers, permissive_kmers, permissive_haplotypes, active_centromeres):
    """Find k-mers correlated with HOR length in training haplotypes only."""

    # Keep sufficiently abundant k-mers before testing their association with length.
    length_candidates = all_kmers.loc[:, all_kmers.max() > 20].copy()

    # Evaluate the subset restricted to accepted k-mers and non-test haplotypes.
    candidate_sets = {
        "UNIQUE": restrict_to_kmers(length_candidates, permissive_kmers, permissive_haplotypes),
    }

    hor_lengths = pd.Series(
        active_centromeres["HOR_LENGTH"].to_numpy(),
        index=[
            "{}_{}".format(config.chromosome, sample_id)
            for sample_id in active_centromeres["ID"]
        ],
    )
    result_tables = []

    # Correlate each k-mer with the target and label results for downstream output.
    for mode, kmer_counts in candidate_sets.items():
        target = hor_lengths.loc[kmer_counts.index].to_numpy()
        print("Filtered data for length prediction: {}".format(kmer_counts.shape))
        result = correlate_length_kmers(kmer_counts, target, config.chromosome)
        cluster_name = "{}_LENGTH_HOR".format(config.chromosome)
        if mode == "UNIQUE":
            cluster_name += "_UNIQUE"
        result["CLUSTER"] = cluster_name
        result_tables.append(result)
        print("Found {} kmers".format(result.shape[0]))
    return pd.concat(result_tables, ignore_index=True)


def run_workflow(config):
    """Run clustering, tag discovery, global length prediction, and normalization."""

    # Suppress the expected single-threaded warning emitted by nested joblib calls.
    warnings.filterwarnings("ignore", message="n_jobs value .* overridden to 1.*")

    # Load active samples, their population annotations, and the full k-mer matrix.
    active_centromeres = load_active_centromeres(config)
    population_colors = load_population_colors(config)
    all_kmers = load_kmer_counts(config, active_centromeres)

    # Retain the variable k-mers used to calculate centromere similarity.
    clustering_kmers = filter_clustering_kmers(all_kmers) # type: ignore
    print("Load full data with shape: {}".format(all_kmers.shape))
    print("Filtered data to receive array with shape: {}".format(clustering_kmers.shape))
    plot_kmer_count_distribution(config, all_kmers, clustering_kmers)

    # Cluster samples from their pairwise k-mer similarity and record display order.
    metadata = build_cluster_metadata(config, active_centromeres, population_colors)
    similarity = load_or_calculate_similarity(config, clustering_kmers)
    linkage = build_linkage(similarity)
    dendrogram_order = plot_clustering(
        config, active_centromeres, metadata, similarity, linkage
    )
    write_sample_order(config, active_centromeres, dendrogram_order)

    # Select informative internal nodes and prepare k-mers for cluster-tag discovery.
    tagging_nodes = select_tagging_nodes(linkage, config)
    filter_kmers(all_kmers, config)
    permissive_kmers = load_permissive_kmers(config)
    present_in = all_kmers.shape[0] - (all_kmers == 0).sum()
    # prefiltering of tagging candidates to remove super-rare or super common k-mers --> performance reasons as there are no cluster > 90% and very small clusters
    tagging_candidates = all_kmers.loc[
        :, (present_in > 2) & (present_in < 0.9 * all_kmers.shape[0])
    ]

    # Exclude held-out test samples when learning tagging, length, and normalization k-mers.
    sample_blacklist = []
    with open("/g/korbel/hain/centromere/revision/HOROSCOPE_manuscript/test_samples_1kGP_seed_10.txt", "r") as f:
        for line in f:
            sample_blacklist.append(line.strip())

    # Map each haplotype ID to its sample ID and retain only training haplotypes.
    permissive_haplotypes = []
    for hp in all_kmers.index:
        sample_name = hp.split("_")[1]
        if not (sample_name in sample_blacklist):
            permissive_haplotypes.append(hp)
    print(f"{len(permissive_haplotypes)} are used for kmer generation (from originally {all_kmers.shape[0]} haplotypes)")

    # Retain both selected k-mers and integer-labelled cluster membership columns.
    cluster_ids = [node.id for node in tagging_nodes]
    permissive_kmers_and_clusters = set(permissive_kmers).union(cluster_ids)
    tagging_data = restrict_to_kmers(
        add_cluster_membership_columns(tagging_candidates, tagging_nodes),
        permissive_kmers_and_clusters,
        permissive_haplotypes,
    )
    plot_tagging_nodes(config, active_centromeres, linkage, tagging_nodes)
    tagging_kmers = find_tagging_kmers(
        tagging_data, tagging_nodes, config.chromosome
    )
    plot_tagging_kmers(
        config,
        tagging_data,
        tagging_kmers,
        linkage,
        tagging_nodes,
        plotting_kmer_counts=tagging_candidates,
    )

    # Find HOR-length predictors using the same training-only k-mer subset.
    centromeres_with_hors = add_hor_measurements(config, active_centromeres)
    total_length_kmers = _find_total_length_kmers(
        config, all_kmers, permissive_kmers, permissive_haplotypes, centromeres_with_hors
    )

    # Identify k-mers suitable for normalizing count-based measurements.
    normalization_data = restrict_to_kmers(all_kmers, permissive_kmers, permissive_haplotypes)
    normalization_kmers = find_normalization_kmers(normalization_data, config)

    # Combine all result types, then persist both k-mers and sample-cluster assignments.
    output_kmers = pd.concat(
        [
            empty_kmer_results(),
            tagging_kmers,
            total_length_kmers,
            normalization_kmers,
        ],
        ignore_index=True,
    )
    write_tagging_kmers(config, output_kmers)
    write_cluster_assignments(config, active_centromeres, tagging_nodes)
    return output_kmers


def parse_arguments():
    """Parse command-line settings for one chromosome workflow run."""

    # Define the command-line interface for selecting data and similarity computation.
    parser = argparse.ArgumentParser(
        description=(
            "Cluster centromeres and find cluster tags, whole-centromere length "
            "predictors, and normalization k-mers."
        )
    )
    parser.add_argument("--chrom", required=True, help="Chromosome, for example chr3")
    parser.add_argument(
        "--recalculate-similarity",
        action="store_true",
        help="Recalculate the pairwise k-mer similarity matrix instead of using its cache.",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=16,
        help="Worker processes for pairwise similarity calculation (default: 16).",
    )

    # Return validated command-line settings to the workflow entry point.
    return parser.parse_args()


def main():
    """Run the workflow from the command line."""

    # Convert CLI options into the shared configuration consumed by all workflow stages.
    arguments = parse_arguments()
    recalculate_similarity = True if arguments.recalculate_similarity else None
    config = make_config(
        arguments.chrom,
        recalculate_similarity=recalculate_similarity,
        workers=arguments.workers,
    )
    tagging_kmer_out_path = "/g/korbel/hain/centromere/tagging_kmers_seed_10"
    config.paths.tagging_kmers = Path(tagging_kmer_out_path)
    os.makedirs(tagging_kmer_out_path, exist_ok=True)
    run_workflow(config)


if __name__ == "__main__":
    # Execute only when invoked as a script, not when imported as a module.
    main()