"""End-to-end centromere k-mer clustering workflow."""

import argparse
import warnings

import pandas as pd # type: ignore

from .clustering import (
    add_cluster_membership_columns,
    build_linkage,
    load_or_calculate_similarity,
    select_tagging_nodes,
)
from .config import make_config
from .data import (
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
from .length import (
    add_hor_measurements,
    find_normalization_kmers,
    correlate_length_kmers,
)
from .output import write_cluster_assignments, write_tagging_kmers
from .plots import (
    plot_clustering,
    plot_kmer_count_distribution,
    plot_tagging_kmers,
    plot_tagging_nodes,
)
from .tagging import empty_kmer_results, find_tagging_kmers


def _find_total_length_kmers(config, all_kmers, permissive_kmers, permissive_haplotypes, active_centromeres):
    """Find k-mers for contig, HOR, and HOR-span length prediction."""

    length_candidates = all_kmers.loc[:, all_kmers.max() > 20].copy()
    candidate_sets = {
        "UNIQUE": restrict_to_kmers(length_candidates, permissive_kmers, permissive_haplotypes),
    }
    targets = (
        ("_HOR", active_centromeres["HOR_LENGTH"].values),
    )
    result_tables = []
    for suffix, target in targets:
        for mode, kmer_counts in candidate_sets.items():
            print("Filtered data for length prediction: {}".format(kmer_counts.shape))
            result = correlate_length_kmers(kmer_counts, target, config.chromosome)
            cluster_name = "{}_LENGTH{}".format(config.chromosome, suffix)
            if mode == "UNIQUE":
                cluster_name += "_UNIQUE"
            result["CLUSTER"] = cluster_name
            result_tables.append(result)
            print("Found {} kmers".format(result.shape[0]))
    return pd.concat(result_tables, ignore_index=True)


def run_workflow(config):
    """Run clustering, tag discovery, global length prediction, and normalization."""

    warnings.filterwarnings("ignore", message="n_jobs value .* overridden to 1.*")
    active_centromeres = load_active_centromeres(config)
    population_colors = load_population_colors(config)
    all_kmers = load_kmer_counts(config, active_centromeres)
    clustering_kmers = filter_clustering_kmers(all_kmers)
    print("Load full data with shape: {}".format(all_kmers.shape))
    print("Filtered data to receive array with shape: {}".format(clustering_kmers.shape))
    plot_kmer_count_distribution(config, all_kmers, clustering_kmers)

    metadata = build_cluster_metadata(config, active_centromeres, population_colors)
    similarity = load_or_calculate_similarity(config, clustering_kmers)
    linkage = build_linkage(similarity)
    dendrogram_order = plot_clustering(
        config, active_centromeres, metadata, similarity, linkage
    )
    write_sample_order(config, active_centromeres, dendrogram_order)

    tagging_nodes = select_tagging_nodes(linkage, config)
    filter_kmers(all_kmers, config)
    permissive_kmers = load_permissive_kmers(config)
    present_in = all_kmers.shape[0] - (all_kmers == 0).sum()
    tagging_candidates = all_kmers.loc[
        :, (present_in > 2) & (present_in < 0.9 * all_kmers.shape[0])
    ]
    sample_blacklist = [] # ["HG00096", ...]
    permissive_haplotypes = []
    for hp in all_kmers.index:
        sample_name = hp.split("_")[1]
        if not (sample_name in sample_blacklist):
            permissive_haplotypes.append(hp)
    print(f"{len(permissive_haplotypes)} are used for kmer generation (from originally {all_kmers.shape[0]} haplotypes)")
    tagging_data = add_cluster_membership_columns(
        restrict_to_kmers(tagging_candidates, permissive_kmers, permissive_haplotypes), tagging_nodes
    )
    plot_tagging_nodes(config, active_centromeres, linkage, tagging_nodes)
    tagging_kmers = find_tagging_kmers(
        tagging_data, tagging_nodes, config.chromosome
    )
    plot_tagging_kmers(config, tagging_data, tagging_kmers, linkage, tagging_nodes)

    centromeres_with_hors = add_hor_measurements(config, active_centromeres)
    total_length_kmers = _find_total_length_kmers(
        config, all_kmers, permissive_kmers, permissive_haplotypes, centromeres_with_hors
    )
    normalization_data = restrict_to_kmers(all_kmers, permissive_kmers, permissive_haplotypes)
    normalization_kmers = find_normalization_kmers(normalization_data, config)

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
    return parser.parse_args()


def main():
    """Run the workflow from the command line."""

    arguments = parse_arguments()
    recalculate_similarity = True if arguments.recalculate_similarity else None
    config = make_config(
        arguments.chrom,
        recalculate_similarity=recalculate_similarity,
        workers=arguments.workers,
    )
    run_workflow(config)


if __name__ == "__main__":
    main()