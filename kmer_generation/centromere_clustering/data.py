"""Data loading and k-mer filtering helpers."""

import numpy as np # type: ignore
import pandas as pd # type: ignore
from pathlib import Path
import subprocess

def load_active_centromeres(config):
    """Load active centromeres for the configured chromosome."""

    samples = pd.read_csv(config.paths.samplesheet, sep="\t")
    return samples[
        (samples["CHROM"] == config.chromosome) & (samples["ACTIVE"] == True)
    ].copy().reset_index(drop=True)


def load_population_colors(config):
    """Return superpopulation colours keyed by population name."""

    colors = pd.read_csv(config.paths.population_colors, sep="\t", comment="#")
    population_colors = {"UNK": [0, 0, 0]}
    for _, row in colors.iterrows():
        population_colors[row["population"]] = [ # type: ignore
            float(value) for value in row["rgb_rel"].split(",")
        ]
    return population_colors


def load_kmer_counts(config, active_centromeres):
    """Load the k-mer count array and label it with its k-mer sequences."""

    ### load kmers from hdf5 file and subset to active centromeres
    df = pd.read_hdf(
        config.kmer_hdf5,
        key=f"/{config.chromosome}"
    ).loc[[f"{config.chromosome}_{x}" for x in active_centromeres["ID"].tolist()]]
    return df
    
    """
    chromosome = config.chromosome
    data = np.load(config.paths.count_arrays / "{}_kmer_counts.npy".format(chromosome))
    index_path = config.paths.count_arrays / "{}_kmer_index.list".format(chromosome)
    with index_path.open() as index_file:
        kmers = [line.strip() for line in index_file]

    sample_labels = [
        "{}_{}".format(chromosome, sample_id)
        for sample_id in active_centromeres["ID"]
    ]
    counts = pd.DataFrame(data, columns=kmers, index=sample_labels)
    return counts.loc[:, (counts.shape[0] - (counts == 0).sum()) > 4]
    """

def filter_clustering_kmers(kmer_counts):
    """Remove rare, ubiquitous, and low-copy k-mers before clustering."""

    present_in = kmer_counts.shape[0] - (kmer_counts == 0).sum()
    return kmer_counts.loc[
        :, (present_in > 2)
        & (present_in < 0.9 * kmer_counts.shape[0])
        & (kmer_counts.max() > 2),
    ].copy()

def load_kmer_list (kmer_list_path):
    """
    Load a list of kmers from a file.

    Parameters
    ----------
    kmer_list_path : str
        Path to the file containing the list of kmers.

    Returns
    -------
    list of str
        List of kmers loaded from the file.
    """
    with open(kmer_list_path) as kmer_list_file:
        return [line.strip() for line in kmer_list_file]

def filter_kmers(all_kmers, config):
    """Filter kmers to remove those outside the centromere and present in other chromosomes."""
    
    kmer_list_path = Path(f"{config.paths.count_arrays}/{config.chromosome}_kmers.list")
    jellyfish_database = Path("/scratch/hain/centromere/kmer_resources/chm13_k61.jf")
    absent_kmers_path = Path(f"{config.paths.count_arrays}/{config.chromosome}_kmers_not_in_jf.list")

    with kmer_list_path.open("w") as kmer_list:
        for kmer in all_kmers.columns:
            kmer_list.write(f"{kmer}\n")

    command = (
        f"jellyfish query -i {jellyfish_database} < {kmer_list_path} "
        f"| paste {kmer_list_path} - "
        "| awk '$2 == 0 {print $1}' "
        f"> {absent_kmers_path}"
    )
    subprocess.run(command, shell=True, check=True)
    
    ### load kmers filtered to remove kmers on CHM13 outside the centromere
    filtered_kmers = load_kmer_list(absent_kmers_path)

    print(f"Initial set contained {all_kmers.shape[1]} kmers")
    print(f"Removal of genomic kmers from outside the centromere yielded {len(filtered_kmers)} centromeric kmers")

    ### remove kmers present in the centromeres of other chromosomes
    set_filtered_kmers = set(filtered_kmers)
    for chrom in [f"chr{x}" for x in range(1, 23)]:
        if chrom == config.chromosome:
            continue
        chrom_kmers = set(load_kmer_list(f"{config.paths.count_arrays}/{chrom}_kmers.list"))
        num_removed = len(set_filtered_kmers & chrom_kmers)
        set_filtered_kmers -= chrom_kmers
        print(f"Removal of centromeric kmers of {chrom}: removed {num_removed} kmers, total {len(set_filtered_kmers)}")
    ### save as list
    with open(f"{config.paths.count_arrays}/{config.chromosome}_kmer_index.filtered.list", "w") as f:
        for kmer in set_filtered_kmers:
            f.write(f"{kmer}\n")
            
    return f"{config.paths.count_arrays}/{config.chromosome}_kmer_index.filtered.list"
    

def load_permissive_kmers(config):
    """Load k-mers passing the external centromere-specificity filter."""

    path = config.paths.count_arrays / "{}_kmer_index.filtered.list".format(
        config.chromosome
    )
    return pd.read_csv(path, sep="\t", header=None, names=["KMER"])["KMER"]


def restrict_to_kmers(kmer_counts, kmers, haplotypes):
    """Restrict a count table to the provided k-mer sequences."""

    return kmer_counts.loc[
        kmer_counts.index.intersection(haplotypes), 
        kmer_counts.columns.intersection(kmers)
    ].copy()


def build_cluster_metadata(config, active_centromeres, population_colors):
    """Collect sample annotations and plotting colours alongside clustering data."""

    metadata = active_centromeres[
        ["SAMPLE", "HP", "ASM", "PROJECT", "ID", "SUPERPOP", "CONTIG_LENGTH"]
    ].copy()
    metadata["ID"] = config.chromosome + "_" + metadata["ID"]
    metadata["SUPERPOPCOL"] = [
        population_colors[population]
        for population in active_centromeres["SUPERPOP"]
    ]
    return metadata


def write_sample_order(config, active_centromeres, dendrogram_order):
    """Save the haplotype order produced by hierarchical clustering."""

    output_path = config.paths.tagging_kmers / "{}_sample_order.txt".format(
        config.chromosome
    )
    with output_path.open("w") as output_file:
        for sample_id in active_centromeres["ID"].iloc[dendrogram_order]:
            output_file.write("{}\n".format(sample_id))