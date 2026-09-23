import os
import argparse
import pandas as pd # type: ignore
import psutil # type: ignore
import numpy as np # type: ignore

WORKDIR = "/scratch/hain/centromere/"
SAMPLESHEET = "/g/korbel/hain/centromere/samplesheets/centromeres_full.tsv"

### only load kmers from the dicey output, do not load the kmer count
def load_kmers_only(path_to_dicey_list):
    kmers = set()
    for line in open(path_to_dicey_list):
        kmers.add(line.rstrip().lstrip().split(" ")[1])
    return kmers

### load kmer counts for a list of kmers (union)
def load_kmer_counts(path_to_dicey_list, union):
    tmp_dict = {}
    # Map k-mers to counts from one sample
    for line in open(path_to_dicey_list):
        tmp_dict[line.rstrip().lstrip().split(" ")[1]] = int(line.rstrip().lstrip().split(" ")[0])

    # Align kmer counts to the union order, use 0 for absent kmers
    tmp_list = []
    for kmer in union:
        try:
            tmp_list.append(tmp_dict[kmer])
        except KeyError:
            tmp_list.append(0)

    del tmp_dict

    return tmp_list

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Aggregate kmer counts for one chromosome")
    parser.add_argument("--chrom", type=str, help="Chromosome")
    # Parse the chromosome to aggregate from the command line.
    args = parser.parse_args()

    ### load samplesheet and subset to active chromosomes from the specified chromosome
    df = pd.read_csv(SAMPLESHEET, sep="\t")
    active_cent = df[(df["CHROM"] == args.chrom) & (df["ACTIVE"] == True)]

    log_fo = open(f"{WORKDIR}aggregate_{args.chrom}.log.txt", "w")

    # Initialize counters and the union set for k-mers.
    merge_counter = 0
    counter = 0
    union = set()
    
    # Collect unique k-mers in batches to limit intermediate memory usage.
    kmers = []
    for idx, row in active_cent.iterrows():

        ### append all kmers of a single sample
        kmers.append(load_kmers_only(f"{WORKDIR}kmers/{row['CHROM']}_{row['SAMPLE']}_{row['HP']}_{row['ASM']}_{row['PROJECT']}.dicey_unique.counts.list"))

        ### periodically merge kmers to limit memory usage
        merge_counter += 1
        if merge_counter == active_cent.shape[0] or merge_counter % 50 == 0:
            Z = union.union(*kmers)
            union = Z
            kmers = []

        counter += 1
        process = psutil.Process(os.getpid())
        memory_info = process.memory_info()
        memory_usage_in_MB = memory_info.rss / (1024 * 1024)
        # Log the progress and resident memory after each sample.
        log_fo.write(f"{counter}/{active_cent.shape[0]} ({round(memory_usage_in_MB / 1024, 1)} GB, {len(union)} kmers)\n")
        log_fo.flush()

    union = list(union)
    del kmers

    # Convert the final k-mer union to an ordered list for array columns.
    fo = open(f"{WORKDIR}count_arrays/{args.chrom}_kmer_index.list", "w")
    for kmer in union:
        fo.write(kmer + "\n")
    fo.close()

    data = []
    counter = 0

    # Build one count vector per active centromere sample.
    for idx, row in active_cent.iterrows():
        data.append(load_kmer_counts(
            f"{WORKDIR}kmers/{row['CHROM']}_{row['SAMPLE']}_{row['HP']}_{row['ASM']}_{row['PROJECT']}.dicey_unique.counts.list",
            union))
        counter += 1
        process = psutil.Process(os.getpid())
        memory_info = process.memory_info()
        # Log progress and memory while constructing the count matrix.
        memory_usage_in_MB = memory_info.rss / (1024 * 1024)
        log_fo.write(f"{counter}/{active_cent.shape[0]} ({round(memory_usage_in_MB / 1024, 1)} GB)\n")
        log_fo.flush()

    np_data = np.array(data)

    del data

    # Convert count vectors into a two-dimensional NumPy matrix.
    np.save(f"{WORKDIR}count_arrays/{args.chrom}_kmer_counts.npy", np_data)

    log_fo.close()

    sample_labels = [args.chrom + "_" + x for x in active_cent["ID"]]
    # Label matrix rows by sample and save a labeled, compressed HDF5 copy.
    df = pd.DataFrame(np_data, columns=union, index=sample_labels)
    df.to_hdf(f"{WORKDIR}count_arrays/{args.chrom}_kmer_counts.h5", mode="w", key=args.chrom, complevel=3)