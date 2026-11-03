# *k*-mer generation and clustering

This directory contains the workflow that derives centromere *k*-mer count
tables, clusters centromere haplotypes, and identifies diagnostic *k*-mers.
Ignored files and directories are intentionally not described here.

| File | Description |
| --- | --- |
| `aggregate_kmers.py` | Loads counted Dicey *k*-mers for active centromeres of one requested chromosome and creates their union. It writes the *k*-mer index plus numpy and HDF5 count matrices. |
| `fasta_to_kmers.py` | Provides `fasta_to_kmers`, which invokes Dicey to generate forward *k*-mers and, when requested, reverse-complement *k*-mers from a FASTA file. It extracts FASTQ sequence records, sorts and deduplicates them with optional counts, and removes intermediate files. |
| `make_kmer_blacklist.ipynb` | Masks centromere intervals in the CHM13 reference genome using Bedtools. It then builds a Jellyfish database of 61-mers outside centromeres for use as a non-centromeric *k*-mer blacklist. |
| `centromere_clustering/clustering.py` | Calculates or loads cached pairwise similarity matrices from *k*-mer presence data. It builds weighted hierarchical linkage trees, selects configured internal nodes, and adds binary membership columns for those nodes. |
| `centromere_clustering/config.py` | Defines curated dendrogram node IDs for each chromosome and the filesystem locations used by the workflow. |
| `centromere_clustering/convert_into_tsv.sh` | Is a SLURM job that extracts selected length-associated *k*-mer columns from chromosome HDF5 count tables. It writes a tab-separated length-*k*-mer table for each configured chromosome. |
| `centromere_clustering/data.py` | Loads active-centromere metadata, population colours, and chromosome *k*-mer count tables. It filters *k*-mers, removes sequences found outside the target centromere or on other chromosomes, builds plotting metadata, and writes dendrogram sample order. |
| `centromere_clustering/length.py` | Measures live HOR length and span from annotation BED files, then identifies *k*-mers correlated with those measurements. It also finds copy-number-one normalization *k*-mers by aligning them to the T2T centromere reference and separating the two chromosome arms. |
| `centromere_clustering/output.py` | Writes the discovered tagging, length, and normalization *k*-mers as a chromosome-specific TSV file. It also writes a table assigning each active sample to every selected dendrogram node. |
| `centromere_clustering/pipeline.py` | Runs the main per-chromosome workflow: load and filter counts, cluster haplotypes, identify cluster-tagging, HOR-length, and normalization *k*-mers, and generate figures and outputs *k*-mers. Can split the "training" data into train and test based on a list of samples, currently hard-coded to use one seed of the 75:25 split (/g/korbel/hain/centromere/revision/HOROSCOPE_manuscript/test_samples_1kGP_seed_10.txt). |
| `centromere_clustering/plots.py` | Generates *k*-mer count-distribution plots, dendrograms with similarity and metadata tracks, and selected-node diagrams. It also creates heatmaps showing the highest-count tagging *k*-mers across the clustering order. |
| `centromere_clustering/run_pipeline_all_chromosomes.sbatch` | Submits an array job for the 24 autosomes and sex chromosomes. Each task runs `pipeline.py` for one chromosome and recalculates its similarity matrix. |
| `centromere_clustering/tagging.py` | Identifies 61-mers that are prevalent within selected dendrogram nodes and rare outside them. It labels candidates as unique or relative tagging *k*-mers and returns their mean within-cluster count and out-of-cluster prevalence. |
