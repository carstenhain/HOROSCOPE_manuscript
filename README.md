# HOROSCOPE Manuscript

Scripts and notebooks used for the manuscript *HOROSCOPE: Decoding human centromere architecture from short reads using k-mer signatures*.

## File Inventory

### CDR

| File | Description |
| --- | --- |
| `CDR/gao_hgsvc_cdr.ipynb` | Calculates CDR positions in HGSVC3 HOR arrays for haplogroup and cancer-breakpoint analyses. |

### Centromere annotation

| File | Description |
| --- | --- |
| `centromere_annotation/humas_annotation.py` | Annotates centromere FASTA sequences with HumAS and STV profiles and recolors HOR calls. |

### Centromere selection

| File | Description |
| --- | --- |
| `centromere_selection/get_centromeres_from_fasta.py` | Extracts centromere sequences from T2T genomes using flanking-region alignments. |

### K-mer generation

| File | Description |
| --- | --- |
| `kmer_generation/README.md` | Documents the k-mer generation and clustering workflows. |
| `kmer_generation/aggregate_kmers.py` | Combines k-mer counts from active centromeres on the same chromosome. |
| `kmer_generation/fasta_to_kmers.py` | Converts FASTA sequences into 61-mers with configurable counting modes. |
| `kmer_generation/make_kmer_blacklist.ipynb` | Builds a blacklist of non-centromeric 61-mers from the CHM13 reference. |

### Centromere clustering package

| File | Description |
| --- | --- |
| `kmer_generation/centromere_clustering/__init__.py` | Marks the centromere clustering workflow as a Python package. |
| `kmer_generation/centromere_clustering/clustering.py` | Calculates Dice similarities, hierarchical clusters, and dendrogram node selections. |
| `kmer_generation/centromere_clustering/config.py` | Defines chromosome-specific cluster identifiers and file paths. |
| `kmer_generation/centromere_clustering/data.py` | Loads active centromeres and constructs and filters k-mer count matrices. |
| `kmer_generation/centromere_clustering/length.py` | Predicts centromere and HOR lengths and finds normalization k-mers. |
| `kmer_generation/centromere_clustering/output.py` | Writes tagging-k-mer results and sample-cluster membership tables. |
| `kmer_generation/centromere_clustering/pipeline.py` | Provides the command-line entry point for the clustering workflow. |
| `kmer_generation/centromere_clustering/plots.py` | Generates dendrogram and k-mer count distribution plots. |
| `kmer_generation/centromere_clustering/tagging.py` | Finds cluster-tagging and relative k-mers from selected dendrogram nodes. |

### Model training

| File | Description |
| --- | --- |
| `model_training/test.py` | Empty placeholder for model-training tests or experiments. |