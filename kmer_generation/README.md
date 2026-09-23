# K-mer Generation and Clustering

## Generate k-mer tables

- `fasta_to_kmers.py`: split a FASTA sequence into k-mers with `dicey`.
- `aggregate_kmers.py`: combine counted k-mers from active centromeres of one chromosome into a count matrix.

## Cluster centromeres

Use `--recalculate-similarity` to rebuild the cached pairwise k-mer similarity matrix, and `--workers N` to set its worker count.

The previous monolithic scripts are now compatibility shims. The implementation is organized in `centromere_clustering/`:

- `config.py`: paths, thresholds, and curated chromosome cluster IDs.
- `data.py`: sample metadata, count-table loading, and k-mer filtering.
- `clustering.py`: pairwise similarity, hierarchical linkage, and selected tree nodes.
- `tagging.py`: cluster-tagging k-mer discovery.
- `length.py`: total centromere/HOR/HOR-span prediction and normalization k-mers.
- `plots.py`: clustering and tagging figures.
- `output.py`: result-table writers.
- `pipeline.py`: command-line orchestration.

Cluster-specific length k-mer discovery is intentionally not part of this workflow.
