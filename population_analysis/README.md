# Population analysis

This directory contains the population-level analyses of centromere architecture predictions from the 1000 Genomes Project (1kGP) and the Human Genome Diversity Project (HGDP).

| File | Description |
| --- | --- |
| `analyze_population_predictions.ipynb` | Loads 1kGP and HGDP entries from `complete_population_predictions.tsv` and joins them with population metadata. It calculates architecture entropy and cluster allele frequencies, assesses prediction missingness, performs multiple factor analysis (MFA), and generates the associated plots. |
| `complete_population_predictions.tsv` | Contains sample-level centromere architecture predictions, cluster assignments, coverage metrics, and population annotations. |
| `gr.bed` | Defines chromosome-wide GRCh38 intervals for the chromosomes included in the population SNP-count export. The Hail script processes each interval independently and writes per-chromosome population AC/AN tables. |
| `igsr_populations.tsv` | Provides IGSR population names, descriptions, geographic coordinates, superpopulation assignments, display colours, and source collections. |
| `KGP_HGDP_SNP_diversity.ipynb` | Calculates biallelic SNP heterozygosity as `2 * AF * (1 - AF)` from per-population alternate-allele count (AC) and allele-number (AN) tables. It summarizes mean heterozygosity for each population and chromosome. |
| `KGP_HGDP_SNP_diversity_gnomad_hail.py` | Reads the gnomAD v3.1.2 HGDP/1KGP subset Hail table and exports AC and AN values for every population to per-chromosome TSV files in `/scratch/hain/pop_freq/`. |
| `KGP_HGDP_SNP_diversity_gnomad_hail.sh` | Slurm submission script for `KGP_HGDP_SNP_diversity_gnomad_hail.py`. It requests one CPU, 500 GB of memory. |
