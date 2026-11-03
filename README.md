# HOROSCOPE Manuscript

Scripts and notebooks used in the manuscript *HOROSCOPE: Decoding human centromere architecture from short reads using *k*-mer signatures*.

## Project structure

### [*k*-mer generation](./kmer_generation/README.md)

Builds and aggregates centromere *k*-mer count matrices, creates the non-centromeric *k*-mer blacklist, and clusters reference centromeres by *k*-mer similarity. The `centromere_clustering` package identifies architecture-specific tagging *k*-mers and predicts centromere and HOR lengths.

### [Centromere selection](./centromere_selection/README.md)

Extracts centromere sequences from T2T assemblies using alignments to flanking regions.

### [Centromere annotation](./centromere_annotation/README.md)

Runs HumAS and STV annotation on centromere sequences to identify alpha-satellite higher-order repeats. The workflow harmonizes contig names and assigns HOR colours according to monomer count.

### [Model training](./model_training/README.md)

Builds assembly-derived truth sets, normalizes *k*-mer counts, and trains models for centromere architecture and HOR-length prediction. It includes held-out-sample and mixed-versus-accurate assembly evaluations, with metrics recorded per chromosome and experiment.

### [Population analysis](./population_analysis/README.md)

Analyses HOROSCOPE centromere architecture predictions across 1000 Genomes Project and Human Genome Diversity Project populations. Computes architecture diversity, missingness, multiple factor analysis results, and population SNP heterozygosity.

### [ascairn comparison](./ascairn_comparison/README.md)

Runs ascairn on selected 1000 Genomes CRAM files through a Nextflow workflow. The analysis notebooks compare ascairn cluster and HOR-length calls with assembly-derived truth and HOROSCOPE predictions.

### [Cancer analysis](./cancer/README.md)

Aggregates and filters HOROSCOPE predictions from the PCAWG and PPCG cancer cohorts. It evaluates centromere changes in chromosome and chromosome-arm copy-number events, including relationships to centromere dip-region positions.

### [CDR](./CDR/README.md)

Calculates relative centromere dip-region (CDR) positions from HGSVC3 HOR annotations. The resulting CDR coordinates and metadata support haplogroup and cancer copy-number analyses.

### [Plotting](./plotting/)

Contains notebooks and generated manuscript figures. The figures summarize *k*-mer, population, ascairn-comparison, CDR, model-evaluation, and cancer-analysis results.

### [Train/test splits](./train_test_splits/)

Defines seeded 1kGP training and held-out sample sets for model development and evaluation. The accompanying notebook documents the strategy used to construct these splits.

### [Reviewer questions](./reviewer_questions/)

Contains analyses, input tables, and supplementary outputs prepared in response to reviewer requests. Subdirectories are organized by revision round and include experiments on assembly accuracy, leave-out evaluation, and *k*-mer origin.