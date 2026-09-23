# HOROSCOPE_manuscript
Script used in the manuscript HOROSCOPE: Decoding human centromere architecture from short reads using k-mer signatures. Includes data generation, analysis and plotting scripts


### CDR
`gao_hgsvc_cdr.ipynb` Notebook for getting CDR for the HGSVC3 samples from the data provided in Gao et al and combine this with HOR annotations to caluclate the relative position of the CDR inside the HOR array to use (a) for haplogroup-specific CDR position and (b) for the inference of the correlation between CDR and breakpoint position in cancer.
### centromere_selection
`get_centromeres_from_fasta.py` Extracts centromeres from T2T genome
### kmer_generation
`fasta_to_kmers.py` Splits fasta file into *k*-mers using dicey, different modes either getting unique *k*-mers or counting unique *k*-mers or unique forward *k*-mers

`aggregate_kmers.py` Aggregate *k*-mers of multiple centromeres from the sample chromosome into a single *k*-mer table

`cluster_centromeres_v2.py` Runs the refactored centromere clustering, k-mer tagging, global length prediction, and normalization workflow. See [kmer_generation/README.md](kmer_generation/README.md) for usage and module layout.
### centromere_annotation
`humas_annotation` Annotates a centromere fasta with humAS and STV and recolors the HOR calls