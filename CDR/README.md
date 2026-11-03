# CDR data

This directory contains a notebook and resources used to calculate and analyse centromere dip region (CDR) positions:

| File | Description |
| --- | --- |
| `centromere dip region annotation.bed` | CDR intervals in centromere cutout coordinates. The three columns are the cutout identifier, CDR start, and CDR end. From Gao et. al. 2026, Nature. |
| `alpha satellite higher order repeat annotation` | Folder with humAS HOR annotations of centromeres (not included). From Gao et. al. 2026, Nature. |
| `gao_hgsvc_cdr.ipynb` | Notebook that derives relative CDR positions from the CDR intervals, joins sample, haplotype, cluster, and superpopulation metadata, and produces analyses and plots. |
| `rel_cdr_position_full.tsv` | Derived relative CDR positions with chromosome, cluster, and superpopulation annotations. |