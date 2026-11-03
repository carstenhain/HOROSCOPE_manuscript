# Model training

This directory contains the code and notebooks used to train and evaluate centromere-cluster and higher-order-repeat (HOR) length predictions from *k*-mer counts.

| File | Description |
| --- | --- |
| `build_datasheets.py` | Builds haplotype-level reference tables by joining centromere records with cluster assignments and HOR lengths. It then aggregates complete haplotype pairs to sample-level architecture calls and mean HOR lengths for model training. |
| `normalization.py` | Normalizes sample *k*-mer counts using chromosome-specific normalization *k*-mers (or a fixed unitig value; not implemented). It supports combined, P-only, Q-only, and unitig normalization modes. |
| `prediction.py` | Calculates the fraction of detected cluster-specific *k*-mers for each sample. It also converts assembly-derived diploid cluster calls into binary labels used to evaluate the predictions. |
| `model_training.ipynb` | Evaluates cluster-presence and trains HOR-length models for a selected chromosome using a predefined held-out 1kGP sample set. It normalizes *k*-mer counts, evaluates cluster calls against assemblies, selects length-correlated *k*-mers, and fits PCA plus linear regression models. |
