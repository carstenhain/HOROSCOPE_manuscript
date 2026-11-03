# Cancer analysis

This directory contains HOROSCOPE-based centromere analyses for Pan-Cancer Analysis of Whole Genomes (PCAWG) and pheochromocytoma/paraganglioma (PPCG) tumor samples. The workflow aggregates cancer predictions, applies copy-number-based filters, and assesses centromere changes associated with chromosome-arm events.

| File | Description |
| --- | --- |
| `aggregate_horoscope_cancer_predictions.ipynb` | Joins HOROSCOPE predictions with cancer metadata for the PCAWG and PPCG cohorts. Prepares the cohort-level prediction tables used by the downstream filtering and arm-level analyses. |
| `pcawg_filtering.ipynb` | Filters PCAWG predictions, including the removal of whole-genome-duplication samples and diagnostics of copy-neutral tumor/germline agreement. |
| `ppcg_filtering.ipynb` | Applies the analogous quality-control and centromere-number calculations to PPCG predictions. |
| `pcawg_arm_level_analysis.ipynb` | Classifies PCAWG chromosome and arm copy-number events and generates event landscapes. Fits Gaussian mixture models, incorporates CDR-position annotations, and examines event-specific centromere changes. |
| `ppcg_arm_level_analysis.ipynb` | Performs the same analysis for the PPCG cohort. |
| `Gao_CDR_position.tsv` | Provides chromosome- and copy-number-event-specific CDR position summaries . The arm-level notebooks use these annotations to stratify centromere-loss analyses. |