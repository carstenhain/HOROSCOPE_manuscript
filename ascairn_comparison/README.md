# ascairn comparison workflow

`main.nf` runs `ascairn type_all` once for every row in `ascairn_samples.csv`.
The sample sheet must have the columns `SAMPLE_NAME` and `ENA_FILE_PATH`.

## Requirements

Launch Nextflow from an environment where `ascairn` is available, for example:

```bash
conda activate ascairn

module load Nextflow
```

The reference FASTA and ascairn resource directory must be readable from the
compute environment. The defaults in `main.nf` match the paths used in
`ascairn_data.ipynb`.

## Run

From this directory:

```bash
nextflow run main.nf -resume
```

Results are copied to `results/<SAMPLE_NAME>/`. Work directories remain in
`work/`; use `-resume` to reuse completed tasks after an interruption.

## Slurm

Submit one Slurm job per sample with the `slurm` profile. Supply the partition
and account used at your site:

```bash
nextflow run main.nf -profile slurm -resume \
    --slurm_partition compute \
    --slurm_account YOUR_ACCOUNT
```

The profile requests 8 CPUs (from `--cpus`), 32 GB of memory, and 24 hours per
sample, and permits at most 50 concurrent tasks. Override the defaults where
necessary, for example `--cpus 16`. The
execution report, timeline, and task trace are written to `pipeline_info/`.

Override paths or compute settings at launch time:

```bash
nextflow run main.nf -resume \\
    --samples ascairn_samples.csv \\
    --outdir /scratch/hain/centromere/revision_ascairn_comparison/output \\
    --resource_dir /scratch/hain/ascairn/ascairn_resource/resource/panel/ascairn_paper_2025 \\
    --reference /scratch/hain/centromere/revision_ascairn_comparison/GRCh38_full_analysis_set_plus_decoy_hla.fa \\
    --cpus 8
```

## Files

| File | Description |
| --- | --- |
| `main.nf` | Defines the Nextflow workflow that runs `ascairn type_all` for each sample in `ascairn_samples.csv`. Slice the data from online crams. Publishes one `<SAMPLE_NAME>.cen_type_all.txt` result per sample to `results/`. |
| `nextflow.config` | Nextflow contig for running main.nf with the ascairn type_all. |
| `ascairn_samples.csv` | Lists the 1kGP samples submitted to the workflow and their ENA CRAM URLs. |
| `1000G_2504_high_coverage.sequence.index` | Is the 1000 Genomes high-coverage sequence index used to identify sample CRAM locations. `ascairn_data.ipynb` combines it with the 698_related_high_coverage samples. |
| `1000G_698_related_high_coverage.sequence.index` | Is the 1000 Genomes related-sample high-coverage sequence index. `ascairn_data.ipynb` combines it with the 504_high_coverage samples. |
| `ascairn_data.ipynb` | Selects samples with both HOROSCOPE data and 1kGP CRAM files, then writes `ascairn_samples.csv`. |
| `ascairn_analysis.ipynb` | Parses ascairn calls and compares their cluster and HOR-length predictions with assembly-derived truth and HOROSCOPE predictions. Calculates  metrics and generates the associated evaluation plots. |
| `results/` | Contains one `.cen_type_all.txt` file per processed sample. Each file records ascairn's two cluster calls, haplotype identifiers, and contig lengths for each chromosome. |