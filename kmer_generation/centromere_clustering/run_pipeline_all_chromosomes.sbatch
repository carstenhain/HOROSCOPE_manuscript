#!/usr/bin/env bash
#SBATCH --job-name=centromere-clustering
#SBATCH --array=1-24
#SBATCH --cpus-per-task=10
#SBATCH --mem=150G
#SBATCH --time=20:00:00
#SBATCH --output=/scratch/hain/centromere/cluster_%j.out
#SBATCH --error=/scratch/hain/centromere/cluster_%j.err

set -euo pipefail

chromosomes=(
    chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12
    chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY
)

chromosome="${chromosomes[$((SLURM_ARRAY_TASK_ID - 1))]}"

python /g/korbel/hain/centromere/revision/HOROSCOPE_manuscript/kmer_generation/centromere_clustering/pipeline.py \
    --chrom "${chromosome}" \
    --workers "${SLURM_CPUS_PER_TASK}" \
    --recalculate-similarity
