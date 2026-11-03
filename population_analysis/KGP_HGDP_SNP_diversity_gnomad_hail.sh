#!/usr/bin/env bash
#SBATCH --job-name=shannon_hail
#SBATCH --cpus-per-task=1
#SBATCH --mem=500G
#SBATCH --time=20:00:00
#SBATCH --output=/scratch/hain/centromere/shannon_%j.out
#SBATCH --error=/scratch/hain/centromere/shannon_%j.err

python /g/korbel/hain/centromere/revision/HOROSCOPE_manuscript/population_analysis/KGP_HGDP_SNP_diversity_gnomad_hail.py