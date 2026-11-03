# Centromere annotation

This directory contains a python script to run humAS for annotating alpha-satellite higher-order repeats (HORs) in a centromere sequence and renames the contig and colors the HOR based on the number of monomers.

| File | Description |
| --- | --- |
| `humas_annotation.py` | Runs the HumAS-HMMER and STV annotation workflows on a centromere FASTA sequence, counts HOR monomers, assigns colours by monomer count, and writes a harmonized BED annotation. |
