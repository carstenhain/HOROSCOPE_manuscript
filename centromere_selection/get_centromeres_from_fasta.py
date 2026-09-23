import pysam # type: ignore
import os
import subprocess
import pandas as pd # type: ignore
import argparse

T2T_INDEX_PATH = "/scratch/hain/centromere/mmi/t2t.asm5.mmi"


def run_cmd(cmd):
	# Run a shell command and fail fast if it exits with a non-zero status.
	subprocess.run(cmd, shell=True, check=True)


def remove_file_if_exists(path):
	# Delete temporary files only when present to avoid cleanup errors.
	if os.path.exists(path):
		os.remove(path)

def get_centromeres (fasta_file_path, flank_fasta_path, centromeric_fasta_path, id):
	### map flanking regions to assembly
	run_cmd(f"minimap2 -ax asm20 -t 10 --secondary=no {fasta_file_path} {flank_fasta_path} | samtools view -bh - | samtools sort - > tmp.{id}.flanks_to_asm.bam")
	run_cmd(f"samtools index tmp.{id}.flanks_to_asm.bam")


	### gathering flank alignment informations
	# Open the flank-alignment BAM and stop early when no alignments were produced.
	try:
		samfile = pysam.AlignmentFile(f"tmp.{id}.flanks_to_asm.bam", "rb")
	except ValueError:
		print(f"Not alignment found for sample {id} - Exiting")
		remove_file_if_exists(f"tmp.{id}.flanks_to_asm.bam")
		remove_file_if_exists(f"tmp.{id}.flanks_to_asm.bam.bai")
		return -1

	chromosome = []
	flank = []
	refname = []
	refstart = []
	refend = []
	reflen = []
	directions = []
	alnlen = []

	mapping_length_cutoff = 180000

	### find alignments to flanks
	# Keep only long flank mappings and collect coordinates/orientation metadata.
	for aln in samfile:

		if (aln.reference_end - aln.reference_start) > mapping_length_cutoff:

			direction = "FOR"
			if aln.is_reverse:
				direction = "REV"

			ref_length = 0
			for contig_name, contig_length in zip(samfile.references, samfile.lengths):
				if contig_name == aln.reference_name:
					ref_length = contig_length

			chromosome.append(aln.query_name.split("_")[0])
			flank.append(aln.query_name.split("_")[1])
			refname.append(aln.reference_name)
			refstart.append(aln.reference_start)
			refend.append(aln.reference_end)
			reflen.append(ref_length)
			directions.append(direction)
			alnlen.append(aln.reference_end - aln.reference_start)

	# Store flank mappings in a dataframe for per-chromosome interval inference.
	df = pd.DataFrame(data={"CHROM":chromosome, "FLANK":flank, "REFNAME":refname, "REFSTART":refstart, "REFEND":refend, "REFLEN":reflen, "DIR":directions, "ALNLEN":alnlen})

	intervals = []

	chromosomes = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]

	# Build centromeric intervals by matching left/right flanks and orientation rules.
	for chrom in chromosomes:
		print(f"Searching for {chrom} centromere")
		chrom_df = df[df["CHROM"] == chrom]
		if chrom_df.shape[0] == 2 and "right" in chrom_df['FLANK'].values and "left" in chrom_df['FLANK'].values:
			#print("Both flanks found")
			### check for same contig and same direction
			if chrom_df["REFNAME"].iloc[0] == chrom_df["REFNAME"].iloc[1] and chrom_df["DIR"].iloc[0] == chrom_df["DIR"].iloc[1]:
				print("\tSame chrom and same direction")
				if chrom_df["DIR"].iloc[0] == "FOR":
					intervals.append([chrom_df["REFNAME"].iloc[0], str(chrom_df[chrom_df["FLANK"] == "left"]["REFSTART"].iloc[0]), str(chrom_df[chrom_df["FLANK"] == "right"]["REFEND"].iloc[0]), id + "_" + chrom_df["CHROM"].iloc[0]])
				else:
					intervals.append([chrom_df["REFNAME"].iloc[0], str(chrom_df[chrom_df["FLANK"] == "right"]["REFSTART"].iloc[0]), str(chrom_df[chrom_df["FLANK"] == "left"]["REFEND"].iloc[0]), id + "_" + chrom_df["CHROM"].iloc[0]])
			else:
				if chrom_df["REFNAME"].iloc[0] == chrom_df["REFNAME"].iloc[1] and chrom_df["DIR"].iloc[0] != chrom_df["DIR"].iloc[1]:
					print("\tSame chrom but different direction")
				if chrom_df["REFNAME"].iloc[0] != chrom_df["REFNAME"].iloc[1]:
					print("\tBroken centromere")
					if chrom_df[chrom_df["FLANK"] == "left"]["DIR"].iloc[0] == "FOR":
						intervals.append([chrom_df[chrom_df["FLANK"] == "left"]["REFNAME"].iloc[0], str(chrom_df[chrom_df["FLANK"] == "left"]["REFSTART"].iloc[0]), str(chrom_df[chrom_df["FLANK"] == "left"]["REFLEN"].iloc[0]), id + "_" + chrom_df["CHROM"].iloc[0] + "_left"])
					if chrom_df[chrom_df["FLANK"] == "left"]["DIR"].iloc[0] == "REV":
						intervals.append([chrom_df[chrom_df["FLANK"] == "left"]["REFNAME"].iloc[0], "0", str(chrom_df[chrom_df["FLANK"] == "left"]["REFEND"].iloc[0]), id + "_" + chrom_df["CHROM"].iloc[0] + "_left"])
					if chrom_df[chrom_df["FLANK"] == "right"]["DIR"].iloc[0] == "FOR":
						intervals.append([chrom_df[chrom_df["FLANK"] == "right"]["REFNAME"].iloc[0], "0", str(chrom_df[chrom_df["FLANK"] == "right"]["REFEND"].iloc[0]), id + "_" + chrom_df["CHROM"].iloc[0] + "_right"])
					if chrom_df[chrom_df["FLANK"] == "right"]["DIR"].iloc[0] == "REV":
						intervals.append([chrom_df[chrom_df["FLANK"] == "right"]["REFNAME"].iloc[0], str(chrom_df[chrom_df["FLANK"] == "right"]["REFSTART"].iloc[0]), str(chrom_df[chrom_df["FLANK"] == "right"]["REFLEN"].iloc[0]), id + "_" + chrom_df["CHROM"].iloc[0] + "_right"])
		else:
			if chrom_df.shape[0] == 1:
				print(f"\tOnly one flank ({chrom_df['FLANK'].tolist()[0]}) found")
			else:
				print(f"\tNo flank found")

	# Write inferred centromeric intervals as BED for sequence extraction.
	fo = open(f"tmp.{id}.centromeric.bed", "w")
	for i in intervals:
		fo.write("\t".join(i) + "\n")
	fo.close()

	# Extract interval sequences from assembly FASTA (handling gzip input when needed).
	if fasta_file_path.endswith(".gz"):
		run_cmd(f"gzip -dkf {fasta_file_path}")
		run_cmd(f"samtools faidx {fasta_file_path.replace('.gz', '')}")
		run_cmd(f"bedtools getfasta -nameOnly -fi {fasta_file_path.replace('.gz', '')} -bed tmp.{id}.centromeric.bed > {centromeric_fasta_path}")
		remove_file_if_exists(fasta_file_path.replace('.gz', ''))
	else:
		run_cmd(f"bedtools getfasta -nameOnly -fi {fasta_file_path} -bed tmp.{id}.centromeric.bed > {centromeric_fasta_path}")

	### map combined fasta to t2t
	run_cmd(f"minimap2 -ax asm5 -t 10 --secondary=no {T2T_INDEX_PATH} {centromeric_fasta_path} | samtools view -bh -F 256 -F 2048 - | samtools sort - > tmp.{id}.cent_to_t2t.bam")

	### gather all centromere with reverse mapped primary alignment and store this information
	# Identify extracted centromeres that align to T2T on reverse strand.
	contigs_to_reverse = []
	for aln in pysam.AlignmentFile(f"tmp.{id}.cent_to_t2t.bam", "rb"):
		query_seq = aln.query_name.split("_")[-1]
		if query_seq in ["left", "right"]:
			query_seq = aln.query_name.split("_")[-2]
		ref_seq = aln.reference_name
		if query_seq == ref_seq:
			if not aln.is_reverse:
				print(f"Centromere cutout {aln.query_name} maps to same chromosomes in forward direction")
			else:
				print(f"Centromere cutout {aln.query_name} maps to same chromosomes in reverse direction")
				contigs_to_reverse.append(aln.query_name)
		else:
			print(f"Centromere cutout maps to different chromosomes (Cutout: {query_seq} | Reference: {ref_seq})")

	### split into single centromere fasta files, reverse indicated centromere
	run_cmd(f"samtools faidx {centromeric_fasta_path}")
	fasta = pysam.FastaFile(centromeric_fasta_path)
	sequences = {}
	# Reverse-complement contigs flagged by T2T mapping; keep others unchanged.
	for reference in fasta.references:
		if reference in contigs_to_reverse:
			complement = str.maketrans('ATCGN', 'TAGCN')
			sequences[reference] = fasta.fetch(reference).upper().translate(complement)[::-1]
		else:
			sequences[reference] = fasta.fetch(reference).upper()

	# Emit per-chromosome FASTA files containing all matching centromeric contigs.
	for chrom in chromosomes:
		this_chrom_sequences = []
		for key in sequences:
			if key.endswith(chrom):
				this_chrom_sequences.append(key)
			if chrom + "_left" in key:
				this_chrom_sequences.append(key)
			if chrom + "_right" in key:
				this_chrom_sequences.append(key)
		if len(this_chrom_sequences) > 0:
			X = centromeric_fasta_path.split("/")
			X[-1] = f"{chrom}_{X[-1]}"
			output_file = "/".join(X)
			with open(output_file, 'w') as f:
				for contig in this_chrom_sequences:
					f.write(f">{contig}\n")
					f.write(f"{sequences[contig]}\n")

	### rewrite the centromeric.bed file and include forward and reverse annotation (relative to T2T)
	# Save final BED with strand annotation relative to T2T orientation.
	cent_bed_fo = open(f"{id}.centromeric.bed", "w")
	for line in open(f"tmp.{id}.centromeric.bed"):
		if line.rstrip().split("\t")[-1] in contigs_to_reverse:
			cent_bed_fo.write(line.rstrip() + "\t1\t-\n")
		else:
			cent_bed_fo.write(line.rstrip() + "\t1\t+\n")
	cent_bed_fo.close()
	
	# Remove temporary intermediates generated during mapping and extraction.
	remove_file_if_exists(centromeric_fasta_path)
	remove_file_if_exists(f"{centromeric_fasta_path}.fai")
	remove_file_if_exists(f"tmp.{id}.flanks_to_asm.bam")
	remove_file_if_exists(f"tmp.{id}.flanks_to_asm.bam.bai")
	remove_file_if_exists(f"tmp.{id}.cent_to_t2t.bam")
	remove_file_if_exists(f"tmp.{id}.centromeric.bed")