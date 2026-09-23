"""Whole-centromere length and normalization k-mer analysis."""

import random
import subprocess
import tempfile

import numpy as np # type: ignore
import pandas as pd # type: ignore
import pysam # type: ignore
from sklearn.cluster import KMeans # type: ignore
from tqdm import tqdm # type: ignore

from tagging import empty_kmer_results, is_kmer


def correlation_analysis(first_values, second_values):
    """Return the absolute Pearson correlation, excluding zero-valued k-mers."""

    first_values = np.asarray(first_values)
    second_values = np.asarray(second_values)
    if np.amin(first_values) <= 0:
        return 0
    correlation = np.abs(np.corrcoef(first_values, second_values)[0, 1])
    return 0 if np.isnan(correlation) else correlation


def count_hor(config, sample_id):
    """Return live-HOR length and span for one centromere sample."""

    humas_directory = config.paths.hor_results_root / "humas_{}".format(
        config.chromosome
    ) / "humas"
    candidate_paths = (
        humas_directory / "AS-HOR-vs-{}_{}.bed".format(config.chromosome, sample_id),
        humas_directory
        / "AS-HOR-vs-{}_{}.contig_harmonized.bed".format(
            config.chromosome, sample_id
        ),
    )
    hor_path = next((path for path in candidate_paths if path.exists()), None)
    if hor_path is None:
        print("No HOR file found for sample {} {}".format(config.chromosome, sample_id))
        return 0, 0

    hor_calls = pd.read_csv(
        hor_path,
        sep="\t",
        names=["CONTIG", "START", "END", "ID", "A", "B", "C", "D", "E"],
        header=None,
    )
    live_hors = hor_calls[hor_calls["ID"].str.contains("L.", regex=False)].copy()
    if live_hors.empty:
        return 0, 0
    live_hors["LENGTH"] = live_hors["END"] - live_hors["START"]
    return live_hors["LENGTH"].sum(), live_hors["END"].max() - live_hors["END"].min()


def add_hor_measurements(config, active_centromeres):
    """Add live-HOR length and span to the active-centromere table."""

    measurements = [
        count_hor(config, sample_id) for sample_id in active_centromeres["ID"]
    ]
    centromeres = active_centromeres.copy()
    centromeres["HOR_LENGTH"] = [measurement[0] for measurement in measurements]
    centromeres["HOR_SPAN"] = [measurement[1] for measurement in measurements]
    return centromeres


def correlate_length_kmers(kmer_counts, length_vector, chromosome):
    """Find k-mers whose copy number correlates with a total-length measure."""
    """
        This function only caluclates correlation between single kmers and the target (e.g. the HOR size) 
        and saves these correlations in a dataframe for further usage.
    """
    
    correlations = []
    kmer_columns = [column for column in kmer_counts.columns if is_kmer(column)]
    for kmer in tqdm(kmer_columns, total=len(kmer_columns)):
        correlations.append(correlation_analysis(kmer_counts[kmer], length_vector))

    kmer_correlations = pd.DataFrame({"KMER": kmer_columns, "R": correlations})
    selected_kmers = kmer_correlations.loc[
        kmer_correlations["R"] > 0.75, "KMER"
    ].tolist() # type: ignore
    return pd.DataFrame(
        {
            "KMER": selected_kmers,
            "CLUSTER": "{}_LENGTH".format(chromosome),
            "MEAN_COUNT_IN": [kmer_counts[kmer].mean() for kmer in selected_kmers],
            "NUM_PRESENT_OUT": 0,
            "TYPE": "LENGTH",
        }
    )


def reverse_complement(sequence):
    """Return the reverse complement of an uppercase DNA sequence."""

    complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return "".join(complement[base] for base in sequence)[::-1]


def _run_command(command, stdout=None):
    subprocess.run(command, check=True, stdout=stdout)


def find_normalization_kmers(kmer_counts, config):
    """Find per-arm CN=1 normalization k-mers using perfect reference alignments."""

    candidate_kmers = kmer_counts.loc[:, (kmer_counts == 1).all()].columns.tolist()
    reference = config.paths.t2t_centromeres / "{}_CHM13Y_HP1_T2T_T2T.fa".format(
        config.chromosome
    )
    _run_command([str(config.paths.bwa), "index", str(reference)])

    with tempfile.TemporaryDirectory(
        dir=str(config.paths.normalization_workdir),
        prefix="{}_normalization_".format(config.chromosome),
    ) as temporary_directory:
        temporary_directory = config.paths.normalization_workdir.__class__(
            temporary_directory
        )
        kmers_fasta = temporary_directory / "kmers.fa"
        raw_sam = temporary_directory / "kmers.sam"
        sorted_sam = temporary_directory / "kmers.sorted.sam"
        filtered_bam = temporary_directory / "kmers.filtered.bam"

        with kmers_fasta.open("w") as output_file:
            for kmer in candidate_kmers:
                output_file.write(">{0}\n{0}\n".format(kmer))

        with raw_sam.open("w") as output_file:
            _run_command(
                [
                    str(config.paths.bwa),
                    "mem",
                    "-t",
                    "8",
                    "-a",
                    str(reference),
                    str(kmers_fasta),
                ],
                stdout=output_file,
            )
        _run_command(
            [str(config.paths.samtools), "sort", str(raw_sam), "-o", str(sorted_sam)]
        )

        with pysam.AlignmentFile(sorted_sam, "r") as input_file:
            with pysam.AlignmentFile(filtered_bam, "wb", template=input_file) as output_file:
                for alignment in input_file:
                    if (
                        alignment.has_tag("NM")
                        and alignment.get_tag("NM") == 0
                        and alignment.cigarstring is not None
                        and "S" not in alignment.cigarstring
                        and "H" not in alignment.cigarstring
                    ):
                        output_file.write(alignment)

        alignments = {"KMER": [], "POS": []}
        with pysam.AlignmentFile(filtered_bam, "rb") as alignment_file:
            for alignment in alignment_file:
                alignments["KMER"].append(alignment.query_name)
                alignments["POS"].append(alignment.reference_start)

    alignment_table = pd.DataFrame(alignments)
    if alignment_table.shape[0] < 1000:
        print("ERROR: Not enough potential normalization kmers found")
        return empty_kmer_results()

    classifications = KMeans(n_clusters=2, n_init="auto").fit_predict(
        alignment_table["POS"].values.reshape(-1, 1) # type: ignore
    )
    alignment_table["CLASSIFICATION"] = classifications
    first_mean = alignment_table.loc[
        alignment_table["CLASSIFICATION"] == 0, "POS"
    ].mean() # type: ignore
    second_mean = alignment_table.loc[
        alignment_table["CLASSIFICATION"] == 1, "POS"
    ].mean() # type: ignore
    p_label, q_label = (0, 1) if first_mean < second_mean else (1, 0)

    random_source = random.Random(42)
    result_tables = []
    for arm, label in (("P", p_label), ("Q", q_label)):
        arm_kmers = alignment_table.loc[
            alignment_table["CLASSIFICATION"] == label, "KMER"
        ].tolist() # type: ignore
        selected_kmers = random_source.sample(arm_kmers, min(5000, len(arm_kmers)))
        kmer_set = set(selected_kmers)
        kmer_set.update(reverse_complement(kmer) for kmer in selected_kmers)
        result_tables.append(
            pd.DataFrame(
                {
                    "KMER": list(kmer_set),
                    "CLUSTER": "{}_NORM_{}".format(config.chromosome, arm),
                    "MEAN_COUNT_IN": 1,
                    "NUM_PRESENT_OUT": 0,
                    "TYPE": "NORM",
                }
            )
        )
    return pd.concat(result_tables, ignore_index=True)