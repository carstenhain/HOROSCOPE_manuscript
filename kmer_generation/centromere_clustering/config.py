"""Configuration for centromere k-mer clustering."""

from pathlib import Path


CLUSTER_IDS_BY_CHROMOSOME = {
    "chr1": (694, 698, 681, 670),
    "chr2": (929, 934, 922),
    "chr3": (999, 1009, 1033, 1064, 1061, 1060, 1051, 1029, 994, 1026, 1048),
    "chr4": (1075, 1079, 1072, 1058, 1081),
    "chr5": (1106, 1107, 1103, 991, 1098),
    "chr6": (1080, 1084, 1081, 1077, 1090, 1075, 1088),
    "chr7": (916, 918, 898, 910, 893, 830, 747, 905),
    "chr8": (1115, 1105, 1098, 1114, 1086, 1100, 1075, 1071, 1102),
    "chr9": (969, 1099, 1112, 992, 983, 938, 1097, 914, 1102, 1101, 1085, 1100, 1110),
    "chr10": (1114, 1116, 968, 1102, 1111, 1110, 1088, 1108),
    "chr11": (1104, 879, 1106, 1099, 1102, 1089, 906, 1094, 1095, 1082, 1085),
    "chr12": (1115, 1072, 1102, 1006, 1076, 1087, 1088, 1105, 1096, 1082, 1089, 1084),
    "chr13": (1025, 1016, 1000, 1009, 981, 1001, 883, 1004, 1002),
    "chr14": (726, 729, 727, 712, 710, 643, 660, 696),
    "chr15": (513, 591, 379, 584, 597),
    "chr16": (1030, 1051, 1075, 1070, 1068, 1071, 1012, 1052, 1074, 1047, 1009),
    "chr17": (1091, 1047, 1089, 1071, 1017, 1076, 1088),
    "chr18": (1098, 755, 1089, 1095, 1084, 1090, 1092),
    "chr19": (1087, 1090, 1088, 1094, 1101, 1098, 912, 1099, 1069, 1080, 1055, 1026, 1107, 1071, 1073, 1110),
    "chr20": (1081, 1091),
    "chr21": (968, 965),
    "chr22": (1031, 1026, 962, 976, 1006, 1027, 1029, 968, 1035, 876, 1036),
    "chrX": (882, 876, 881),
    "chrY": (245, 260, 249, 264),
}


class WorkflowPaths:
    """Filesystem locations used by the workflow."""

    def __init__(self):
        self.samplesheet = Path("/g/korbel/hain/centromere/samplesheets/centromeres_full.tsv")
        self.population_colors = Path("/g/korbel/hain/centromere/hgsvc3_resources/hgsvc3_pop_colors.tsv")
        self.count_arrays = Path("/scratch/hain/centromere/count_arrays")
        self.figures = Path("/g/korbel/hain/centromere/figures_new")
        self.tagging_kmers = Path("/g/korbel/hain/centromere/tagging_kmers_new")
        self.hor_results_root = Path("/scratch/hain")
        self.normalization_workdir = Path("/scratch/hain/centromere")
        self.t2t_centromeres = Path("/scratch/hain/centromere/centromeres")
        self.bwa = Path("/g/korbel/hain/conda-envs/jupyter_wave/bin/bwa")
        self.samtools = Path("/g/korbel/hain/conda-envs/jupyter_wave/bin/samtools")


class WorkflowConfig:
    """Parameters for one chromosome clustering run."""

    def __init__(
        self,
        chromosome,
        paths,
        cluster_ids,
        recalculate_similarity,
        kmer_hdf5,
        workers=16,
        tree_min=1.2,
        tree_max=10.0,
    ):
        self.chromosome = chromosome
        self.paths = paths
        self.cluster_ids = cluster_ids
        self.recalculate_similarity = recalculate_similarity
        self.kmer_hdf5 = kmer_hdf5
        self.workers = workers
        self.tree_min = tree_min
        self.tree_max = tree_max


def make_config(
    chromosome: str,
    *,
    paths=None,
    recalculate_similarity=None,
    kmer_hdf5=None,
    workers: int = 16,
) -> WorkflowConfig:
    """Build a configuration while preserving the old per-chromosome defaults.

    ``kmer_hdf5`` may be a string or ``Path`` pointing to the chromosome k-mer
    count table. When omitted, the standard count-array location is used.
    """

    cluster_ids = CLUSTER_IDS_BY_CHROMOSOME.get(chromosome, (1000,))
    if recalculate_similarity is None:
        recalculate_similarity = chromosome not in CLUSTER_IDS_BY_CHROMOSOME
    paths = paths or WorkflowPaths()
    if kmer_hdf5 is None:
        kmer_hdf5 = paths.count_arrays / "{}_kmer_counts.non_rare.h5".format(
            chromosome
        )

    return WorkflowConfig(
        chromosome=chromosome,
        paths=paths,
        cluster_ids=cluster_ids,
        recalculate_similarity=recalculate_similarity,
        kmer_hdf5=Path(kmer_hdf5),
        workers=workers,
    )