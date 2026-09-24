import pandas as pd # type: ignore
import sys
import warnings
from pathlib import Path

PROJECT_ROOT = Path("/g/korbel/hain/centromere/revision/HOROSCOPE_manuscript")
CLUSTERING_SOURCE = PROJECT_ROOT / "kmer_generation" / "centromere_clustering"
sys.path.insert(0, str(CLUSTERING_SOURCE))

from config import CLUSTER_IDS_BY_CHROMOSOME, make_config # type: ignore
from length import correlation_analysis, count_hor, reverse_complement # type: ignore

warnings.filterwarnings("ignore")

def gather_cluster_and_hor_haplotype (active_centromeres:pd.DataFrame, samples:list, chrom:str, clustering_config):
    ### build table with cluster and centromere length per haplotype
    
    # empty data for dataframe
    sample_information_data = {"SAMPLE":[], "HP":[]}

    # initialize all haplotypes - two per sample
    for sample in samples:
        for i in range(0, 2):
            sample_information_data["SAMPLE"].append(sample)
            sample_information_data["HP"].append(f"H{1 + i}")

    # transform to dataframe
    sample_information = pd.DataFrame(data=sample_information_data)

    # prepare empty columns for cluster and length for each centromere
    sample_information[f"{chrom}_CLUSTER"] = "NA"
    sample_information[f"{chrom}_HOR_LENGTH"] = -1

    # load cluster assignments for this centromere and remove columns for clusters outside of the manually chose clusters
    chrom_cluster_assignments = pd.read_csv(f"/g/korbel/hain/centromere/tagging_kmers/{chrom}_sample_cluster_assignment.tsv", sep="\t", index_col="SAMPLE")[[f"{chrom}_{x}" for x in CLUSTER_IDS_BY_CHROMOSOME[chrom]]]

    # set cluster ids per haplotype
    for idx, row in chrom_cluster_assignments.iterrows():
        row_sample = idx.split("_")[0] # pyright: ignore[reportAttributeAccessIssue]
        row_hp = idx.split("_")[1].replace("P", "") # pyright: ignore[reportAttributeAccessIssue]
        row_cluster = "UNK"
        for cluster in chrom_cluster_assignments.columns:
            if row[cluster] == 1:
                row_cluster = cluster
                break
        sample_information.loc[(sample_information["SAMPLE"] == row_sample) & (sample_information["HP"] == row_hp), f"{chrom}_CLUSTER"] = row_cluster
    
    # load HOR size
    for _, row in active_centromeres.iterrows():
        row_sample = row['ID'].split("_")[0]
        row_hp = row['ID'].split("_")[1].replace("P", "")   
        hor_len, hor_span = count_hor(clustering_config, row['ID'])
        sample_information.loc[(sample_information["SAMPLE"] == row_sample) & (sample_information["HP"] == row_hp), f"{chrom}_HOR_LENGTH"] = hor_len
        
    return sample_information


def squash_sample_haplotypes(sample_information: pd.DataFrame, samples: list, chrom: str) -> pd.DataFrame:
    ### Squash haplotype-level information into sample-level information
    ### This function takes the haplotype-level data (output of gather_cluster_and_hor_haplotype) and computes mean lengths,
    ### mean HOR lengths and cluster IDs for each sample.
    ### Returns a DataFrame with one row per sample, contains only samples with two complete haplotypes.
    
    
    squashed_sample_information = pd.DataFrame(data={"SAMPLE":sample_information["SAMPLE"].unique()})
    squashed_sample_information.set_index("SAMPLE", inplace=True)

    # prepare columns
    squashed_sample_information[f"{chrom}_MEAN_HOR_LENGTH"] = -1.0
    squashed_sample_information[f"{chrom}_CLUSTER_ID"] = "NA"

    for sample in samples:
        # slice for this sample
        Z = sample_information[sample_information["SAMPLE"] == sample]

        # check if both haplotypes exist for this sample
        if not (Z[f"{chrom}_CLUSTER"] == "NA").any() and not (Z[f"{chrom}_CLUSTER"] == "UNK").any() and not (Z[f"{chrom}_HOR_LENGTH"] == -1).any():
            squashed_sample_information.loc[sample, f"{chrom}_MEAN_HOR_LENGTH"] = Z[f"{chrom}_HOR_LENGTH"].mean()
            try:
                cluster_numbers = sorted([int(x.split("_")[1]) for x in Z[f"{chrom}_CLUSTER"].tolist()])
                squashed_sample_information.loc[sample, f"{chrom}_CLUSTER_ID"] = "_".join([f"{chrom}_{str(x)}" for x in cluster_numbers])
            except Exception as e:
                print(f"Error processing sample {sample}: {e}")
                print(Z)
                break
        else:
            pass # this are the samples with only one complete centromere - all values are -1 for those samples

    print(f"Number of samples: {len(squashed_sample_information)}")
    print(f"Number of samples with complete information: {(squashed_sample_information[f'{chrom}_MEAN_HOR_LENGTH'] != -1).sum()}")

    return squashed_sample_information[squashed_sample_information[f'{chrom}_MEAN_HOR_LENGTH'] != -1] # type: ignore