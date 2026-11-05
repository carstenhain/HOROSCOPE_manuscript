"""Discovery of cluster-tagging k-mers."""

import numpy as np # type: ignore
import pandas as pd # type: ignore


RESULT_COLUMNS = ["KMER", "CLUSTER", "MEAN_COUNT_IN", "NUM_PRESENT_OUT", "TYPE"]


def is_kmer(column_name):
    """Return whether a data-frame column name is a 61-mer sequence."""

    return isinstance(column_name, str) and len(column_name) == 61


def empty_kmer_results():
    """Create an empty result table with the workflow output schema."""

    return pd.DataFrame(columns=RESULT_COLUMNS)


def find_tagging_kmers(kmer_counts, tagging_nodes, chromosome):
    """Find unique and relative k-mers that label selected dendrogram nodes."""

    kmer_columns = [column for column in kmer_counts.columns if is_kmer(column)]
    tagging_kmers = empty_kmer_results()

    for node in tagging_nodes:
        inside = kmer_counts[kmer_counts[node.id] == 1]
        outside = kmer_counts[kmer_counts[node.id] == 0]

        print(
            "Gathering tagging kmers for node {} with {} samples".format(
                node.id, inside.shape[0]
            )
        )
        ### Compute statistics for k-mer presence inside and outside the cluster
        statistics = pd.DataFrame(
            {
                "Inside_Present": (inside[kmer_columns] != 0).sum(),
                "Outside_Present": (outside[kmer_columns] != 0).sum()
            }
        )
        
        inside_cutoff = 0.98 * inside.shape[0]             # set to used threshold, test: minimal influence
        outside_cutoff = 1 if inside.shape[0] <= 10 else 5 # set a stricter threshold for small clusters
        
        unique_kmers = statistics[
            (statistics["Inside_Present"] >= inside_cutoff) &
            (statistics["Outside_Present"] <= outside_cutoff)
        ]

        node_results = empty_kmer_results()
        for selected_kmers, label in (
            (unique_kmers, "UNIQUE"), # removed approach for RELATIVE kmers
        ):
            selected_columns = selected_kmers.index.tolist()
            candidate_results = pd.DataFrame(
                {
                    "KMER": selected_columns,
                    "CLUSTER": "{}_{}".format(chromosome, node.id),
                    "MEAN_COUNT_IN": inside[selected_columns].mean().values,
                    "NUM_PRESENT_OUT": outside[selected_columns].clip(upper=1).sum().values,
                    "TYPE": label,
                }
            )
            prefilter_count = candidate_results.shape[0]
            candidate_results = candidate_results.sort_values(
                by="MEAN_COUNT_IN", ascending=True
            ).iloc[: min(20000, unique_kmers.shape[0])]
            print(f"Filtering changed the number of kmers from {prefilter_count} to {candidate_results.shape[0]}")
            
            node_results = pd.concat([node_results, candidate_results], ignore_index=True) # type: ignore

        if node_results.empty:
            print("\tNo tagging kmers found for this cluster")
            continue

        maximum_count = node_results["MEAN_COUNT_IN"].max()
        highest_count_kmer = node_results.loc[
            node_results["MEAN_COUNT_IN"] == maximum_count, "KMER"
        ].iloc[0] # type: ignore
        unique_count = (node_results["TYPE"] == "UNIQUE").sum()
        print(
            f"\tFound {node_results.shape[0]} tagging kmers ({unique_count} UNIQUE) "
            f"(best kmer: {highest_count_kmer}, CN {maximum_count})"
        )
        tagging_kmers = pd.concat([tagging_kmers, node_results], ignore_index=True)

    return tagging_kmers