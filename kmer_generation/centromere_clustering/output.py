"""Workflow result-table writers."""

import pandas as pd # type: ignore


def write_tagging_kmers(config, tagging_kmers):
    """Write discovered tagging, length, and normalization k-mers."""

    path = config.paths.tagging_kmers / "{}_tagging_kmers.tsv".format(config.chromosome)
    tagging_kmers.to_csv(path, sep="\t", index=False)


def write_cluster_assignments(config, active_centromeres, tagging_nodes):
    """Write binary sample membership for every selected tagging node."""

    assignments = pd.DataFrame({"SAMPLE": active_centromeres["ID"].tolist()})
    sample_ids = active_centromeres["ID"].tolist()
    for node in tagging_nodes:
        node_members = set(node.pre_order(lambda child: child.id))
        assignments["{}_{}".format(config.chromosome, node.id)] = [
            int(index in node_members) for index in range(len(sample_ids))
        ]
    path = config.paths.tagging_kmers / "{}_sample_cluster_assignment.tsv".format(
        config.chromosome
    )
    assignments.to_csv(path, sep="\t", index=False)