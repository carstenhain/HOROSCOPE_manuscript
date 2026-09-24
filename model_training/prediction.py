import pandas as pd # type: ignore

def count_present_kmer_fraction (normed_kmer_df:pd.DataFrame, samples:list, chrom:str, clusters:list, unique_cutoff=0.1):
    ### gather fraction of cluster unique kmers present per sample
    
    # initialize empty dataframe    
    result_df = pd.DataFrame(data={"SAMPLE":samples})
    
    for cluster in clusters:
                
        #print(f"Cluster {chrom}_{cluster} has {normed_kmer_df[normed_kmer_df['CLUSTER'] == f'{chrom}_{cluster}'].shape[0]} kmers")
        
        Z = normed_kmer_df[normed_kmer_df['CLUSTER'] == f'{chrom}_{cluster}'][samples]
        result_df[f"{chrom}_{cluster}"] = 1 - Z[Z > unique_cutoff].isna().sum().values / Z.shape[0] # type: ignore
        
    return result_df

def build_asm_prediction_dataframe (sample_hp_cluster_df, samples, clusters, chrom):
    # initialize empty dataframe
    asm_df = pd.DataFrame(data={"SAMPLE":samples})

    for cluster in clusters:
        
        col_cluster = []
        for sample in samples:
            if cluster in sample_hp_cluster_df.loc[sample, f"{chrom}_CLUSTER_ID"]:
                col_cluster.append(1)
            else:
                col_cluster.append(0)

        asm_df[f"{cluster}"] = col_cluster
        
    return asm_df