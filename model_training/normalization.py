import numpy as np # type: ignore
import pandas as pd # type: ignore

def normalize (kmer_df, chrom, mode, samples):
    
    # get normalization value for each column
    if mode == "ALL":
        normalization = kmer_df[kmer_df["CLUSTER"].isin([f"{chrom}_NORM_Q", f"{chrom}_NORM_P"])][samples].median()
    elif mode == "P":
        normalization = kmer_df[kmer_df["CLUSTER"].isin([f"{chrom}_NORM_P"])][samples].median()
    elif mode == "Q":
        normalization = kmer_df[kmer_df["CLUSTER"].isin([f"{chrom}_NORM_Q"])][samples].median()
    elif mode == "UNITIG":
        ### placeholder for chromosomes without normalization kmers
        normalization = pd.Series(23 * np.ones(len(samples)), index=samples)
    
    # apply normalization to each KM column
    out_df = pd.concat([kmer_df[["KMER", "CLUSTER"]], kmer_df[samples] / normalization], axis=1)
    
    return out_df