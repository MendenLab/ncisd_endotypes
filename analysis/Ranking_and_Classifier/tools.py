import pandas as pd 
import scanpy
import numpy as np


def load_brain( cfg, norm='DESeq2' ):

    adata = scanpy.read( cfg.data_path )

    index = adata.to_df().index

    if norm == 'DESeq2' :
        df = adata.to_df()
    elif norm == 'TMM' :
        df = pd.read_csv('./data-cohort/brain_kiel/Normed_counts_BRAIN_Kiel_TMM.csv',
                         index_col=0)
        df = df.loc[index]

    pattern = adata.obs["Pattern"].to_numpy()
    endotypes = adata.obs["Endotypes"].to_numpy()

    aa = pd.DataFrame( index=index )
    aa['pattern'] = pattern 
    aa['endotypes'] = endotypes

    return df, aa

