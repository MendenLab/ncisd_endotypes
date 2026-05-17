from analysis.utils import add_colors, add_endotypes, normalisation

import numpy as np
import scanpy as sc
import pandas as pd
from datetime import date

import scipy.stats as stats

import os

import rpy2.robjects.packages as rpackages
from rpy2.robjects import numpy2ri, pandas2ri

numpy2ri.activate()
pandas2ri.activate()

base = rpackages.importr("base")
ape = rpackages.importr("ape")
factoextra = rpackages.importr("factoextra")
dplyr = rpackages.importr('dplyr')
ggtree = rpackages.importr('ggtree')
ggplot2 = rpackages.importr('ggplot2')

fileformat = '.pdf'


def main(save_folder, input_dir):
    adata = sc.read(
        os.path.join(input_dir,
                     'Lesion_RNAseq_20210720_patient_meta_data_v04__CH_20220209_TherapyResponse_v04_v220715.h5'))
    adata = add_endotypes.load_molecular_subtypes(adata=adata, data_root=input_dir)

    # Sample: PseudoID 10350, MUC4259
    # new diag: lichenoid drug reaction
    # new pattern: pattern 1
    adata.obs.loc[adata.obs['diag'] == 'undefined', 'Pattern'] = '1'
    adata.obs['diag'] = adata.obs['diag'].astype(str)
    adata.obs.loc[adata.obs['diag'] == 'undefined', 'diag'] = 'lichenoid drug reaction'
    adata.obs['diag'] = adata.obs['diag'].astype('category')

    adata, _ = add_colors.diag_order_lesion(adata=adata)
    adata, _ = add_colors.sdiag_order_lesion(adata=adata)
    adata, _ = add_colors.pattern_order_lesion(adata=adata)
    adata, _ = add_colors.endotype_order_lesion(adata=adata, obs_name='Molecular Subtype res0.9')

    merge_methods = ['mean', 'median']

    # Read out dataframe
    # Should we use only our 7% of genes?
    df = adata.to_df(layer='counts')[list(adata.var_names[adata.var['geneselection']])]
    # Apply normalisation and batch correction
    array_gex = normalisation.edger_normalise(X=np.asarray(df), coldata=adata.obs, batch_keys=['batchID', 'Sex_x'])
    df_normed = pd.DataFrame(array_gex, columns=adata.var_names[adata.var['geneselection']],
                             index=adata.obs.index)
    # Scale data
    df_normed_zscores = stats.zscore(df_normed, axis=0)
    adata.obs['Molecular Subtype res0.9'] = adata.obs['Molecular Subtype res0.9'].astype(int)
    df_normed_zscores['cluster'] = adata.obs['Molecular Subtype res0.9'] + 1
    df_normed_zscores['cluster'] = df_normed_zscores['cluster'].astype(str)

    df_normed_zscores['cluster'] = df_normed_zscores['cluster'].replace(
        {'1': 'E12', '2': 'E7', '3': 'E13', "4": 'E3', "5": 'E1', "6": 'E9', "7": 'E8',
         "8": 'E10', "9": 'E11', "10": 'E5', "11": 'E6', "12": 'E2', "13": 'E4'})

    for merge in merge_methods:
        # Merge samples of a cluster into one column by using mean, median, ect
        if merge == 'mean':
            df_merged = df_normed_zscores.groupby(['cluster']).mean()  # Used in manuscript
        elif merge == 'median':
            df_merged = df_normed_zscores.groupby(['cluster']).median()

        df_merged.to_csv(os.path.join(save_folder, 'Normed_Dataframe_{}_gene_selection.csv'.format(merge)))


if __name__ == '__main__':
    output_dir = os.path.join('/Volumes/CH__data/Projects/Eyerich_AG_projects',
                              'BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer/analysis/output',
                              'Dendrogram', str(date.today()))
    os.makedirs(output_dir, exist_ok=True)

    data_root = os.path.join('/Volumes/CH__data/Projects/Eyerich_AG_projects',
                             'BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer/data/h5_files/LESION')

    main(save_folder=output_dir, input_dir=data_root)
