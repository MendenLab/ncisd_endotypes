import matplotlib.pyplot as plt
import scanpy as sc
import os
import pandas as pd

import numpy as np
from datetime import date


def main(save_folder):
    file_name = '20210720_patient_meta_data_v04__CH'
    h5file_name = os.path.join('/Volumes', 'CH__data', 'Projects', 'Eyerich_AG_projects',
                               'BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer', 'analysis', 'Molecular_subtypes',
                               'input', 'h5_files', 'NON_LESION', 'LNL_RNAseq_{}.h5'.format(file_name))
    adata = sc.read(h5file_name)

    # add obs where lesion and non-lesion samples are annotated with same disease name
    adata.obs['disease'] = 'Unknown'
    disease_cat = list(adata.obs.diag.cat.categories)
    disease_cat.remove('non lesional')
    for diag in disease_cat:
        mask_diag = adata.obs['diag'] == diag
        id_index = list(adata.obs['Pseudo ID'][mask_diag])

        adata.obs['disease'][adata.obs['Pseudo ID'].isin(id_index)] = diag

    # Plot HVG genes in Heatmap of
    # TODO replace with disease discribing genes
    diag = 'eczema'

    adata_disease = adata[adata.obs['disease'] == diag].copy()

    var_value = np.var(adata_disease.X, axis=0)
    min_disp = np.min(var_value)  # 0.5
    max_disp = np.max(var_value)
    mean_value = np.mean(adata_disease.X, axis=0)
    min_mean = np.min(mean_value)  # 0.0125
    max_mean = np.max(mean_value)  # 3
    sc.pp.highly_variable_genes(adata_disease,  flavor='cell_ranger', batch_key='batchID', n_top_genes=10,
                                min_disp=min_disp, max_disp=max_disp, min_mean=min_mean, max_mean=max_mean)

    adata_pso_hvg = adata_disease[:, adata_disease.var['highly_variable']].copy()

    sc.pp.scale(adata_pso_hvg, zero_center=True)
    adata_pso_hvg.raw = adata_pso_hvg

    # Plot Heatmap with top 20 genes (of one patient: Pseudo ID 10656)
    # adata_hvg_patient = adata_pso_hvg[adata_pso_hvg.obs['Pseudo ID'] == '10417'].copy()

    sc.pl.heatmap(adata_pso_hvg, var_names=adata_disease.var_names[adata_disease.var['highly_variable']],
                  groupby='sampleType', use_raw=True, figsize=(6, 2), show=False, cmap='viridis')
    plt.ylabel(ylabel='biopsy type', fontsize=18)
    plt.yticks(fontsize=16)
    plt.xticks(fontsize=16)
    plt.savefig(os.path.join(save_folder, 'Figure_1A_GEx_heatmap.pdf'), bbox_inches='tight')
    plt.close('all')

    # Save data to excel file
    var_names = adata_disease.var_names[adata_disease.var['highly_variable']]
    obs_cols = ['sampleType']
    var_symbols = var_names
    var = adata_pso_hvg.raw.var
    dim_names = var.index

    X = adata_pso_hvg.raw.X
    idx = dim_names.get_indexer(list(var_names))
    mutable_idxer = [slice(None), slice(None)]
    mutable_idxer[1] = idx
    matrix = X[tuple(mutable_idxer)]

    # Create dataframe
    df = pd.DataFrame(index=adata_pso_hvg.obs_names)
    df = pd.concat(
        [df, pd.DataFrame(matrix, columns=var_symbols, index=adata_pso_hvg.obs_names)],
        axis=1,
    )
    df = pd.concat([df, adata_pso_hvg.obs[obs_cols]], axis=1)

    # sort by keys
    keys = ['sampleType'] + list(np.unique(var_names))
    df = df[keys]

    categorical = df['sampleType'].astype('category')
    categorical.name = 'sampleType'

    df = df[var_names].set_index(categorical)
    # Save to excel sheet
    df.to_excel(os.path.join(save_folder, 'Figure_1A_GEx_heatmap.xlsx'))


if __name__ == '__main__':
    output_dir = os.path.join('/Volumes', 'CH__data', 'Projects', 'Eyerich_AG_projects',
                              'BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer', 'analysis', 'Molecular_subtypes',
                              'output', 'Figure_1A_GEx', str(date.today()))
    os.makedirs(output_dir, exist_ok=True)

    main(save_folder=output_dir)
