import os
import pickle
import numpy as np
import pandas as pd


def load_pickle(data_root):
    try:
        with open(os.path.join(data_root, 'Optimalres_0.5__OptimalGPTK_7.pkl'), 'rb') as ff:
            df_res05 = pickle.load(ff)
        with open(os.path.join(data_root, 'Optimalres_0.6__OptimalGPTK_7.pkl'), 'rb') as ff:
            df_res06 = pickle.load(ff)
        with open(os.path.join(data_root, 'Optimalres_0.7__OptimalGPTK_7.pkl'), 'rb') as ff:
            df_res07 = pickle.load(ff)
        with open(os.path.join(data_root, 'Optimalres_0.8__OptimalGPTK_7.pkl'), 'rb') as ff:
            df_res08 = pickle.load(ff)
        with open(os.path.join(data_root, 'Optimalres_0.9__OptimalGPTK_7.pkl'), 'rb') as ff:
            df_res09 = pickle.load(ff)
    except ModuleNotFoundError:
        df_res05 = pd.read_pickle(os.path.join(data_root, 'Optimalres_0.5__OptimalGPTK_7.pkl'))
        df_res06 = pd.read_pickle(os.path.join(data_root, 'Optimalres_0.6__OptimalGPTK_7.pkl'))
        df_res07 = pd.read_pickle(os.path.join(data_root, 'Optimalres_0.7__OptimalGPTK_7.pkl'))
        df_res08 = pd.read_pickle(os.path.join(data_root, 'Optimalres_0.8__OptimalGPTK_7.pkl'))
        df_res09 = pd.read_pickle(os.path.join(data_root, 'Optimalres_0.9__OptimalGPTK_7.pkl'))

    return df_res05, df_res06, df_res07, df_res08, df_res09


def load_molecular_subtypes(adata, data_root):
    df_res05, df_res06, df_res07, df_res08, df_res09 = load_pickle(data_root)
    subtypes_res05 = df_res05.iloc[0]['Molecular subtypes']
    subtypes_res06 = df_res06.iloc[0]['Molecular subtypes']
    subtypes_res07 = df_res07.iloc[0]['Molecular subtypes']
    subtypes_res08 = df_res08.iloc[0]['Molecular subtypes']
    subtypes_res09 = df_res09.iloc[0]['Molecular subtypes']

    assert np.all(df_res05.iloc[0]['MUC IDs'] == adata.obs.index), 'Index are not in the same order as in adata'

    adata.obs['Molecular Subtype res0.5'] = subtypes_res05
    adata.obs['Molecular Subtype res0.5'] = adata.obs['Molecular Subtype res0.5'].astype('category')
    adata.obs['Molecular Subtype res0.6'] = subtypes_res06
    adata.obs['Molecular Subtype res0.6'] = adata.obs['Molecular Subtype res0.6'].astype('category')
    adata.obs['Molecular Subtype res0.7'] = subtypes_res07
    adata.obs['Molecular Subtype res0.7'] = adata.obs['Molecular Subtype res0.7'].astype('category')
    adata.obs['Molecular Subtype res0.8'] = subtypes_res08
    adata.obs['Molecular Subtype res0.8'] = adata.obs['Molecular Subtype res0.8'].astype('category')
    adata.obs['Molecular Subtype res0.9'] = subtypes_res09
    adata.obs['Molecular Subtype res0.9'] = adata.obs['Molecular Subtype res0.9'].astype('category')

    # Add embedding
    with open(os.path.join(data_root, 'MolecularSubtypes_embedding.pkl'), 'rb') as ff:
        embedding_geneselection = pickle.load(ff)

    adata.obsm['X_geneselection_umap'] = np.asarray(embedding_geneselection)
    adata.var['geneselection'] = adata.var_names.isin(df_res09.iloc[0]['Gene names'])

    return adata


def reorder_mucids(adata, df_subtype):
    df_clusters = pd.DataFrame.from_dict({
        'MUC IDs': df_subtype.iloc[0]['MUC IDs'], 'Molecular subtypes': df_subtype.iloc[0]['Molecular subtypes']})
    index = list(adata.obs.index[adata.obs.index.isin(df_clusters["MUC IDs"])])
    df_clusters["MUC IDs"] = df_clusters["MUC IDs"].astype("category")
    df_clusters["MUC IDs"] = df_clusters["MUC IDs"].cat.set_categories(index)
    df_clusters = df_clusters.sort_values(["MUC IDs"])

    df_clusters.index = df_clusters["MUC IDs"]

    return df_clusters


def load_molecular_subtypes_lesion_nonlesion(adata, data_root):
    df_res05, df_res06, df_res07, df_res08, df_res09 = load_pickle(data_root)

    df_clusters_res05 = reorder_mucids(adata=adata, df_subtype=df_res05)
    df_clusters_res06 = reorder_mucids(adata=adata, df_subtype=df_res06)
    df_clusters_res07 = reorder_mucids(adata=adata, df_subtype=df_res07)
    df_clusters_res08 = reorder_mucids(adata=adata, df_subtype=df_res08)
    df_clusters_res09 = reorder_mucids(adata=adata, df_subtype=df_res09)

    assert np.all(df_clusters_res05['MUC IDs'] == adata.obs.index[
        adata.obs.index.isin(df_clusters_res05["MUC IDs"])]), 'Index are not in the same order as in adata'

    mask_res05 = adata.obs.index.isin(df_clusters_res05["MUC IDs"])
    adata.obs['Molecular Subtype res0.5'] = 'non lesional'
    adata.obs.loc[mask_res05, 'Molecular Subtype res0.5'] = df_clusters_res05['Molecular subtypes'].astype(str)
    adata.obs['Molecular Subtype res0.5'] = adata.obs['Molecular Subtype res0.5'].astype('category')

    mask_res06 = adata.obs.index.isin(df_clusters_res06["MUC IDs"])
    adata.obs['Molecular Subtype res0.6'] = 'non lesional'
    adata.obs.loc[mask_res06, 'Molecular Subtype res0.6'] = df_clusters_res06['Molecular subtypes'].astype(str)
    adata.obs['Molecular Subtype res0.6'] = adata.obs['Molecular Subtype res0.6'].astype('category')

    mask_res07 = adata.obs.index.isin(df_clusters_res07["MUC IDs"])
    adata.obs['Molecular Subtype res0.7'] = 'non lesional'
    adata.obs.loc[mask_res07, 'Molecular Subtype res0.7'] = df_clusters_res07['Molecular subtypes'].astype(str)
    adata.obs['Molecular Subtype res0.7'] = adata.obs['Molecular Subtype res0.7'].astype('category')

    mask_res08 = adata.obs.index.isin(df_clusters_res08["MUC IDs"])
    adata.obs['Molecular Subtype res0.8'] = 'non lesional'
    adata.obs.loc[mask_res08, 'Molecular Subtype res0.8'] = df_clusters_res08['Molecular subtypes'].astype(str)
    adata.obs['Molecular Subtype res0.8'] = adata.obs['Molecular Subtype res0.8'].astype('category')

    mask_res09 = adata.obs.index.isin(df_clusters_res09["MUC IDs"])
    adata.obs['Molecular Subtype res0.9'] = 'non lesional'
    adata.obs.loc[mask_res09, 'Molecular Subtype res0.9'] = df_clusters_res09['Molecular subtypes'].astype(str)
    adata.obs['Molecular Subtype res0.9'] = adata.obs['Molecular Subtype res0.9'].astype('category')

    adata.var['geneselection'] = adata.var_names.isin(df_res09.iloc[0]['Gene names'])

    return adata
