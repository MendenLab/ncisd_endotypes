from scripts.feature_engineering import imputation, encoding
from sklearn import preprocessing

import os
from datetime import date

import scanpy as sc
import numpy as np
import pandas as pd
import itertools
import collections
import umap

import matplotlib.pyplot as plt
import plotly.graph_objects as go


class Hoverover_Plot:

    def __init__(self):
        # plot params
        self.dotsize = 10
        self.dotsize_responder_type = 14
        self.figure_size = (8, 8)
        self.xy_fontsize = 16
        self.xy_ticks = 12
        self.title_fontsize = 18
        self.legend_fontsize = 14
        self.text_fontsize = 12
        self.file_format = '.png'

    def plotly_interactive_umap(self, df, obs: str, save_folder: str, x_lab: str, y_lab: str, colors: list):
        """Plot interactive Volcano plot using plotly

        Parameters
        ----------
        colors
        df : pd.DataFrame
        obs : str
        save_folder : str
            path to save folder
        x_lab : str
            x-axis label
        y_lab : str
            y-axis label

        Returns
        -------

        """
        # Create custom text
        custom_text_mucid = list(df.index)
        custom_text_sdiag = df['sdiag'].values
        custom_text_diag = df['diag'].values
        custom_text_pattern = df['Pattern'].values
        custom_text = ["MUC ID: {}<br>Diagnosis: {}<br>Sub-diagnosis: {}<br>Pattern: {}".format(
            custom_text_mucid[val], custom_text_diag[val], custom_text_sdiag[val],
            custom_text_pattern[val]) for val in range(0, len(custom_text_diag))]

        df_plotly = pd.DataFrame.from_dict(
            {'UMAP1': df['x'], 'UMAP2': df['y'],
             'labels': df[obs], 'diag': df['diag'], 'Pattern': df['Pattern'],
             'text': custom_text})

        labels = list(df[obs].cat.categories)

        # plot
        fig = go.Figure()
        for ind_label, cluster_id in enumerate(labels):
            inds = np.where(np.asarray(labels) == cluster_id)[0]
            mask_cluster = df[obs] == cluster_id

            color = list(colors)

            # Plot samples without responder information
            fig.add_trace(
                go.Scatter(x=df_plotly['UMAP1'][mask_cluster], y=df_plotly['UMAP2'][mask_cluster],
                           mode='markers', marker=dict(size=self.dotsize, color=color[inds[0]]),
                           opacity=1, hoverinfo='text', name=cluster_id,
                           text=df_plotly.loc[mask_cluster, 'text'].values))  # hover text goes here,

        fig.update_layout(xaxis_title=x_lab, yaxis_title=y_lab, autosize=False, width=800, height=800)

        fig.write_html(os.path.join(
            save_folder, "{}_interactive.html".format(obs)))


def encode_impute_data(bulk_data, meta_data_legend, meta_data_data_type, save_folder):
    diagnosis_labels = bulk_data.obs['diag']

    dichotome_cols = meta_data_legend.loc[meta_data_legend['Explanation '].isin(
        meta_data_data_type['nominal']), 'Abbreviation'].tolist()
    # add manually
    dichotome_cols.extend(['Sex_x', 'Therapeutic target'])
    nominal_cols = [e for e in dichotome_cols if e in bulk_data.obs]  # 22
    # + Lab_totalIgE
    for for_col in nominal_cols:
        if bulk_data.obs[for_col].dtype == 'object':
            # Com_renal, Com_hepa, Com_ast, Com_dm, Com_rca, Com_iaa, Com_ah, Com_uv, Com_arths
            bulk_data.obs[for_col] = bulk_data.obs[for_col].astype(str).str.strip()
        if np.any(bulk_data.obs[for_col] == ''):
            bulk_data.obs.loc[bulk_data.obs[for_col] == '', for_col] = np.nan
        bulk_data.obs[for_col] = bulk_data.obs[for_col].astype(str).astype('category')
    # 2.2.1.2 Apply encoding -  dropped Com_uv as attribute was uniform
    df_encoded_nominals = encoding.encode_nominal_feature(
        data=bulk_data.obs[nominal_cols], dummy_variable_encoding=True)  # 43
    # 2.2.1.3 Handle missing values categories in encoded features
    df_encoded_nominals = encoding.imputation_replace_zero_nan(
        df_data=df_encoded_nominals, original_nominal_feature_names=nominal_cols)  # 21
    # set index
    df_encoded_nominals.index = bulk_data.obs.index
    print('Number of nominal categories: ', len(nominal_cols))  # 22
    print('Number of encoded nominal categories: ', df_encoded_nominals.shape[1])  # 21

    # 2.2.2 Ordinal Encoding
    ordinal_cols = meta_data_legend.loc[meta_data_legend['Explanation '].isin(
        meta_data_data_type['ordinal']), 'Abbreviation'].tolist()  # 26
    # add manually
    ordinal_cols.extend(['The_sor', 'The_dor', 'Hty_onset', 'Hty_pruritus'])
    # replace empty space '': MUC2870_Hist_Dist_Lympho
    for col in ordinal_cols:
        if bulk_data.obs[col].dtype == 'object':
            # print(col)
            bulk_data.obs[col] = bulk_data.obs[col].astype(str).str.strip()
        for row in bulk_data.obs.index:
            if bulk_data.obs.loc[row, col] == '':
                # print("{}_{}".format(row, col))
                bulk_data.obs.loc[row, col] = np.nan
    # remove columns which are not in bulk.obs
    ordinal_cols = [e for e in ordinal_cols if e in bulk_data.obs]  # 26
    # Apply ordinal encoding on the ordinal features
    df_encoded_ordinals = encoding.encode_ordinal_categories(data=bulk_data.obs[ordinal_cols].astype(float))
    # set index
    df_encoded_ordinals.index = bulk_data.obs.index
    print('Number of ordinal categories: ', len(ordinal_cols))  # 26
    print('Number of encoded ordinal categories: ', df_encoded_ordinals.shape[1])  # 26

    # 2.2.3 Treat continuous variable -> Decision 16.11.2021: leave them as they are
    continuous_cols = meta_data_legend.loc[meta_data_legend['Explanation '].isin(
        meta_data_data_type['continuous']), 'Abbreviation'].tolist()  # 21
    # add manually
    # ('Hist_Neutro_quant', 'Hist_Eosino' , 'Hist_ID_quant', 'Hist_LCV', 'Lab_Hb',  'Lab_leuco', 'Lab_granulo',
    #  'Lab_eosino', 'Lab_crea', 'Lab_GPT', 'Lab_totalIgE', 'Hty_chronic', 'The_leuco', 'The_Hb', 'The_lympho',
    #  'The_granulo', 'The_eosino', 'The_crea', 'The_GPT')
    continuous_cols.extend(['age', 'Hist_Cap', 'Hist_Exocytosis', 'Hist_Lympho', 'Hist_Mucin',
                            'Hist_Ortho_quant', 'Hist_Para_quant',
                            'Hist_Aka_quant', 'Hist_Spongiosis',
                            'The_leuco', 'The_Hb', 'The_lympho', 'The_granulo', 'The_eosino',
                            'The_crea', 'The_GPT', 'The_weeks'])
    continuous_cols = list(np.unique(continuous_cols)) # 32
    # remove columns which are not in bulk.obs
    continuous_cols = [e for e in continuous_cols if e in bulk_data.obs]  # 32
    for col in continuous_cols:
        if bulk_data.obs[col].dtype == 'object':
            bulk_data.obs[col] = bulk_data.obs[col].astype(str).str.strip()
        for row in bulk_data.obs.index:
            if bulk_data.obs.loc[row, col] == '':
                bulk_data.obs.loc[row, col] = np.nan
    # remove categories from these columns
    df_continous = bulk_data.obs[continuous_cols].copy()
    # where is this value: '>5000' -> MUC2711
    df_continous.loc[df_continous['Lab_totalIgE'] == '>5000', 'Lab_totalIgE'] = 5000
    # MUC7393: '0-1'
    df_continous.loc[df_continous['The_weeks'] == '0-1', 'The_weeks'] = 1
    df_continous = df_continous.astype('float')
    print('Number of continuous variables: ', df_continous.shape[1])  # 32

    # 2.2.5 Combine features again
    # Check if I capture all columns
    check_cols = list(itertools.chain(nominal_cols, ordinal_cols, continuous_cols))  # 86
    diff_adata_cols = [item for item, count in collections.Counter(check_cols).items() if count > 1]
    assert len(diff_adata_cols) == 0, \
        'Caution: your number of features to encode is not the same as your in your adata object: '.format(
            diff_adata_cols)

    # Concatenate all encoded features 292 rows x 92 columns
    df_encoded_features = pd.concat([df_encoded_nominals, df_encoded_ordinals, df_continous], axis=1)

    # 2.2.6 Remove categories with single category value such as 'Com_uv'
    drop_cols = [temp_col for temp_col in df_encoded_features.columns
                 if len(df_encoded_features[temp_col].dropna().unique()) <= 1]
    # df_encoded_features = feature_engineering.remove_lowvariance_features(df_data=df_encoded_features)
    df_encoded_features = df_encoded_features.drop(drop_cols, axis=1)

    # 2.5 Handle missing values
    # https://www.sciencedirect.com/science/article/pii/S0895435618308710
    # 2.5.1 Get number of complete cases (samples without missing values)
    df_complete_samples = df_encoded_features.copy()
    df_complete_samples.dropna(axis='rows', inplace=True)  # 8
    cc_diseases = []
    for ind_ in df_complete_samples.index:
        cc_diseases.append(diagnosis_labels[diagnosis_labels.index == ind_].values[0])
        print(diagnosis_labels[diagnosis_labels.index == ind_])

    # Proportion of complete cases: 23.36% == 8 samples
    prop_cs = df_encoded_features.shape[0] * df_complete_samples.shape[0] / 100
    print('Percentage of complete cases: ', prop_cs)
    # 2.5.2 Get number of features without missing values
    df_complete_features = df_encoded_features.copy()
    df_complete_features.dropna(axis='columns', inplace=True)  # 1 (age)
    # Proportion of complete features: 0.92% == 1 feature (age)
    prop_cf = df_encoded_features.shape[1] * df_complete_features.shape[1] / 100
    print('Percentage of complete features: ', prop_cf)

    # Impute missing features
    fig, ax = plt.subplots(figsize=(8, 8))
    df_encoded_features.isnull().sum().plot.hist(ax=ax, bins=12)
    ax.set_xlabel('Number of missing observations')
    plt.tight_layout()
    fig.savefig(os.path.join(save_folder, 'Frequency_missing_values.png'))
    plt.close(fig=fig)
    print("count of NULL values before imputation\n")
    print(df_encoded_features.isnull().sum())

    return df_encoded_features, nominal_cols, ordinal_cols, continuous_cols


def prepare_attributes(adata, save_folder):
    meta_data_legend = pd.read_excel(
        os.path.join(
            '/Volumes/CH__data/Projects/Eyerich_AG_projects/BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer',
            'raw_data', 'clinical_data', "patient_meta_data_final.xlsx"), sheet_name='legend')

    meta_data_data_type = pd.read_excel(
        os.path.join(
            '/Volumes/CH__data/Projects/Eyerich_AG_projects/BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer',
            'raw_data', 'clinical_data', "Tables.xlsx"), sheet_name='Table_S1')

    # Encode data
    df_encoded, nominal_cols, ordinal_cols, continuous_cols = encode_impute_data(
        bulk_data=adata, meta_data_legend=meta_data_legend,
        meta_data_data_type=meta_data_data_type, save_folder=save_folder)

    print('Number of variables after filtering: ', df_encoded.shape[1])

    # Drop columns with more than 70% missing values
    df_percentage = df_encoded.isnull().mean() * 100
    variables_to_drop = list(df_percentage[df_percentage > 70].index)
    print('Attributes to drop due to missing values > 70%: ', variables_to_drop)
    df_encoded = df_encoded.drop(variables_to_drop, axis=1)
    continuous_cols = [col for col in continuous_cols if col not in variables_to_drop]
    nominal_cols = [col for col in nominal_cols if col not in variables_to_drop]
    ordinal_cols = [col for col in ordinal_cols if col not in variables_to_drop]

    # Normalise ordinal and continuous data before imputation as the distance
    # measure in KNN imputation suffers from large values
    scaler = preprocessing.MinMaxScaler()
    ord_continuous_cols = continuous_cols + ordinal_cols
    df_ord_continuous_encoded_normed = pd.DataFrame(
        scaler.fit_transform(df_encoded.loc[:, ord_continuous_cols]),
        columns=ord_continuous_cols, index=df_encoded.index)
    df_encoded.loc[:, ord_continuous_cols] = df_ord_continuous_encoded_normed

    # Impute data - apply KNN only on numeric data (continuous variables)
    # set k=1 as we would get otherwise values inbetween categories such as 1.5 instead of 1 or 2
    df = imputation.impute_knn(func_data_numeric=df_encoded, ordinal_categories=ordinal_cols, k=1)

    return df


def get_embedding_clinical_attributes(scaled_data):
    """

    Parameters
    ----------
    scaled_data : pandas.Dataframe
    save_folder : str
    labels : pandas.Series or numpy.array or list
    data_type :
    mapping:

    Returns
    -------

    """
    # todo optimise min_dist and n_neighbors
    reducer = umap.UMAP(random_state=42)
    embedding = reducer.fit_transform(scaled_data)

    return embedding


def main(input_dir, save_folder):
    print('Start script to create hoverover plot for GEx embedding')
    bulk_adata = sc.read(
        os.path.join(input_dir, 
                     'Lesion_RNAseq_20210720_patient_meta_data_v04__CH_20220209_TherapyResponse_v04_v220715_TherapyResponseCorrected__Endotypes_230620.h5'))

    # Add correct clinical attributes
    df_corrected_metadata = pd.read_excel(os.path.join(
        '/Volumes', 'CH__data', 'Projects', 'Eyerich_AG_projects', 'BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer',
        'raw_data', 'clinical_data', 'patient_meta_data_final.xlsx'))
    # rename Sex.x, Helmholtz.identifyer', Diag_diag.x
    df_corrected_metadata = df_corrected_metadata.rename(
        columns={'Helmholtz.identifyer': 'Helmholtz_identifyer', 'Diag_diag.x': 'Diag_diag_x', 'Sex.x': 'Sex_x'})
    df_corrected_metadata.index = df_corrected_metadata['Helmholtz_identifyer']

    # common columns
    common_cols = np.intersect1d(bulk_adata.obs.columns, df_corrected_metadata.columns)
    assert np.all(df_corrected_metadata.loc[bulk_adata.obs.index, common_cols].index == bulk_adata.obs.index)
    # TODO might have to adjust some variables -> Peter
    bulk_adata.obs.loc[:, common_cols] = df_corrected_metadata.loc[bulk_adata.obs.index, common_cols]
    # rename row undefined to UD in Pattern
    bulk_adata.obs = bulk_adata.obs.replace({"Pattern": {"undefined": "UD"}})

    labels = pd.DataFrame(columns=['Pattern'], data=bulk_adata.obs['Pattern'])
    rename_patterns = {"Pattern": {"1": 1, "2a": 2, "2b": 3, "3": 4, "4a": 5, "4b": 6, "5": 7, "UD": 8}}
    labels = labels.replace(rename_patterns)
    labels['Pattern'] = labels['Pattern'].astype('category')

    # 2.2.8 Label Encoding
    label_train_encoded = encoding.encode_label(label=labels['Pattern'].values)
    labels = label_train_encoded.transform(labels.values)

    # Attribute preparation
    df = prepare_attributes(adata=bulk_adata, save_folder=save_folder)
    # get embeddings
    umap_emb = get_embedding_clinical_attributes(scaled_data=df)

    df_umap = pd.DataFrame.from_dict({
        'x': umap_emb[:, 0], 'y': umap_emb[:, 1], 'diag': bulk_adata.obs['diag'], 'sdiag': bulk_adata.obs['sdiag'],
        'Pattern': labels, 'Endotypes': bulk_adata.obs['Endotypes']})
    df_umap.index = bulk_adata.obs.index

    plot_obj = Hoverover_Plot()
    plot_obj.plotly_interactive_umap(
        df=df_umap, obs='Endotypes', save_folder=save_folder, x_lab="UMAP1", y_lab="UMAP2",
        colors=bulk_adata.uns['Endotypes_colors'])


if __name__ == '__main__':
    output_dir = os.path.join('/Volumes', 'CH__data', 'Projects', 'Eyerich_AG_projects',
                              'BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer', 'analysis', 'Molecular_subtypes',
                              'output', 'Figure_UMAP_Hoverover_clinical_attributes', str(date.today()))
    os.makedirs(output_dir, exist_ok=True)

    data_root = os.path.join(
        '/Volumes', 'CH__data', 'Projects', 'Eyerich_AG_projects', 'BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer',
        'analysis', 'Molecular_subtypes', 'input', 'h5_files', 'LESION')

    main(save_folder=output_dir, input_dir=data_root)
