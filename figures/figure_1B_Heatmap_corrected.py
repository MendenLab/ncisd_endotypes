from scripts.feature_engineering import encoding, imputation
from scripts.utils import add_colors
from sklearn import preprocessing

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import scanpy as sc
import numpy as np
import os

import pandas as pd
from datetime import date

import itertools
import collections

import seaborn as sns


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

    # Concatenate all encoded features: 342 x 88 (24.03.2024)
    df_encoded_features = pd.concat([df_encoded_nominals, df_encoded_ordinals, df_continous], axis=1)

    # 2.2.6 Remove categories with single category value such as 'Com_uv' == zero variance
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

    # Plot missing features
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
        bulk_data=adata, meta_data_legend=meta_data_legend, meta_data_data_type=meta_data_data_type,
        save_folder=save_folder)

    # save used category names in dataframe
    df = pd.DataFrame({'nominal': pd.Series(nominal_cols), 'ordinal': pd.Series(ordinal_cols),
                       'continuous': pd.Series(continuous_cols)})
    df.to_excel(os.path.join(save_folder, 'Figure_1B_used_attributes.xlsx'))

    print('Number of variables after filtering: ', df_encoded.shape[1])

    # # Drop columns with more than 70% missing values
    # df_percentage = df_encoded.isnull().mean() * 100
    # variables_to_drop = list(df_percentage[df_percentage > 70].index)
    # print('Attributes to drop due to missing values > 70%: ', variables_to_drop)
    # df_encoded = df_encoded.drop(variables_to_drop, axis=1)
    # continuous_cols = [col for col in continuous_cols if col not in variables_to_drop]
    # nominal_cols = [col for col in nominal_cols if col not in variables_to_drop]
    # ordinal_cols = [col for col in ordinal_cols if col not in variables_to_drop]

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


def impute_data(bulk_adata):
    df_features = bulk_adata.obs[
        ['age', 'Sex_x',
         'Com_arths', 'Com_uv', 'Com_ah', 'Com_iaa', 'Com_rca', 'Com_dm', 'Com_ast', 'Com_hepa', 'Com_renal',
         'Com_tonsil', 'Com_BMI', 'Com_ibd', 'Com_nail', 'Com_muc', 'Com_scalp',
         'Morph_erythema', 'Morph_Elev', 'Morph_papu', 'Morph_vesicles', 'Morph_pustules', 'Morph_scaling',
         'Morph_dryness',
         'Hty_onset', 'Hty_fh', 'Hty_pruritus', 'Hty_exanth', 'Hty_chronic', 'Hty_self', 'Hty_photos',
         'Hty_progr', 'Hty_sport', 'Hty_alc', 'Hty_smoker',
         'Lab_csa', 'Lab_leuco', 'Lab_Hb', 'Lab_lympho', 'Lab_granulo', 'Lab_eosino', 'Lab_crea', 'Lab_GPT',
         'Lab_totalIgE', 'Lab_specificIgE',
         'The_dor', 'The_sor', 'The_leuco', 'The_Hb', 'The_lympho', 'The_granulo', 'The_eosino', 'The_crea',
         'The_GPT', 'The_weeks',
         'Hist_SE', 'Hist_Hyper_quant', 'Hist_Para_quant', 'Hist_Para_quali', 'Hist_Ortho_quant', 'Hist_Aka_quant',
         'Hist_Aka_quali', 'Hist_Granu', 'Hist_Serum', 'Hist_Spongiosis', 'Hist_Cap', 'Hist_Lympho', 'Hist_Dist_Lympho',
         'Hist_Exocytosis', 'Hist_Microabscess', 'Hist_Bacteria', 'Hist_Neutro_quant', 'Hist_Neutro_quali',
         'Hist_Eosino', 'Hist_ID_quant', 'Hist_Keratoses', 'Hist_Nr_Keratoses', 'Hist_ID_distribution', 'Hist_ID_quali',
         'Hist_LCV', 'Hist_Mucin']]

    df_features_knn_imputation = imputation.impute_knn(func_data_numeric=df_features, k=1)

    return df_features_knn_imputation


def separated_heatmap(df, df_colors, save_folder, figsize=(6, 12), cmap='RdYlBu'):

    fig = plt.figure(figsize=figsize)  # (w, h)

    # Split heatmap into clinical attribute groups
    # 0.74 / 89 = 0.008, use 0.006 * rows, space between attr groups: 0.004
    # 0. Age
    df_age_sex = df.loc[['age'], :]  # 1 rows
    axmatrix = fig.add_axes([0.2, 0.738, 0.65, 0.006])  # (left, bottom, width, height)
    sns.heatmap(df_age_sex,
                vmin=0, xticklabels=False, yticklabels=True, ax=axmatrix,
                cmap=cmap, center=0.5, cbar=False)
    axmatrix.set_yticklabels(axmatrix.get_yticklabels(), rotation=0)
    axmatrix.set_xticklabels(axmatrix.get_xticklabels(), rotation=90)
    axmatrix.yaxis.set_label_position("right")
    axmatrix.yaxis.tick_right()
    # axmatrix.set_yticklabels([])  # xlabels[:6], rotation=0
    axmatrix.set_xticklabels([])
    axmatrix.set_xlabel('')
    axmatrix.spines['top'].set_visible(True)
    axmatrix.spines['right'].set_visible(True)
    axmatrix.spines['bottom'].set_visible(True)
    axmatrix.spines['left'].set_visible(True)

    # 1. Sex
    df_age_sex = df.loc[['Sex'], :]  # 1 rows
    axmatrix = fig.add_axes([0.2, 0.728, 0.65, 0.006])  # (left, bottom, width, height)
    sns.heatmap(df_age_sex,
                vmin=0, xticklabels=False, yticklabels=True, ax=axmatrix,
                cmap=cmap, center=0.5, cbar=False)
    axmatrix.set_yticklabels(axmatrix.get_yticklabels(), rotation=0)
    axmatrix.set_xticklabels(axmatrix.get_xticklabels(), rotation=90)
    axmatrix.yaxis.set_label_position("right")
    axmatrix.yaxis.tick_right()
    # axmatrix.set_yticklabels([])  # xlabels[:6], rotation=0
    axmatrix.set_xticklabels([])
    axmatrix.set_xlabel('')
    axmatrix.spines['top'].set_visible(True)
    axmatrix.spines['right'].set_visible(True)
    axmatrix.spines['bottom'].set_visible(True)
    axmatrix.spines['left'].set_visible(True)

    # 2. Morphology
    df_morphology = df.loc[list(df.index[df.index.str.contains('Morph_')]), :]  # 7 rows
    axmatrix = fig.add_axes([0.2, 0.682, 0.65, 0.042])  # (left, bottom, width, height)
    sns.heatmap(df_morphology, vmin=0, xticklabels=False, yticklabels=True, ax=axmatrix,
                cmap=cmap, center=0.5, cbar=False)
    axmatrix.set_yticklabels(axmatrix.get_yticklabels(), rotation=0)
    axmatrix.set_xticklabels(axmatrix.get_xticklabels(), rotation=90)
    axmatrix.yaxis.set_label_position("right")
    axmatrix.yaxis.tick_right()
    axmatrix.set_xticklabels([])
    axmatrix.set_xlabel('')
    axmatrix.spines['top'].set_visible(True)
    axmatrix.spines['right'].set_visible(True)
    axmatrix.spines['bottom'].set_visible(True)
    axmatrix.spines['left'].set_visible(True)

    # 3. Severity
    df_severity = df.loc[list(df.index[df.index.str.contains('Diag_')]), :]  # 4 rows
    axmatrix = fig.add_axes([0.2, 0.654, 0.65, 0.024])  # (left, bottom, width, height)
    sns.heatmap(df_severity, vmin=0, xticklabels=False, yticklabels=True, ax=axmatrix,
                cmap=cmap, center=0.5, cbar=False)
    axmatrix.set_yticklabels(axmatrix.get_yticklabels(), rotation=0)
    axmatrix.set_xticklabels(axmatrix.get_xticklabels(), rotation=90)
    axmatrix.yaxis.set_label_position("right")
    axmatrix.yaxis.tick_right()
    axmatrix.set_xticklabels([])
    axmatrix.set_xlabel('')
    axmatrix.spines['top'].set_visible(True)
    axmatrix.spines['right'].set_visible(True)
    axmatrix.spines['bottom'].set_visible(True)
    axmatrix.spines['left'].set_visible(True)

    # 4. Therapy
    df_therapy = df.loc[list(df.index[df.index.str.contains('The_|Therapeutic ')]), :]  # 14 rows
    axmatrix = fig.add_axes([0.2, 0.566, 0.65, 0.084])  # (left, bottom, width, height)
    sns.heatmap(df_therapy, vmin=0, xticklabels=False, yticklabels=True, ax=axmatrix,
                cmap=cmap, center=0.5, cbar=False)
    axmatrix.set_yticklabels(axmatrix.get_yticklabels(), rotation=0)
    axmatrix.set_xticklabels(axmatrix.get_xticklabels(), rotation=90)
    axmatrix.yaxis.set_label_position("right")
    axmatrix.yaxis.tick_right()
    axmatrix.set_xticklabels([])
    axmatrix.set_xlabel('')
    axmatrix.spines['top'].set_visible(True)
    axmatrix.spines['right'].set_visible(True)
    axmatrix.spines['bottom'].set_visible(True)
    axmatrix.spines['left'].set_visible(True)

    # 5. Lab
    df_lab = df.loc[list(df.index[df.index.str.contains('Lab_')]), :]  # 10 rows
    axmatrix = fig.add_axes([0.2, 0.502, 0.65, 0.06])  # (left, bottom, width, height)
    sns.heatmap(df_lab, vmin=0, xticklabels=False, yticklabels=True, ax=axmatrix,
                cmap=cmap, center=0.5, cbar=False)
    axmatrix.set_yticklabels(axmatrix.get_yticklabels(), rotation=0)
    axmatrix.set_xticklabels(axmatrix.get_xticklabels(), rotation=90)
    axmatrix.yaxis.set_label_position("right")
    axmatrix.yaxis.tick_right()
    axmatrix.set_xticklabels([])
    axmatrix.set_xlabel('')
    axmatrix.spines['top'].set_visible(True)
    axmatrix.spines['right'].set_visible(True)
    axmatrix.spines['bottom'].set_visible(True)
    axmatrix.spines['left'].set_visible(True)

    # 6. Histology
    df_histology = df.loc[list(df.index[df.index.str.contains('Hist_')]), :]  # 26 rows
    axmatrix = fig.add_axes([0.2, 0.342, 0.65, 0.156])  # (left, bottom, width, height)
    sns.heatmap(df_histology, vmin=0, xticklabels=False, yticklabels=True, ax=axmatrix,
                cmap=cmap, center=0.5, cbar=False)
    axmatrix.set_yticklabels(axmatrix.get_yticklabels(), rotation=0)
    axmatrix.set_xticklabels(axmatrix.get_xticklabels(), rotation=90)
    axmatrix.yaxis.set_label_position("right")
    axmatrix.yaxis.tick_right()
    axmatrix.set_xticklabels([])
    axmatrix.set_xlabel('')
    axmatrix.spines['top'].set_visible(True)
    axmatrix.spines['right'].set_visible(True)
    axmatrix.spines['bottom'].set_visible(True)
    axmatrix.spines['left'].set_visible(True)

    # 7. Comorbidity
    df_comorbidity = df.loc[list(df.index[df.index.str.contains('Com_')]), :]  # 14 rows
    axmatrix = fig.add_axes([0.2, 0.254, 0.65, 0.084])  # (left, bottom, width, height)
    sns.heatmap(df_comorbidity, vmin=0, xticklabels=False, yticklabels=True, ax=axmatrix,
                cmap=cmap, center=0.5, cbar=False)
    axmatrix.set_yticklabels(axmatrix.get_yticklabels(), rotation=0)
    axmatrix.set_xticklabels(axmatrix.get_xticklabels(), rotation=90)
    axmatrix.yaxis.set_label_position("right")
    axmatrix.yaxis.tick_right()
    axmatrix.set_xticklabels([])
    axmatrix.set_xlabel('')
    axmatrix.spines['top'].set_visible(True)
    axmatrix.spines['right'].set_visible(True)
    axmatrix.spines['bottom'].set_visible(True)
    axmatrix.spines['left'].set_visible(True)

    # 8. History
    df_history = df.loc[list(df.index[df.index.str.contains('Hty_')]), :]  # 12 rows
    axmatrix = fig.add_axes([0.2, 0.178, 0.65, 0.072])  # (left, bottom, width, height)
    sns.heatmap(df_history, vmin=0, xticklabels=False, yticklabels=True, ax=axmatrix,
                cmap=cmap, center=0.5, cbar=False)
    axmatrix.set_yticklabels(axmatrix.get_yticklabels(), rotation=0)
    axmatrix.set_xticklabels(axmatrix.get_xticklabels(), rotation=90)
    axmatrix.yaxis.set_label_position("right")
    axmatrix.yaxis.tick_right()
    axmatrix.set_xticklabels([])
    axmatrix.set_xlabel('')
    axmatrix.spines['top'].set_visible(True)
    axmatrix.spines['right'].set_visible(True)
    axmatrix.spines['bottom'].set_visible(True)
    axmatrix.spines['left'].set_visible(True)

    # Add colorbar
    axcolor = fig.add_axes(
        [0.2, 0.15, 0.65, 0.02])  # (left, bottom, width, height) fig.add_axes([0.9, 0.1, 0.03, 0.72])
    cbar = plt.colorbar(mappable=axmatrix.collections[0], cax=axcolor, ticks=[0, 0.5, 1], orientation="horizontal")
    # cbar.ax.get_xaxis().labelpad = 15
    cbar.ax.set_xlabel('scaled values', rotation=0, fontsize=12)

    # Add heatmap showing the colors of the clusters
    ax_colors = fig.add_axes([0.2, 0.75, 0.65, 0.02])  # (left, bottom, width, height)
    cmap_colors = LinearSegmentedColormap.from_list('', df_colors.iloc[0, :].values, len(df_colors.iloc[0, :].values))
    gradient = np.linspace(0, 1, len(df_colors.columns))
    gradient = np.vstack((gradient, gradient))
    ax_colors.imshow(gradient, aspect='auto', cmap=cmap_colors)
    ax_colors.get_xticklabels([])
    ax_colors.get_yticklabels([])
    ax_colors.get_xaxis().set_visible(False)
    ax_colors.get_yaxis().set_visible(False)

    plt.savefig(os.path.join(save_folder, "separated_Heatmap.pdf"), bbox_inches='tight')
    plt.close()


def main_all_patients(save_folder):
    adata = sc.read('/Users/christina.hillig/PycharmProjects/Eyerich_ERCGrant/Molecular_subtypes/input/data_freeze/Lesion_RNAseq_20210720_patient_meta_data_v04__CH_20220209_TherapyResponse_v04_v220715.h5')
    adata, sorter = add_colors.diag_order_lesion(adata=adata)
    adata, _ = add_colors.sdiag_order_lesion(adata=adata)
    adata, _ = add_colors.pattern_order_lesion(adata=adata)

    df_diag = pd.DataFrame(adata.obs['diag'].copy(), columns=['diag'])

    # Attribute preparation
    df = prepare_attributes(adata=adata, save_folder=save_folder)

    print('Number of clinical attributes: ', len(df.columns))

    # Sort df by diagnosis
    # Create the dictionary that defines the order for sorting
    sorter_index = dict(zip(sorter, range(len(sorter))))
    df_diag['Rank'] = df_diag['diag'].map(sorter_index)
    df_diag['Rank'] = df_diag['Rank'].astype(int)
    df_diag.sort_values(by='Rank', inplace=True)
    df = df.loc[list(df_diag.index)]

    # Plot Heatmap
    # Prepare a vector of color mapped to the 'diag' column
    my_palette = dict(zip(list(adata.obs.diag.cat.categories), list(adata.uns['diag_colors'])))
    row_colors = adata.obs.loc[list(df_diag.index), 'diag'].map(my_palette)

    df = df.transpose()
    # remove NaN rows which are NaN in all columns
    df = df.dropna(how='all')  # Com_uv
    # sort attributes
    sns.set(font_scale=1.2)
    sns.clustermap(df, metric="euclidean", method="ward", col_colors=row_colors, row_cluster=False, col_cluster=False,
                   cmap="RdYlBu", figsize=(12, 16), dendrogram_ratio=.2, annot_kws={"size": 8}, xticklabels=False,
                   yticklabels=True, cbar_kws={"shrink": 2}, cbar_pos=(0.1, .2, .03, .4))
    plt.xlabel('')
    plt.yticks(fontsize=18)
    # plt.tight_layout()
    plt.savefig(os.path.join(save_folder, 'Heatmap.pdf'))
    plt.close('all')


def main(save_folder):
    adata = sc.read(os.path.join(
        '/Volumes', 'CH__data', 'Projects', 'Eyerich_AG_projects', 'BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer',
        'analysis', 'Molecular_subtypes', 'input', 'h5_files', 'LESION',
        'Lesion_RNAseq_20210720_patient_meta_data_v04__CH_20220209_TherapyResponse_v04_v220715_TherapyResponseCorrected__Endotypes_230620.h5'))

    df_corrected_metadata = pd.read_excel(os.path.join(
        '/Volumes', 'CH__data', 'Projects', 'Eyerich_AG_projects', 'BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer',
        'raw_data', 'clinical_data', 'patient_meta_data_final.xlsx'))
    # rename Sex.x, Helmholtz.identifyer', Diag_diag.x
    df_corrected_metadata = df_corrected_metadata.rename(
        columns={'Helmholtz.identifyer': 'Helmholtz_identifyer', 'Diag_diag.x': 'Diag_diag_x', 'Sex.x': 'Sex_x'})
    df_corrected_metadata.index = df_corrected_metadata['Helmholtz_identifyer']

    # common columns
    # set(df_corrected_metadata.columns) - set(adata.obs.columns)
    common_cols = np.intersect1d(adata.obs.columns, df_corrected_metadata.columns)
    assert np.all(df_corrected_metadata.loc[adata.obs.index, common_cols].index == adata.obs.index)
    # TODO might have to adjust some variables -> Peter
    adata.obs.loc[:, common_cols] = df_corrected_metadata.loc[adata.obs.index, common_cols]

    adata, sorter = add_colors.diag_order_lesion(adata=adata)
    adata, _ = add_colors.sdiag_order_lesion(adata=adata)
    adata, _ = add_colors.pattern_order_lesion(adata=adata)
    df_diag = pd.DataFrame(adata.obs['diag'].copy(), columns=['diag'])

    # Attribute preparation
    df = prepare_attributes(adata=adata, save_folder=save_folder)

    print('Number of clinical attributes: ', len(df.columns))  # 89

    # Sort df by diagnosis
    # Create the dictionary that defines the order for sorting
    sorter_index = dict(zip(sorter, range(len(sorter))))
    df_diag['Rank'] = df_diag['diag'].map(sorter_index)
    df_diag['Rank'] = df_diag['Rank'].astype(int)
    df_diag.sort_values(by='Rank', inplace=True)
    df = df.loc[list(df_diag.index)]

    # rename and sort attributes
    attr = ['Com_arths', 'Com_ah', 'Com_iaa', 'Com_rca', 'Com_dm',
            'Com_ast', 'Com_hepa', 'Com_renal', 'Com_tonsil',
            'Com_ibd', 'Com_nail', 'Com_muc', 'Com_scalp',
            'Hty_fh', 'Hty_exanth', 'Hty_self', 'Hty_photos',
            'Hty_progr', 'Hty_smoker', 'Lab_csa', 'Lab_specificIgE',
            'Sex', 'Therapeutic target_IL-17', 'Therapeutic target_IL-23',
            'Therapeutic target_MTX', 'Therapeutic target_TNF', 'Diag_PASI',
            'Diag_SCORAD', 'Diag_PGA', 'Diag_DLQI', 'Morph_erythema', 'Morph_Elev',
            'Morph_papu', 'Morph_vesicles', 'Morph_pustules', 'Morph_scaling',
            'Morph_dryness', 'Hty_sport', 'Hty_alc', 'Hist_SE', 'Hist_Para_quali',
            'Hist_Aka_quali', 'Hist_Granu', 'Hist_Serum', 'Hist_Dist_Lympho',
            'Hist_Microabscess', 'Hist_Bacteria', 'Hist_Neutro_quali',
            'Hist_Keratoses', 'Hist_Nr_Keratoses', 'Hist_ID_distribution',
            'Hist_ID_quali', 'The_sor', 'The_dor', 'Hty_onset', 'Hty_pruritus',
            'Com_BMI', 'Hist_Aka_quant', 'Hist_Cap', 'Hist_Eosino',
            'Hist_Exocytosis', 'Hist_Hyper_quant', 'Hist_ID_quant', 'Hist_LCV',
            'Hist_Lympho', 'Hist_Mucin', 'Hist_Neutro_quant', 'Hist_Ortho_quant',
            'Hist_Para_quant', 'Hist_Spongiosis', 'Hty_chronic', 'Lab_GPT',
            'Lab_Hb', 'Lab_crea', 'Lab_eosino', 'Lab_granulo', 'Lab_leuco',
            'Lab_lympho', 'Lab_totalIgE', 'The_GPT', 'The_Hb', 'The_crea',
            'The_eosino', 'The_granulo', 'The_leuco', 'The_lympho', 'The_weeks', 'age']

    df.columns = attr
    attr_sorted = [
        'age', 'Sex',
        'Morph_erythema', 'Morph_Elev', 'Morph_papu', 'Morph_vesicles',
        'Morph_pustules', 'Morph_scaling', 'Morph_dryness',
        'Diag_PASI', 'Diag_SCORAD', 'Diag_PGA', 'Diag_DLQI',
        'Therapeutic target_IL-17', 'Therapeutic target_IL-23', 'Therapeutic target_TNF', 'Therapeutic target_MTX',
        'The_GPT', 'The_Hb', 'The_crea', 'The_eosino', 'The_granulo', 'The_leuco', 'The_lympho', 'The_weeks',
        'The_dor', 'The_sor',
        'Lab_GPT','Lab_Hb', 'Lab_crea', 'Lab_eosino', 'Lab_granulo', 'Lab_leuco', 'Lab_lympho', 'Lab_totalIgE',
        'Lab_csa', 'Lab_specificIgE',
        'Hist_SE', 'Hist_Para_quali', 'Hist_Aka_quali', 'Hist_Granu', 'Hist_Serum', 'Hist_Dist_Lympho',
        'Hist_Microabscess', 'Hist_Bacteria', 'Hist_Neutro_quali', 'Hist_Keratoses', 'Hist_Nr_Keratoses',
        'Hist_ID_distribution', 'Hist_ID_quali', 'Hist_Aka_quant', 'Hist_Cap', 'Hist_Eosino', 'Hist_Exocytosis',
        'Hist_Hyper_quant', 'Hist_ID_quant', 'Hist_LCV', 'Hist_Lympho', 'Hist_Mucin', 'Hist_Neutro_quant',
        'Hist_Ortho_quant', 'Hist_Para_quant', 'Hist_Spongiosis',
        'Com_arths', 'Com_ah', 'Com_iaa', 'Com_rca', 'Com_dm', 'Com_ast', 'Com_hepa', 'Com_renal', 'Com_tonsil',
        'Com_ibd', 'Com_nail', 'Com_muc', 'Com_scalp', 'Com_BMI',
        'Hty_fh', 'Hty_exanth', 'Hty_self', 'Hty_photos', 'Hty_progr', 'Hty_smoker',  'Hty_alc', 'Hty_chronic',
        'Hty_onset', 'Hty_pruritus', 'Hty_sport']
    df = df[attr_sorted]

    df = df.transpose()  # shape: 88 x 342
    # remove NaN rows which are NaN in all columns
    df = df.dropna(how='all')  # Com_uv

    print('Number of final clinical attributes: ', len(df.index))
    df.transpose().to_excel(os.path.join(save_folder, 'Clinical_attributes.xlsx'))

    # Plot clustermap of clinical attributes
    # Prepare a vector of color mapped to the 'cyl' column
    my_palette = dict(zip(list(adata.obs.diag.cat.categories), list(adata.uns['diag_colors'])))
    col_colors = adata.obs.loc[list(df_diag.index), 'diag'].map(my_palette)
    sns.set(font_scale=1.2)

    kws = dict(cbar_kws=dict(ticks=[0, 0.50, 1], orientation='horizontal', shrink=2),
               annot_kws={"size": 6}, figsize=(6, 8))
    # Greys, 'RdYlBu'
    g = sns.clustermap(
        df, metric="euclidean", method="ward", col_colors=col_colors, row_cluster=False, col_cluster=False,
        cmap="RdYlBu", dendrogram_ratio=.01, xticklabels=False, yticklabels=True, cbar_pos=(0.042, .03, .566, .02),
        linewidths=0, **kws)
    ax = g.ax_heatmap
    ax.set_xlabel('')
    ax.yaxis.set_tick_params(labelsize=8)
    g.savefig(os.path.join(save_folder, 'Heatmap.pdf'), bbox_inches='tight')
    plt.close('all')

    df_colors = pd.DataFrame(col_colors).T
    separated_heatmap(df=df, df_colors=df_colors, save_folder=save_folder, figsize=(6, 12), cmap='RdYlBu')


if __name__ == '__main__':
    output_dir = os.path.join('/Volumes', 'CH__data', 'Projects', 'Eyerich_AG_projects',
                              'BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer', 'analysis', 'Molecular_subtypes',
                              'output', 'Figure_1B_Heatmap_corrected', str(date.today()))
    os.makedirs(output_dir, exist_ok=True)

    main(save_folder=output_dir)
