from scripts.feature_engineering import imputation, scaling, encoding
from scripts.utils import add_colors
from sklearn import preprocessing

from datetime import date

import matplotlib.pyplot as plt
import matplotlib as mpl
import umap

import seaborn as sns
import scanpy as sc
import pickle

import os
# %matplotlib notebook

import pandas as pd
import numpy as np
import itertools
import collections
from math import pi


def encode_impute_data(bulk_data, meta_data_legend, save_folder):
    diagnosis_labels = bulk_data.obs['diag']

    dichotome_cols = [meta_data_legend['Abbreviation'][for_ind]
                      for for_ind, for_x in enumerate(meta_data_legend['Explanation ']) if 'dichotome' in for_x]
    dichotome_cols.extend(['Sex_x', 'Therapeutic target'])
    nominal_cols = [e for e in dichotome_cols if e in bulk_data.obs]  # 23
    # + Lab_totalIgE
    for for_col in nominal_cols:
        bulk_data.obs[for_col] = bulk_data.obs[for_col].astype('category')
    # 2.2.1.2 Apply encoding
    df_encoded_nominals = encoding.encode_nominal_feature(
        data=bulk_data.obs[nominal_cols], dummy_variable_encoding=True)
    # 2.2.1.3 Handle missing values categories in encoded features
    df_encoded_nominals = encoding.imputation_replace_zero_nan(
        df_data=df_encoded_nominals, original_nominal_feature_names=nominal_cols)  # 44
    # set index
    df_encoded_nominals.index = bulk_data.obs.index

    # 2.2.2 Ordinal Encoding
    ordinal_cols_1 = [meta_data_legend['Abbreviation'][for_ind]
                      for for_ind, for_x in enumerate(meta_data_legend['Explanation ']) if 'Morphology' in for_x]  # 7
    ordinal_cols_2 = [meta_data_legend['Abbreviation'][for_ind]
                      for for_ind, for_x in enumerate(meta_data_legend['Explanation ']) if
                      '(qualitative)' in for_x]  # 14
    ordinal_cols = list(itertools.chain(ordinal_cols_1, ordinal_cols_2))  # 21
    # remove columns which are not in bulk.obs
    ordinal_cols = [e for e in ordinal_cols if e in bulk_data.obs]  # 19
    # some have to be added manually .. such as The_ except The_weeks
    ordinal_cols.extend(['Hist_Nr_Keratoses', 'Hty_sport', 'Hty_alc', 'Hty_onset', 'Hty_pruritus', 'Diag_PGA',
                         'The_dor', 'The_sor', 'The_leuco', 'The_Hb', 'The_lympho',
                         'The_granulo', 'The_eosino', 'The_crea', 'The_GPT'])  # 24
    # Apply ordinal encoding on the ordinal features
    df_encoded_ordinals = encoding.encode_ordinal_categories(data=bulk_data.obs[ordinal_cols])
    # set index
    df_encoded_ordinals.index = bulk_data.obs.index

    # 2.2.3 Treat continuous variable -> Decision 16.11.2021: leave them as they are
    continuous_cols_1 = [meta_data_legend['Abbreviation'][for_ind]
                         for for_ind, for_x in enumerate(meta_data_legend['Explanation ']) if
                         '(quantitative)' in for_x]  # 13
    continuous_cols_2 = [meta_data_legend.loc[for_ind, 'Abbreviation']
                         for for_ind, for_x in enumerate(meta_data_legend['Explanation '])
                         if '(number)' in for_x and not pd.isnull(meta_data_legend.loc[for_ind, 'Unit'])]  # 9
    # some have to be added manually ..
    continuous_cols = list(
        itertools.chain(continuous_cols_1, continuous_cols_2, ['Com_BMI', 'Hty_chronic', 'age', 'The_weeks']))  # 25
    # remove columns which are not in bulk.obs
    continuous_cols = [e for e in continuous_cols if e in bulk_data.obs]  # 24
    # remove categories from these columns
    df_continous = bulk_data.obs[continuous_cols].copy()
    df_continous = df_continous.astype('float')

    # 2.2.5 Combine features again
    # Check if I capture all columns
    check_cols = list(itertools.chain(nominal_cols, ordinal_cols, continuous_cols))  # 70
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
            'raw_data', 'clinical_data', "20210720_patient_meta_data_v04__CH.xlsx"), sheet_name='legend')

    # Encode data
    df_encoded, nominal_cols, ordinal_cols, continuous_cols = encode_impute_data(
        bulk_data=adata, meta_data_legend=meta_data_legend, save_folder=save_folder)

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


def plot_umap_clinical_attributes(scaled_data, save_folder, labels, data_type='Pattern', mapping=None):
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
    # print(embedding.shape)

    try:
        cmap = plt.cm.Paired  # define the colormap
        cmap_bounds = np.arange(-0.5, max(labels) + 1.5)
        cmap_norm = mpl.colors.BoundaryNorm(cmap_bounds, cmap.N)
    except ValueError:
        cmap = plt.cm.tab20  # define the colormap
        cmap_bounds = np.arange(-0.5, max(labels) + 1.5)
        cmap_norm = mpl.colors.BoundaryNorm(cmap_bounds, cmap.N)

    fig, ax = plt.subplots(ncols=2, constrained_layout=True, gridspec_kw={'width_ratios': [10, 1]})
    ax[0].scatter(embedding[:, 0], embedding[:, 1], c=labels, cmap=cmap, norm=cmap_norm)
    ax[0].set_xlabel('UMAP1')
    ax[0].set_ylabel('UMAP2')
    # create a second axes for the colorbar
    # ax2 = fig.add_axes([0.7, 0.1, 0.03, 0.8])
    cb = mpl.colorbar.ColorbarBase(ax[1], cmap=cmap, norm=cmap_norm, spacing='proportional',
                                   ticks=np.linspace(0, max(labels), max(labels)+1),
                                   boundaries=cmap_bounds, format='%1i')
    if mapping.__class__ == 'NoneType':
        cb.set_ticklabels(['1', '2a', '2b', '3', '4a', '4b', '5'])
    else:
        od = collections.OrderedDict(sorted(mapping.items()))
        cb.set_ticklabels(list(od.values()))
    # remove borders from plot
    ax[0].spines['top'].set_visible(False)
    ax[0].spines['right'].set_visible(False)
    ax[0].spines['bottom'].set_visible(True)
    ax[0].spines['left'].set_visible(True)

    # plt.gca().set_aspect('equal', 'datalim')
    # plt.title('UMAP projection of clinical attributes')
    # plt.tight_layout()
    fig.savefig(os.path.join(save_folder, 'UMAP_embedding_clinical_attributes_{}.pdf'.format(data_type)))
    plt.close(fig=fig)

    return embedding


def load_molecular_subtypes(adata, data_root):
    try:
        with open(os.path.join(data_root, 'Optimalres_0.5__OptimalGPTK_7.pkl'), 'rb') as ff:
            df_res05 = pickle.load(ff)
        subtypes_res05 = df_res05.iloc[0]['Molecular subtypes']
        with open(os.path.join(data_root, 'Optimalres_0.8__OptimalGPTK_7.pkl'), 'rb') as ff:
            df_res08 = pickle.load(ff)
        subtypes_res08 = df_res08.iloc[0]['Molecular subtypes']
        with open(os.path.join(data_root, 'Optimalres_0.9__OptimalGPTK_7.pkl'), 'rb') as ff:
            df_res09 = pickle.load(ff)
        subtypes_res09 = df_res09.iloc[0]['Molecular subtypes']
    except ModuleNotFoundError:
        df_res05 = pd.read_pickle(os.path.join(data_root, 'Optimalres_0.5__OptimalGPTK_7.pkl'))
        subtypes_res05 = df_res05.iloc[0]['Molecular subtypes']
        df_res08 = pd.read_pickle(os.path.join(data_root, 'Optimalres_0.8__OptimalGPTK_7.pkl'))
        subtypes_res08 = df_res08.iloc[0]['Molecular subtypes']
        df_res09 = pd.read_pickle(os.path.join(data_root, 'Optimalres_0.9__OptimalGPTK_7.pkl'))
        subtypes_res09 = df_res09.iloc[0]['Molecular subtypes']

    assert np.all(df_res05.iloc[0]['MUC IDs'] == adata.obs.index), 'Index are not in the same order as in adata'

    adata.obs['Molecular Subtype res0.5'] = subtypes_res05
    adata.obs['Molecular Subtype res0.5'] = adata.obs['Molecular Subtype res0.5'].astype('category')
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


def plot_umap(df_umap, palette, save_folder, hue):
    fig, ax = plt.subplots(figsize=(6, 6))
    sns.scatterplot(data=df_umap, x='x', y='y', hue=hue, palette=palette, edgecolor='k', linewidth=0, s=100)
    ax.set_ylabel('UMAP2', fontsize=18)
    ax.set_xlabel('UMAP1', fontsize=18)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.get_legend().remove()
    # Put a legend to the right of the current axis
    # ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=2)
    plt.tight_layout()
    plt.savefig(os.path.join(save_folder, 'Figure_1C_UMAP_{}.pdf'.format(hue)))
    plt.close('all')


def plot_lollipop_plot(df_patients, save_folder):
    my_range = range(0, len(df_patients.transpose().index))
    for ind, patient in enumerate(df_patients.index):
        markerline, stemlines, baseline = plt.stem(
            df_patients.loc[patient, :], linefmt='grey', orientation='horizontal')
        markerline.set_markerfacecolor('black')
        markerline.set_markeredgecolor('black')
        plt.xlabel('Scaled values', fontsize=18)
        plt.ylabel('Clinical attributes', fontsize=18)
        plt.yticks(my_range, df_patients.transpose().index, fontsize=16)
        plt.xticks(fontsize=16)
        plt.gca().invert_yaxis()
        sns.despine()
        plt.tight_layout()
        plt.savefig(os.path.join(save_folder, 'Lollipop_patient_{}.pdf'.format(patient)))
        plt.close('all')


def plot_areachart_plot(df_patients, save_folder):
    # For vertical area chart: https://stackoverflow.com/questions/50802556/how-to-plot-a-vertical-area-plot-with-pandas
    for ind, patient in enumerate(df_patients.index):
        df_patients.iloc[ind, :].transpose().plot.area(color='lightgrey')
        plt.xticks(rotation=315, ha='left', fontsize=18)
        plt.yticks(fontsize=18)
        sns.despine()
        plt.tight_layout()
        plt.savefig(os.path.join(save_folder, 'AreaChart_patient_{}.pdf'.format(patient)))
        plt.close('all')

    df_patients.transpose().plot.area(stacked=False)
    plt.xticks(rotation=315, ha='left', fontsize=18)
    plt.yticks(fontsize=18)
    sns.despine()
    plt.tight_layout()
    plt.savefig(os.path.join(save_folder, 'AreaChart_patient__{}.pdf'.format("_".join(df_patients.index))))
    plt.close('all')


def plot_radar_plot(df_patients, categories, save_folder):
    N = len(categories)
    # We are going to plot the first line of the data frame.
    # But we need to repeat the first value to close the circular graph:
    letter = ['E', 'D', 'F']
    for ind, patient in enumerate(df_patients['Helmholtz_identifyer']):
        # values
        values = df_patients.loc[ind].drop('Helmholtz_identifyer').values.flatten().tolist()
        values += values[:1]

        # What will be the angle of each axis in the plot? (we divide the plot / number of variable)
        angles = [n / float(N) * 2 * pi for n in range(N)]
        angles += angles[:1]

        # Initialise the spider plot
        ax = plt.subplot(111, polar=True)
        # Draw one axe per variable + add labels
        plt.xticks(angles[:-1], categories, color='k', size=18)
        # Draw ylabels
        ax.set_rlabel_position(0)
        plt.yticks([0.2, 0.4, 0.6, 0.8], ["0.2", "0.4", "0.6", "0.8"], color="grey", size=12)
        ax.set_ylim(0, 1)
        # Plot data
        ax.plot(angles, values, linewidth=1, linestyle='solid', color='grey')
        # Fill area
        ax.fill(angles, values, 'grey', alpha=0.4)

        plt.savefig(os.path.join(
            save_folder, 'Figure_1{}_Radarplot_patient_{}.pdf'.format(letter[ind], patient)),
            bbox_inches='tight')
        plt.close('all')


def main(save_folder):
    adata = sc.read(os.path.join(
        '/Volumes', 'CH__data', 'Projects', 'Eyerich_AG_projects', 'BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer',
        'analysis', 'Molecular_subtypes', 'input', 'h5_files', 'LESION',
        'Lesion_RNAseq_20210720_patient_meta_data_v04__CH_20220209_TherapyResponse_v04_v220715_TherapyResponseCorrected__Endotypes_230620.h5'))
    # # New: add molecular subthypes
    # adata = load_molecular_subtypes(adata=adata, data_root=data_root)
    # adata, _ = add_colors.diag_order_lesion(adata=adata)
    # adata, _ = add_colors.sdiag_order_lesion(adata=adata)
    # adata, _ = add_colors.pattern_order_lesion(adata=adata)
    # adata, _ = add_colors.endotype_order_lesion(adata=adata, obs_name='Molecular Subtype res0.9')

    labels = pd.DataFrame(columns=['Pattern'], data=adata.obs['Pattern'])
    rename_patterns = {"Pattern": {"1": 1, "2a": 2, "2b": 3, "3": 4, "4a": 5, "4b": 6, "5": 7, "UD": 8}}
    labels = labels.replace(rename_patterns)
    labels['Pattern'] = labels['Pattern'].astype('category')
    diagnosis_labels = adata.obs['diag']

    # 2.2.8 Label Encoding
    label_train_encoded = encoding.encode_label(label=labels['Pattern'].values)
    labels = label_train_encoded.transform(labels.values)

    # Attribute preparation
    df = prepare_attributes(adata=adata, save_folder=save_folder)

    # convert columns to categories
    ticks_pattern = dict({0: "1", 1: "2a", 2: "2b", 3: "3", 4: "4a", 5: "4b", 6: "5", 7: "UD"})
    # get embeddings
    # TODO adjust color
    umap_emb = plot_umap_clinical_attributes(
        scaled_data=df, save_folder=save_folder, labels=labels, data_type='Pattern_clinical_attributes',
        mapping=ticks_pattern)

    df_umap = pd.DataFrame.from_dict({
        'x': umap_emb[:, 0], 'y': umap_emb[:, 1], 'diag': diagnosis_labels, 'Pattern': labels,
        'Endotypes': adata.obs['Endotypes']})

    df_umap[['x', 'y', 'diag']].to_excel(os.path.join(save_folder, 'Figure_1C_Clinical_attributes_UMAP.xlsx'))

    plot_umap(df_umap=df_umap, palette=list(adata.uns['diag_colors']),
              save_folder=save_folder, hue='diag')

    plot_umap(df_umap=df_umap, palette=list(adata.uns['Pattern_colors']),
              save_folder=save_folder, hue='Pattern')

    plot_umap(df_umap=df_umap, palette=list(adata.uns['Endotypes_colors']),
              save_folder=save_folder, hue='Endotypes')

    # df columns = 'patient', 'clinical attribute', 'value'
    # similar: 'MUC20721', 'MUC7393', different: 'MUC4241'
    pso_patients = ['MUC20721', 'MUC7393', 'MUC4241']
    # pso_patients = adata.obs[adata.obs['diag'] == 'psoriasis'].index[:3]
    df_patients = pd.DataFrame(data=df.loc[pso_patients, :])
    df_patients = df_patients[['age', 'Sex_x_m', 'Morph_pustules', 'Diag_PGA', 'Lab_lympho',
                               'Hist_Aka_quant', 'Hty_sport', 'Com_iaa_2.0']]
    # df_patients = pd.melt(df_patients, id_vars=['Pseudo ID'])
    # # use 5 clinical attributes which describe the most variability (check PCA of clinical attributes?)
    # df_patients.columns = ['patient', 'clinical attribute', 'value']

    df_patients.to_excel(os.path.join(save_folder, 'Figure_1DEF_Clinical_attributes_AreaChart.xlsx'))

    # Area Chart
    plot_areachart_plot(df_patients, save_folder)

    # Make Lollipop plot
    plot_lollipop_plot(df_patients=df_patients, save_folder=save_folder)

    # Plot Radarplot
    df_patients = df_patients.reset_index()
    df_patients.to_excel(os.path.join(save_folder, 'Figure_1DEF_Clinical_attributes_Radarplot.xlsx'))

    # number of variable
    # categories = df_patients.columns[1:]
    # rename categories
    categories = ['age', 'Sex', 'Morphology', 'Severity', 'Labor', 'Histology', 'History', 'Comorbidity']
    plot_radar_plot(df_patients=df_patients, categories=categories, save_folder=save_folder)


if __name__ == '__main__':
    output_dir = os.path.join('/Volumes', 'CH__data', 'Projects', 'Eyerich_AG_projects',
                              'BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer', 'analysis', 'Molecular_subtypes',
                              'output', 'Figure_1C', str(date.today()))
    os.makedirs(output_dir, exist_ok=True)

    main(save_folder=output_dir)
