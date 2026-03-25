"""
@description: Find enriched clinical attributes in Endotypes,
distinguish between categorical nominal, categorical ordinal and Continuous attributes

@author : Christina Hillig
"""

from scripts.feature_engineering import imputation, encoding
from scripts.utils import check_distributions, perform_test_nonparametric

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

from sklearn import preprocessing
import seaborn as sns
import scanpy as sc
import pandas as pd
import numpy as np
import itertools
import collections

import scipy.cluster.hierarchy as sch
from scipy.spatial.distance import pdist

import os
from datetime import date

import statsmodels.stats.multicomp as mc
from scikit_posthocs import posthoc_dunn
from scipy.stats import f_oneway, chi2_contingency, fisher_exact
from itertools import combinations
from statsmodels.sandbox.stats.multicomp import multipletests

import rpy2.robjects.packages as rpackages
from rpy2.robjects import numpy2ri, pandas2ri

numpy2ri.activate()
pandas2ri.activate()

base = rpackages.importr("base")
dollar = base.__dict__["$"]
clinfun = rpackages.importr("clinfun")


def diag_order_lesion(adata):
    diag_order = [
        'lichen planus', 'lupus erythematosus', 'eczema', 'prurigo simplex subacuta', 'bullous pemphigoid', 'psoriasis',
        'pityriasis rubra pilaris', 'morphea', 'venous ulcer', 'systemic sclerosis', 'granuloma annulare',
        'sarcoidosis', 'psoriasis pustulosa', 'pyoderma gangrenosum', 'cutaneous lymphoma',
        'cutaneous side effects of biologics', 'darier disease', 'keratosis lichenoides chronica', 'erythrodermia',
        'parapsoriasis', 'undefined']
    adata.obs['diag'] = adata.obs['diag'].cat.reorder_categories(diag_order)
    adata.uns['diag_colors'] = [
        'orange', 'goldenrod', 'maroon', 'firebrick', 'salmon', 'blue', 'dodgerblue', 'forestgreen', 'darkolivegreen',
        'seagreen', 'springgreen', 'palegreen', 'darkviolet', 'mediumorchid', 'dimgrey', 'darkgray', 'silver', 'grey',
        'slategrey', 'lightgrey', 'whitesmoke']

    return adata, diag_order


def pattern_order_lesion(adata):
    pattern_order = ['1', '2a', '2b', '3', '4a', '4b', '5', 'undefined']
    adata.obs['Pattern'] = adata.obs['Pattern'].cat.reorder_categories(pattern_order)
    adata.uns['Pattern_colors'] = ['darkorange', 'darkred', 'tomato', 'royalblue', 'darkgreen',
                                   'green', 'darkviolet', 'grey']

    return adata, pattern_order


def encode_impute_data(bulk_data, meta_data_legend, save_folder):
    diagnosis_labels = bulk_data.obs['diag']

    dichotome_cols = [meta_data_legend['Abbreviation'][for_ind]
                      for for_ind, for_x in enumerate(meta_data_legend['Explanation ']) if 'dichotome' in for_x]
    dichotome_cols.extend(['Sex_x', 'Therapeutic target'])
    nominal_cols = [e for e in dichotome_cols if e in bulk_data.obs]  # 24
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
    print('Number of nominal categories: ', len(nominal_cols))  # 24
    print('Number of encoded nominal categories: ', df_encoded_nominals.shape[1])  # 27

    # 2.2.2 Ordinal Encoding
    ordinal_cols_1 = [meta_data_legend['Abbreviation'][for_ind]
                      for for_ind, for_x in enumerate(meta_data_legend['Explanation ']) if 'Morphology' in for_x]  # 7
    ordinal_cols_2 = [meta_data_legend['Abbreviation'][for_ind]
                      for for_ind, for_x in enumerate(meta_data_legend['Explanation ']) if
                      '(qualitative)' in for_x]  # 14
    ordinal_cols = list(itertools.chain(ordinal_cols_1, ordinal_cols_2))  # 21
    # remove columns which are not in bulk.obs: 'Hty_negative', 'Hty_positive'
    ordinal_cols = [e for e in ordinal_cols if e in bulk_data.obs]  # 19
    # some have to be added manually .. such as The_ except The_weeks
    ordinal_cols.extend(['Hist_Nr_Keratoses', 'Hty_sport', 'Hty_alc', 'Hty_onset', 'Hty_pruritus',
                         'The_dor', 'The_sor', 'The_leuco', 'The_Hb', 'The_lympho',
                         'The_granulo', 'The_eosino', 'The_crea', 'The_GPT'])  # 33
    # Apply ordinal encoding on the ordinal features
    df_encoded_ordinals = encoding.encode_ordinal_categories(data=bulk_data.obs[ordinal_cols])
    # set index
    df_encoded_ordinals.index = bulk_data.obs.index
    print('Number of nominal categories: ', len(ordinal_cols))  # 33
    print('Number of encoded ordinal categories: ', df_encoded_ordinals.shape[1])  # 33

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
    continuous_cols = [e for e in continuous_cols if e in bulk_data.obs]  # 25
    # remove categories from these columns
    df_continous = bulk_data.obs[continuous_cols].copy()
    df_continous = df_continous.astype('float')
    print('Number of continuous variables: ', df_continous.shape[1])  # 25

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


def encode_labels(adata):
    labels = pd.DataFrame(columns=['Pattern'], data=adata.obs['Pattern'])
    rename_patterns = {"Pattern": {"1": 1, "2a": 2, "2b": 3, "3": 4, "4a": 5, "4b": 6, "5": 7, "undefined": 8}}
    labels = labels.replace(rename_patterns)
    labels['Pattern'] = labels['Pattern'].astype('category')

    # 2.2.8 Label Encoding
    label_train_encoded = encoding.encode_label(label=labels.values)
    encoded_labels = label_train_encoded.transform(labels.values)

    return encoded_labels


def prepare_attributes(adata, input_dir_attr, save_folder):
    meta_data_legend = pd.read_excel(
        os.path.join(input_dir_attr, "20210720_patient_meta_data_v04__CH.xlsx"), sheet_name='legend')

    # Encode data
    df_encoded,  nominal_cols, ordinal_cols, continuous_cols = encode_impute_data(
        bulk_data=adata, meta_data_legend=meta_data_legend, save_folder=save_folder)

    print('Number of variables after filtering: ', df_encoded.shape[1])

    # Drop columns with more than 70% missing values
    df_percentage = df_encoded.isnull().mean() * 100
    variables_to_drop = list(df_percentage[df_percentage > 70].index)
    print('Attributes to drop due to missing values > 70%: ', variables_to_drop)
    df_encoded = df_encoded.drop(variables_to_drop, axis=1)  # 73
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

    return df, nominal_cols, ordinal_cols, continuous_cols


def save_encoded_features(df_scaled, save_folder, nominal_cols, ordinal_cols, continuous_cols):
    # Add row with kind of category for each attribute
    df_encoded_imp__scaled_features = df_scaled.copy()
    number_rows = len(df_encoded_imp__scaled_features)
    df_encoded_imp__scaled_features.loc[number_rows] = list(
        np.repeat('Dichotomous', df_encoded_imp__scaled_features[nominal_cols].shape[1])) + list(
        np.repeat('Ordinal', df_encoded_imp__scaled_features[ordinal_cols].shape[1])) + list(
        np.repeat('Continuous', df_encoded_imp__scaled_features[continuous_cols].shape[1]))
    df_encoded_imp__scaled_features = df_encoded_imp__scaled_features.rename(index={number_rows: 'Datatype'})

    df_encoded_imp__scaled_features.to_excel(os.path.join(save_folder, 'encoded_imp__scaled_features.xlsx'))


def perform_test_continuous_variables(df, continuous_cols, encoded__labels, save_folder):
    save_folder_tmp = os.path.join(save_folder, 'Continuous')
    os.makedirs(save_folder_tmp, exist_ok=True)
    # 1. Read out continuous attributes and into normal distributed and non-parameteric
    df_continuous = df[continuous_cols]
    df_temp = df_continuous.copy()

    unique_labels = np.unique(encoded__labels)
    names = [(str(s), 'Rest') for s in unique_labels]

    dict_all_continuous_attributes = {
        "Anova": pd.DataFrame(index=list(df_continuous.columns), columns=names),
        "Kruskal-Wallis": pd.DataFrame(index=list(df_continuous.columns), columns=names)}
    df_kruskal = pd.DataFrame()
    df_anova = pd.DataFrame()
    # Create 1 vs Rest groups
    for group in names:
        df_temp['group'] = list(encoded__labels)
        df_temp.loc[df_temp['group'] != group[0], 'group'] = 'Rest'
        df_temp['group'] = df_temp['group'].astype('category')

        # 0. Check which attribute is normally distributed using Shapiro test
        df_var_normal_test = pd.DataFrame(
            columns=['Normal_distribution', 'var_homogeneity_stats', 'var_homogeneity_pval'],
            index=df_continuous.columns)
        for ind, col in enumerate(df_continuous.columns):
            df_var_normal_test.loc[col, 'Normal_distribution'] = check_distributions.check_normal_distribution(
                df=df_continuous[col])
            stats_temp, p_temp = check_distributions.check_variance_homogeneity(
                data=df_temp, feature=col, predictor='group', normal=True)
            df_var_normal_test.loc[col, 'var_homogeneity_stats'] = stats_temp
            df_var_normal_test.loc[col, 'var_homogeneity_pval'] = p_temp

        # 1.2 group per encoded group labels
        data_temp = [df_continuous.loc[ids, :] for ids in df_temp.groupby('group').groups.values()]

        # Compare in a one vs rest fashion
        pval_list = {'Anova': [], 'Kruskal-Wallis': []}
        # Loop through variables
        for ind_var_normal in df_var_normal_test.index:
            variable_classes = [df_ind.loc[:, ind_var_normal].values for df_ind in data_temp]

            # plot class distributions
            res = dict(zip(np.arange(0, len(np.unique(encoded__labels))), variable_classes))

            # Apply first layer Anova or Kruskal Wallis test
            if df_var_normal_test.loc[ind_var_normal, 'Normal_distribution']:
                # Apply Anova
                if df_var_normal_test.loc[ind_var_normal, 'var_homogeneity_pval'] > 0.05:
                    # Perform oneway-Anova
                    _, p_temp = f_oneway(variable_classes)
                    pval_list['Anova'].append((p_temp, ind_var_normal))
                else:
                    print('Test for Normal + unequal Variance: Not yet implemented')
                    # Welch one-factor ANOVA with ad hoc Welch t comparisons
            else:
                # Apply Kruskal-Wallis
                # Perform non-parametric test eg Kruskal-Wallis
                _, p_temp = perform_test_nonparametric.perform_kruskal_wallis_test(
                    data=variable_classes, features=None, predictor=None)
                pval_list['Kruskal-Wallis'].append((p_temp, ind_var_normal))

        # Save Anova and Kruskal-Wallis test results to dataframe
        d_swap = {v: k for k, v in dict(pval_list['Anova']).items()}
        df_anova = pd.concat([df_anova, pd.DataFrame(d_swap, index=[group]).T], axis=1)
        d_swap = {v: k for k, v in dict(pval_list['Kruskal-Wallis']).items()}
        df_kruskal = pd.concat([df_kruskal, pd.DataFrame(d_swap, index=[group]).T], axis=1)

        # Loop through applied tests
        for test_key in pval_list.keys():
            if len(pval_list[test_key]) > 0:
                for df_val, attr in pval_list[test_key]:

                    # Save padj values for each pairwise comparision in dataframe
                    dict_all_continuous_attributes[test_key].loc[attr, [group]] = df_val

    df_anova.to_excel(os.path.join(save_folder_tmp, 'Stats.xlsx'), sheet_name='Anova')
    # Save padj values to excel file
    with pd.ExcelWriter(os.path.join(save_folder_tmp, 'Stats.xlsx'),
                        mode="a", engine="openpyxl") as writer:
        df_kruskal.to_excel(writer, sheet_name='Kruskal-Wallis', index=True)

    return dict_all_continuous_attributes


def perform_test_ordinal_category(df_scaled, ordinal_cols, encoded__labels, save_folder):
    save_folder_tmp = os.path.join(save_folder, 'Ordinal')
    os.makedirs(save_folder_tmp, exist_ok=True)

    unique_labels = np.unique(encoded__labels)
    names = [(str(s), 'Rest') for s in unique_labels]

    df_ordinal = df_scaled[ordinal_cols]
    df_temp_ordinal = df_ordinal.copy()

    df_all_ordinal_attributes = pd.DataFrame(index=list(df_ordinal.columns), columns=names)
    df_kruskal = pd.DataFrame()

    # Create 1 vs Rest groups
    for group in names:
        df_temp_ordinal['group'] = list(encoded__labels)
        df_temp_ordinal.loc[df_temp_ordinal['group'] != group[0], 'group'] = 'Rest'
        df_temp_ordinal['group'] = df_temp_ordinal['group'].astype('category')

        data_temp = [df_ordinal.loc[ids, :] for ids in df_temp_ordinal.groupby('group').groups.values()]

        # Compare in a one vs rest fashion
        pval_list = {'Kruskal': []}
        # Loop through variables
        for ordinal_attr in ordinal_cols:
            variable_classes = [df_ind.loc[:, ordinal_attr].values for df_ind in data_temp]

            # Apply Kruskal-Wallis
            # Perform non-parametric test eg Kruskal-Wallis
            _, p_temp = perform_test_nonparametric.perform_kruskal_wallis_test(
                data=variable_classes, features=None, predictor=None)
            pval_list['Kruskal'].append((p_temp, ordinal_attr))

        # Save Kruskal-Wallis test results to dataframe
        d_swap = {v: k for k, v in dict(pval_list['Kruskal']).items()}
        df_kruskal = pd.concat([df_kruskal, pd.DataFrame(d_swap, index=['Kruskal pval']).T], axis=1)

        # Using the padj values -> find out in which class the attribute is significantly different
        test_key = "Kruskal"
        if len(pval_list[test_key]) > 0:
            for df_pval, attr in pval_list[test_key]:
                # Save padj values for each pairwise comparision in dataframe
                df_all_ordinal_attributes.loc[attr, [group]] = df_pval

    df_kruskal.to_excel(os.path.join(save_folder_tmp, 'Stats.xlsx'), sheet_name='Kruskal-Wallis')

    return df_all_ordinal_attributes


def perform_test_nominal_category(df, nominal_cols, encoded__labels, save_folder):
    save_folder_tmp = os.path.join(save_folder, 'Nominal')
    os.makedirs(save_folder_tmp, exist_ok=True)
    nominal_cols_encoded = df.columns[df.columns.str.contains("|".join(nominal_cols))]

    unique_labels = np.unique(encoded__labels)
    names = [(str(s), 'Rest') for s in unique_labels]

    df_nominal = df.loc[:, nominal_cols_encoded]
    df_temp_nominal = df_nominal.copy()

    df_all_nominal_attributes = pd.DataFrame(index=list(df_nominal.columns), columns=names)
    df_chisquare = pd.DataFrame()

    # Create 1 vs Rest groups
    for group in names:
        df_temp_nominal['group'] = list(encoded__labels)
        df_temp_nominal.loc[df_temp_nominal['group'] != group[0], 'group'] = 'Rest'
        df_temp_nominal['group'] = df_temp_nominal['group'].astype('category')

        data_temp = [df_nominal.loc[ids, :] for ids in df_temp_nominal.groupby('group').groups.values()]

        # Compare in a one vs rest fashion
        pval_list = {'Chi-square': []}
        # Loop through variables
        for nominal_attr in nominal_cols_encoded:
            variable_classes = [df_ind.loc[:, nominal_attr].values for df_ind in data_temp]

            # plot class distributions
            res = dict(zip(np.arange(0, len(np.unique(encoded__labels))), variable_classes))
            df_attr_classes = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in res.items()]))
            df_melted = df_attr_classes.melt()

            # Remove nan values from dataframe
            df_melted_tmp = df_melted.loc[~df_melted.loc[:, 'value'].isna(), :]

            # create a contingency table of the nominal category and class
            cont_table = pd.crosstab(df_melted_tmp['variable'], df_melted_tmp['value'])

            # perform chi-square test of independence
            chi2, p_temp, dof, expected = chi2_contingency(cont_table)

            pval_list['Chi-square'].append((p_temp, nominal_attr))

        # Save Chi-square test results to dataframe
        d_swap = {v: k for k, v in dict(pval_list['Chi-square']).items()}
        df_chisquare = pd.concat([pd.DataFrame(d_swap, index=['Chi-square pval']).T], axis=1)

        # Using the padj values -> find out in which class the attribute is significantly different
        test_key = "Chi-square"
        if len(pval_list[test_key]) > 0:
            # loop through padj list and read out df containing padj values and attribute
            for pvals, attr in pval_list[test_key]:
                # Save padj values for each pairwise comparision in dataframe
                df_all_nominal_attributes.loc[attr, [group]] = pvals

    df_chisquare.to_excel(os.path.join(save_folder_tmp, 'Stats.xlsx'), sheet_name='Chi-square')

    return df_all_nominal_attributes


def plot_summary(df, save_folder, key, figsize=(16, 6)):
    df_tmp_ = -np.log10(df.dropna(axis='index', how='all').astype(float))
    cmap_reds = plt.get_cmap('Reds')
    num_colors = 40  # int(np.ceil(np.amax(df_tmp_) / (-np.log10(0.05))))
    colors = ['lightgrey'] * (num_colors - 1) + [cmap_reds(i / num_colors) for i in range(4, num_colors + 3)]
    cmap = LinearSegmentedColormap.from_list('', colors, len(colors))

    if df_tmp_.shape[0] > 0:
        fig, ax = plt.subplots(figsize=figsize)
        sns.heatmap(df_tmp_, vmin=0, xticklabels=True, yticklabels=True, ax=ax,
                    cmap=cmap, center=-np.log10(0.05), cbar=False, linewidths=1.3)
        ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
        cbar = plt.colorbar(ax.collections[0], ticks=np.arange(0, np.amax(df_tmp_), -np.log10(0.05)))
        cbar.ax.get_yaxis().labelpad = 15
        cbar.ax.set_ylabel('-log10(p-values)', rotation=90, fontsize=12)
        plt.tight_layout()

        plt.savefig(os.path.join(save_folder, "Summary_Heatmap_{}.pdf".format(key)))

        plt.close(fig=fig)


def plot_clustered_summary(df, save_folder, key, pval_cut=0.05, figsize=(16, 12)):
    df = df.dropna(axis='index', how='all').astype(float)
    df_copy = df.copy()
    df_copy = df_copy.dropna(axis='index', how='all')
    df_copy[df_copy > pval_cut] = np.nan
    df_copy = df_copy.dropna(axis='index', how='all')
    # Keep only those rows which dont only have NaN values
    index = list(df_copy.index)
    # drop non-significant drows
    df = df.loc[index]
    df_tmp_ = -np.log10(df)
    cmap_reds = plt.get_cmap('Reds')
    num_colors = 40  # int(np.ceil(np.amax(df_tmp_) / (-np.log10(0.05))))
    colors = ['lightgrey'] * (num_colors - 1) + [cmap_reds(i / num_colors) for i in range(4, num_colors + 3)]
    cmap = LinearSegmentedColormap.from_list('', colors, len(colors))

    df_colors = pd.DataFrame(
        data=np.asarray([["#006ddb"], ["#b6dbff"], ["#004949"], ["#009292"], ["#ff6db6"], ["#490092"],
                         ["#b66dff"], ["#000000"], ["#920000"], ["#E69F00"], ["#D55E00"], ["#8B4513"],
                         ["#999999"]]).T, columns=df.columns)

    df_tmp_ = df_tmp_.T

    # TODO swap x and y axis
    # Compute and plot first dendrogram. (left)
    fig = plt.figure(figsize=figsize)  # (w, h)
    # Top dendrogram
    ax1 = fig.add_axes([0.115, 0.1, 0.08, 0.7])  # (left, bottom, width, height) fig.add_axes([0.2, 0.825, 0.5, 0.08])
    # compute distance and cluster rows (comparisons)
    c_dist = pdist(df_tmp_)
    Y = sch.linkage(c_dist, method='ward')
    Z1 = sch.dendrogram(Y, orientation='left')
    ax1.set_xticks([])
    ax1.set_yticks([])
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['bottom'].set_visible(False)
    ax1.spines['left'].set_visible(False)

    # Compute and plot second dendrogram. (top)
    ax2 = fig.add_axes([0.2, 0.805, 0.6, 0.08])  # (left, bottom, width, height) fig.add_axes([0.095, 0.1, 0.1, 0.7])
    # compute distance and cluster columns (samples)
    c_dist = pdist(df_tmp_.T)
    Y = sch.linkage(c_dist, method='ward')
    Z2 = sch.dendrogram(Y)
    ax2.set_xticks([])
    ax2.set_yticks([])
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.spines['bottom'].set_visible(False)
    ax2.spines['left'].set_visible(False)

    idx1 = Z2['leaves']
    df_tmp_clustered = df_tmp_.iloc[:, idx1]
    idx2 = Z1['leaves']
    df_tmp_clustered = df_tmp_clustered.iloc[idx2, :]

    axmatrix = fig.add_axes([0.2, 0.1, 0.6, 0.7])  # fig.add_axes([0.2, 0.1, 0.5, 0.72])  # (left, bottom, width, height)
    sns.heatmap(df_tmp_clustered, vmin=0, xticklabels=True, yticklabels=True, ax=axmatrix,
                cmap=cmap, center=-np.log10(0.05), cbar=False, linewidths=1.3)
    axmatrix.set_yticklabels(axmatrix.get_yticklabels(), rotation=0)
    axmatrix.set_xticklabels(axmatrix.get_xticklabels(), rotation=90)
    axmatrix.set_yticklabels([])
    axmatrix.yaxis.set_label_position("right")
    axmatrix.yaxis.tick_right()

    # Add colorbar
    axcolor = fig.add_axes([0.82, 0.1, 0.01, 0.7])  # (left, bottom, width, height)
    cbar = plt.colorbar(mappable=axmatrix.collections[0], cax=axcolor,
                        ticks=np.arange(0, np.amax(df_tmp_clustered), -np.log10(0.05)))
    cbar.ax.get_yaxis().labelpad = 15
    cbar.ax.set_ylabel('-log10(p-values)', rotation=90, fontsize=12)
    axmatrix.spines['top'].set_visible(True)
    axmatrix.spines['right'].set_visible(True)
    axmatrix.spines['bottom'].set_visible(True)
    axmatrix.spines['left'].set_visible(True)

    # Add heatmap showing the colors of the clusters
    ax_colors = fig.add_axes([0.802, 0.1, 0.01, 0.7])  # (left, bottom, width, height)
    # Reorder to HC cluster
    df_colors = df_colors.iloc[:, idx2]
    cmap_colors = LinearSegmentedColormap.from_list('', df_colors.iloc[0, :].values, len(df_colors.iloc[0, :].values))
    gradient = np.linspace(0, 1, len(df_colors.columns))
    gradient = np.vstack((gradient, gradient))
    ax_colors.imshow(gradient.T, aspect='auto', cmap=cmap_colors)
    ax_colors.get_xticklabels([])
    ax_colors.get_yticklabels([])
    ax_colors.get_xaxis().set_visible(False)
    ax_colors.get_yaxis().set_visible(False)

    plt.savefig(os.path.join(save_folder, "Summary_clustered_Heatmap_{}_rotated.pdf".format(key)), bbox_inches='tight')
    plt.close()


def plot_clustered_comparisons_summary(df, save_folder, key, xlabels, pval_cut=0.05, figsize=(16, 12)):
    df = df.dropna(axis='index', how='all').astype(float)
    df_copy = df.copy()
    df_copy = df_copy.dropna(axis='index', how='all')
    df_copy[df_copy > pval_cut] = np.nan
    df_copy = df_copy.dropna(axis='index', how='all')
    # Keep only those rows which dont only have NaN values
    index = list(df_copy.index)
    # drop non-significant drows
    df = df.loc[index]
    df_tmp_ = -np.log10(df)
    cmap_reds = plt.get_cmap('Reds')
    num_colors = 40  # int(np.ceil(np.amax(df_tmp_) / (-np.log10(0.05))))
    colors = ['lightgrey'] * (num_colors - 1) + [cmap_reds(i / num_colors) for i in range(4, num_colors + 3)]
    cmap = LinearSegmentedColormap.from_list('', colors, len(colors))

    df_colors = pd.DataFrame(
        data=np.asarray([["#ff6db6"], ["#490092"],
                         ["#b66dff"], ["#000000"], ["#920000"], ["#E69F00"], ["#D55E00"], ["#8B4513"],
                         ["#999999"], ["#006ddb"], ["#b6dbff"], ["#004949"], ["#009292"]]).T, columns=df.columns)

    # Dendrogram
    # compute distance and cluster columns (samples)
    c_dist = pdist(df_tmp_)
    Y = sch.linkage(c_dist, method='ward')

    # Compute and plot first dendrogram. (left)
    fig = plt.figure(figsize=figsize)  # (w, h)
    # Compute and plot second dendrogram. (left)
    ax2 = fig.add_axes([0.16, 0.1, 0.03, 0.72])  # (left, bottom, width, height) fig.add_axes([0.095, 0.1, 0.1, 0.7])
    # Temporarily override the default line width:
    with plt.rc_context({'lines.linewidth': 0.6}):
        Z1 = sch.dendrogram(Y, orientation='left', color_threshold=0, above_threshold_color='k')
        # sch.dendrogram(Y, orientation='left', color_threshold=0, above_threshold_color='k',
        #                truncate_mode='lastp',  p=5, show_contracted=True)
    ax2.set_xticks([])
    ax2.set_yticks([])
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.spines['bottom'].set_visible(False)
    ax2.spines['left'].set_visible(False)

    idx1 = Z1['leaves']
    df_tmp_clustered = df_tmp_.iloc[idx1, :]

    # Split heatmap into Pattern 2a-like, 3-like, and 1-like columns
    # Pattern 2a-like
    df_pattern_2a = df_tmp_clustered[[('E5', 'Rest'),  ('E6', 'Rest'),  ('E7', 'Rest'),
                                      ('E8', 'Rest'), ('E9', 'Rest'), ('E10', 'Rest')]]
    axmatrix = fig.add_axes([0.2, 0.1, 0.28, 0.72])  # (left, bottom, width, height)
    sns.heatmap(df_pattern_2a,
                vmin=0, xticklabels=True, yticklabels=True, ax=axmatrix,
                cmap=cmap, center=-np.log10(0.05), cbar=False, linewidths=1.3)
    axmatrix.set_yticklabels(axmatrix.get_yticklabels(), rotation=0)
    axmatrix.set_xticklabels(axmatrix.get_xticklabels(), rotation=90)
    axmatrix.set_yticklabels([])
    axmatrix.set_yticks([])
    axmatrix.set_xticklabels(xlabels[:6], rotation=0)
    axmatrix.spines['top'].set_visible(True)
    axmatrix.spines['right'].set_visible(True)
    axmatrix.spines['bottom'].set_visible(True)
    axmatrix.spines['left'].set_visible(True)

    # Pattern 3-like
    df_pattern_3 = df_tmp_clustered[[('E11', 'Rest'), ('E12', 'Rest'), ('E13', 'Rest')]]
    axmatrix = fig.add_axes([0.50, 0.1, 0.14, 0.72])  # (left, bottom, width, height)
    sns.heatmap(df_pattern_3, vmin=0, xticklabels=True, yticklabels=True, ax=axmatrix,
                cmap=cmap, center=-np.log10(0.05), cbar=False, linewidths=1.3)
    axmatrix.set_yticklabels(axmatrix.get_yticklabels(), rotation=0)
    axmatrix.set_xticklabels(axmatrix.get_xticklabels(), rotation=90)
    axmatrix.set_yticklabels([])
    axmatrix.set_yticks([])
    axmatrix.set_xticklabels(xlabels[6:9], rotation=0)
    axmatrix.spines['top'].set_visible(True)
    axmatrix.spines['right'].set_visible(True)
    axmatrix.spines['bottom'].set_visible(True)
    axmatrix.spines['left'].set_visible(True)

    # Pattern 1-like
    df_pattern_1 = df_tmp_clustered[[('E1', 'Rest'),  ('E2', 'Rest'),  ('E3', 'Rest'), ('E4', 'Rest')]]
    axmatrix = fig.add_axes([0.66, 0.1, 0.19, 0.72])  # (left, bottom, width, height)
    sns.heatmap(df_pattern_1, vmin=0, xticklabels=True, yticklabels=True, ax=axmatrix,
                cmap=cmap, center=-np.log10(0.05), cbar=False, linewidths=1.3)
    axmatrix.set_yticklabels(axmatrix.get_yticklabels(), rotation=0)
    axmatrix.set_xticklabels(axmatrix.get_xticklabels(), rotation=90)
    axmatrix.yaxis.set_label_position("right")
    axmatrix.yaxis.tick_right()
    axmatrix.set_xticklabels(xlabels[9:], rotation=0)
    axmatrix.spines['top'].set_visible(True)
    axmatrix.spines['right'].set_visible(True)
    axmatrix.spines['bottom'].set_visible(True)
    axmatrix.spines['left'].set_visible(True)

    # Add colorbar
    axcolor = fig.add_axes([0.2, 0.05, 0.65, 0.02])  # (left, bottom, width, height) fig.add_axes([0.9, 0.1, 0.03, 0.72])
    cbar = plt.colorbar(mappable=axmatrix.collections[0], cax=axcolor,
                        ticks=[0, -np.log10(0.05)] + list(
                            np.arange(4, np.amax(df_tmp_clustered), 3)) + [np.amax(df_tmp_clustered)],
                        orientation="horizontal")
    # cbar.ax.get_xaxis().labelpad = 15
    cbar.ax.set_xlabel('-log10(p-values)', rotation=0, fontsize=12)

    # # Add heatmap showing the colors of the clusters
    # ax_colors = fig.add_axes([0.2, 0.825, 0.5, 0.02])  # (left, bottom, width, height)
    # cmap_colors = LinearSegmentedColormap.from_list('', df_colors.iloc[0, :].values, len(df_colors.iloc[0, :].values))
    # gradient = np.linspace(0, 1, len(df_colors.columns))
    # gradient = np.vstack((gradient, gradient))
    # ax_colors.imshow(gradient, aspect='auto', cmap=cmap_colors)
    # ax_colors.get_xticklabels([])
    # ax_colors.get_yticklabels([])
    # ax_colors.get_xaxis().set_visible(False)
    # ax_colors.get_yaxis().set_visible(False)

    plt.savefig(os.path.join(save_folder, "Summary_clustered_attr_Heatmap_{}_pval_{}.pdf".format(key, pval_cut)),
                bbox_inches='tight')
    plt.close()


def plot_clustered_samples_summary(df, save_folder, key, figsize=(16, 10)):
    df_tmp_ = -np.log10(df.dropna(axis='index').astype(float))
    cmap_reds = plt.get_cmap('Reds')
    num_colors = 40  # int(np.ceil(np.amax(df_tmp_) / (-np.log10(0.05))))
    colors = ['lightgrey'] * (num_colors - 1) + [cmap_reds(i / num_colors) for i in range(4, num_colors + 3)]
    cmap = LinearSegmentedColormap.from_list('', colors, len(colors))

    # Compute and plot first dendrogram.
    fig = plt.figure(figsize=figsize)
    ax1 = fig.add_axes([0.1, 0.81, 0.645, 0.1])  # (left, bottom, width, height)
    # compute distance and cluster columns
    c_dist = pdist(df_tmp_.T)
    Y = sch.linkage(c_dist, method='ward')
    Z1 = sch.dendrogram(Y)
    ax1.set_xticks([])
    ax1.set_yticks([])
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['bottom'].set_visible(False)
    ax1.spines['left'].set_visible(False)

    idx1 = Z1['leaves']
    df_tmp_clustered = df_tmp_.iloc[:, idx1]

    axmatrix = fig.add_axes([0.1, 0.1, 0.8, 0.7])  # (left, bottom, width, height)
    sns.heatmap(df_tmp_clustered, vmin=0, xticklabels=True, yticklabels=True, ax=axmatrix,
                cmap=cmap, center=-np.log10(0.05), cbar=False, linewidths=1.3)
    axmatrix.set_yticklabels(axmatrix.get_yticklabels(), rotation=0)
    axmatrix.set_xticklabels(axmatrix.get_xticklabels(), rotation=90)
    cbar = plt.colorbar(axmatrix.collections[0], ticks=np.arange(0, np.amax(df_tmp_clustered), -np.log10(0.05)))
    cbar.ax.get_yaxis().labelpad = 15
    cbar.ax.set_ylabel('-log10(p-values)', rotation=90, fontsize=12)
    sns.despine(fig=fig)
    plt.savefig(os.path.join(save_folder, "Summary_clustered_comparisons_Heatmap_{}.pdf".format(key)),
                bbox_inches='tight')
    plt.close()


def plot_selected_attributes(df, save_folder, key, xlabels, pval_cut=0.05, figsize=(16, 12)):
    df = df.dropna(axis='index', how='all').astype(float)
    df_copy = df.copy()
    df_copy = df_copy.dropna(axis='index', how='all')
    df_copy[df_copy > pval_cut] = np.nan
    df_copy = df_copy.dropna(axis='index', how='all')
    # Keep only those rows which don't only have NaN values
    index = list(df_copy.index)
    # drop non-significant drows
    df = df.loc[index]
    df_tmp_ = -np.log10(df)
    cmap_reds = plt.get_cmap('Reds')
    num_colors = 40
    colors = ['lightgrey'] * (num_colors - 1) + [cmap_reds(i / num_colors) for i in range(4, num_colors + 3)]
    cmap = LinearSegmentedColormap.from_list('', colors, len(colors))

    df_colors = pd.DataFrame(
        data=np.asarray([["#ff6db6"], ["#490092"],
                         ["#b66dff"], ["#000000"], ["#920000"], ["#E69F00"], ["#D55E00"], ["#8B4513"],
                         ["#999999"], ["#006ddb"], ["#b6dbff"], ["#004949"], ["#009292"]]).T, columns=df.columns)

    fig = plt.figure(figsize=figsize)  # (w, h)
    # Compute and plot second dendrogram. (left)
    ax2 = fig.add_axes([0.16, 0.1, 0.03, 0.72])  # (left, bottom, width, height) fig.add_axes([0.095, 0.1, 0.1, 0.7])
    ax2.set_xticks([])
    ax2.set_yticks([])
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.spines['bottom'].set_visible(False)
    ax2.spines['left'].set_visible(False)

    # Split heatmap into Pattern 2a-like, 3-like, and 1-like columns
    # Pattern 2a-like
    df_pattern_2a = df_tmp_[[('E5', 'Rest'),  ('E6', 'Rest'),  ('E7', 'Rest'),
                                      ('E8', 'Rest'), ('E9', 'Rest'), ('E10', 'Rest')]]
    axmatrix = fig.add_axes([0.2, 0.1, 0.28, 0.72])  # (left, bottom, width, height)
    sns.heatmap(df_pattern_2a,
                vmin=0, xticklabels=True, yticklabels=True, ax=axmatrix,
                cmap=cmap, center=-np.log10(0.05), cbar=False, linewidths=1.3)
    axmatrix.set_yticklabels(axmatrix.get_yticklabels(), rotation=0)
    axmatrix.set_xticklabels(axmatrix.get_xticklabels(), rotation=90)
    axmatrix.set_yticklabels([])
    axmatrix.set_yticks([])
    axmatrix.set_xticklabels(xlabels[:6], rotation=0)
    axmatrix.spines['top'].set_visible(True)
    axmatrix.spines['right'].set_visible(True)
    axmatrix.spines['bottom'].set_visible(True)
    axmatrix.spines['left'].set_visible(True)

    # Pattern 3-like
    df_pattern_3 = df_tmp_[[('E11', 'Rest'), ('E12', 'Rest'), ('E13', 'Rest')]]
    axmatrix = fig.add_axes([0.50, 0.1, 0.14, 0.72])  # (left, bottom, width, height)
    sns.heatmap(df_pattern_3, vmin=0, xticklabels=True, yticklabels=True, ax=axmatrix,
                cmap=cmap, center=-np.log10(0.05), cbar=False, linewidths=1.3)
    axmatrix.set_yticklabels(axmatrix.get_yticklabels(), rotation=0)
    axmatrix.set_xticklabels(axmatrix.get_xticklabels(), rotation=90)
    axmatrix.set_yticklabels([])
    axmatrix.set_yticks([])
    axmatrix.set_xticklabels(xlabels[6:9], rotation=0)
    axmatrix.spines['top'].set_visible(True)
    axmatrix.spines['right'].set_visible(True)
    axmatrix.spines['bottom'].set_visible(True)
    axmatrix.spines['left'].set_visible(True)

    # Pattern 1-like
    df_pattern_1 = df_tmp_[[('E1', 'Rest'),  ('E2', 'Rest'),  ('E3', 'Rest'), ('E4', 'Rest')]]
    axmatrix = fig.add_axes([0.66, 0.1, 0.19, 0.72])  # (left, bottom, width, height)
    sns.heatmap(df_pattern_1, vmin=0, xticklabels=True, yticklabels=True, ax=axmatrix,
                cmap=cmap, center=-np.log10(0.05), cbar=False, linewidths=1.3)
    axmatrix.set_yticklabels(axmatrix.get_yticklabels(), rotation=0)
    axmatrix.set_xticklabels(axmatrix.get_xticklabels(), rotation=90)
    axmatrix.yaxis.set_label_position("right")
    axmatrix.yaxis.tick_right()
    axmatrix.set_xticklabels(xlabels[9:], rotation=0)
    axmatrix.spines['top'].set_visible(True)
    axmatrix.spines['right'].set_visible(True)
    axmatrix.spines['bottom'].set_visible(True)
    axmatrix.spines['left'].set_visible(True)

    # Add colorbar
    axcolor = fig.add_axes([0.2, 0.02, 0.65, 0.02])  # (left, bottom, width, height) fig.add_axes([0.9, 0.1, 0.03, 0.72])
    cbar = plt.colorbar(mappable=axmatrix.collections[0], cax=axcolor,
                        ticks=[0, -np.log10(0.05)] + list(
                            np.arange(4, np.amax(df_tmp_), 3)) + [np.amax(df_tmp_)],
                        orientation="horizontal")
    # cbar.ax.get_xaxis().labelpad = 15
    cbar.ax.set_xlabel('-log10(p-values)', rotation=0, fontsize=12)

    # # Add heatmap showing the colors of the clusters - TODO split into three columns
    # ax_colors = fig.add_axes([0.2, 0.83, 0.65, 0.02])  # (left, bottom, width, height)
    # cmap_colors = LinearSegmentedColormap.from_list('', df_colors.iloc[0, :].values, len(df_colors.iloc[0, :].values))
    # gradient = np.linspace(0, 1, len(df_colors.columns))
    # gradient = np.vstack((gradient, gradient))
    # ax_colors.imshow(gradient, aspect='auto', cmap=cmap_colors)
    # ax_colors.get_xticklabels([])
    # ax_colors.get_yticklabels([])
    # ax_colors.get_xaxis().set_visible(False)
    # ax_colors.get_yaxis().set_visible(False)

    plt.savefig(os.path.join(save_folder, "Summary_attr_Heatmap_{}_pval_{}.pdf".format(key, pval_cut)),
                bbox_inches='tight')
    plt.close()


def main(input_dir, input_meta_dir, save_folder):
    adata = sc.read(os.path.join(
        input_dir,
        'Lesion_RNAseq_20210720_patient_meta_data_v04__CH_20220209_TherapyResponse_v04_v220715_TherapyResponseCorrected__Endotypes_230605.h5'))

    obs_name = 'Endotypes'
    pval_cut = 0.05
    figsize = (8, 12)

    # Encode labels
    if obs_name == 'Pattern':
        encoded__labels = encode_labels(adata)
    else:
        encoded__labels = adata.obs[obs_name].values

    # Attribute preparation
    df, nominal_cols, ordinal_cols, continuous_cols = prepare_attributes(
        adata=adata, input_dir_attr=input_meta_dir, save_folder=save_folder)

    # 1 vs Rest
    # 1. Continuous variables
    dict_all_continuous_attributes = perform_test_continuous_variables(
        df=df, continuous_cols=continuous_cols, encoded__labels=encoded__labels, save_folder=save_folder)

    # 2. Ordinal (Dichotomus) categories
    df_all_ordinal_attributes = perform_test_ordinal_category(
        df_scaled=df, ordinal_cols=ordinal_cols, encoded__labels=encoded__labels, save_folder=save_folder)

    # 3. Nominal categories
    df_all_nominal_attributes = perform_test_nominal_category(
        df=df, nominal_cols=nominal_cols, encoded__labels=encoded__labels, save_folder=save_folder)

    # Plot summary plot
    df_all_continuous_attributes = pd.DataFrame()
    for key_name in dict_all_continuous_attributes:
        dict_all_continuous_attributes = pd.concat(
            [df_all_continuous_attributes, dict_all_continuous_attributes[key_name]])

        plot_summary(df=dict_all_continuous_attributes[key_name], save_folder=save_folder,
                     key='continuous_{}'.format(key_name))

    plot_summary(df=df_all_ordinal_attributes, save_folder=save_folder, key='ordinal')
    plot_summary(df=df_all_nominal_attributes, save_folder=save_folder, key='nominal')

    # concatenate dfs
    df_merged = pd.concat([df_all_continuous_attributes, df_all_ordinal_attributes,
                           df_all_nominal_attributes], axis=0)
    # reorder columns like sorted in dendrogram
    adata.obs['Endotypes'] = adata.obs['Endotypes'].cat.reorder_categories(
        ['E5', 'E6', 'E7', 'E8', 'E9', 'E10', 'E11', 'E12', 'E13', 'E1', 'E2', 'E3', 'E4'])
    sorted_colnames = [(str(s), 'Rest') for s in list(adata.obs[obs_name].cat.categories)]
    df_merged = df_merged.loc[:, sorted_colnames]

    plot_summary(df=df_merged, save_folder=save_folder, key='merged', figsize=(12, 8))
    # Apply unsupervised clustering on clinical attributes - which attributes do cluster together?
    plot_clustered_summary(df=df_merged, save_folder=save_folder, pval_cut=pval_cut, key='merged', figsize=(18, 12))
    plot_clustered_comparisons_summary(df=df_merged, save_folder=save_folder,
                                       xlabels=list(adata.obs[obs_name].cat.categories),
                                       pval_cut=pval_cut, key='merged', figsize=figsize)

    # Create new heatmap plot without clustering of rows
    plot_selected_attributes(df=df_merged, save_folder=save_folder,
                             xlabels=list(adata.obs[obs_name].cat.categories),
                             pval_cut=pval_cut, key='all', figsize=figsize)

    df_merged.to_excel(os.path.join(save_folder, 'Figure_S4J_Heatmap_enriched_clinical_attributes_Endotypes.xlsx'))


if __name__ == '__main__':
    general_dir = os.path.join(
        '/Volumes', 'CH__data', 'Projects', 'Eyerich_AG_projects', 'BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer')
    output_dir = os.path.join(general_dir, 'analysis', 'Molecular_subtypes', 'output', 'Figure_S4', str(date.today()))
    os.makedirs(output_dir, exist_ok=True)

    data_root = os.path.join(general_dir, 'analysis', 'Molecular_subtypes', 'input', 'h5_files', 'LESION')

    meta_dir = os.path.join(general_dir, 'raw_data', 'clinical_data')

    main(save_folder=output_dir, input_dir=data_root, input_meta_dir=meta_dir)
