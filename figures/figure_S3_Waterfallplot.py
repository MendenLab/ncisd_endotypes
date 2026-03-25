"""
@description: Find enriched clinical attributes in Endotypes,
distinguish between categorical nominal, categorical ordinal and Continuous attributes

@author : Christina Hillig
"""

from scripts.feature_engineering import imputation, encoding
from scripts.utils import check_distributions, perform_test_nonparametric

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.collections import PatchCollection

from sklearn import preprocessing
import seaborn as sns
import scanpy as sc
import pandas as pd
import numpy as np
import itertools
import collections

import os
from datetime import date

from scipy.stats import ttest_ind, f_oneway, chi2_contingency
from statsmodels.sandbox.stats.multicomp import multipletests

plt.rcParams["font.family"] = "Arial"


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
    # drop columns which are not clinical traits
    drop_cols_other = ['Therapeutic target']
    nominal_cols = [ele for ele in nominal_cols if ele not in drop_cols_other]
    df_nominals = bulk_data.obs[nominal_cols]  # .drop(drop_cols_other, axis=1)
    # 2.2.1.2 Apply encoding -  dropped Com_uv as attribute was uniform
    df_encoded_nominals = encoding.encode_nominal_feature(
        data=df_nominals, dummy_variable_encoding=True)  # 43
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
    # drop columns which are not clinical traits
    drop_cols_other = ['Diag_PASI', 'Diag_SCORAD', 'Diag_PGA', 'Diag_DLQI', 'The_sor', 'The_dor']
    ordinal_cols = [ele for ele in ordinal_cols if ele not in drop_cols_other]
    df_ordinals = bulk_data.obs[ordinal_cols]  # .drop(drop_cols_other, axis=1)
    # Apply ordinal encoding on the ordinal features
    df_encoded_ordinals = encoding.encode_ordinal_categories(data=df_ordinals.astype(float))
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
                            'The_crea', 'The_GPT'])
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

    # 2.2.6 Remove categories with single category value such as 'Com_uv' -> no longer the case
    drop_cols = [temp_col for temp_col in df_encoded_features.columns
                 if len(df_encoded_features[temp_col].dropna().unique()) <= 1]
    # df_encoded_features = feature_engineering.remove_lowvariance_features(df_data=df_encoded_features)
    df_encoded_features = df_encoded_features.drop(drop_cols, axis=1)

    print('Number of encoded attributes: ', df_encoded_features.shape[1])

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

    return (df_encoded_features, nominal_cols, ordinal_cols, continuous_cols,
            list(df_encoded_nominals.columns), list(df_encoded_ordinals.columns))


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
        os.path.join(
            '/Volumes/CH__data/Projects/Eyerich_AG_projects/BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer',
            'raw_data', 'clinical_data', "patient_meta_data_final.xlsx"), sheet_name='legend')

    meta_data_data_type = pd.read_excel(
        os.path.join(
            '/Volumes/CH__data/Projects/Eyerich_AG_projects/BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer',
            'raw_data', 'clinical_data', "Tables.xlsx"), sheet_name='Table_S1')

    # Encode data
    (df_encoded,  nominal_cols, ordinal_cols, continuous_cols,
     nominal_encoded_cols, ordinal_encoded_cols) = encode_impute_data(
        bulk_data=adata, meta_data_legend=meta_data_legend,
        meta_data_data_type=meta_data_data_type, save_folder=save_folder)

    print('Number of variables after encoding: ', df_encoded.shape[1])

    # Drop columns with more than 70% missing values
    df_percentage = df_encoded.isnull().mean() * 100
    # df_encoded.isnull().sum().sort_values(ascending=False)/len(df_encoded)*100
    # Test if variables are missing completely at random (MCAR)
    df_tmp = df_encoded.copy()
    cat_cols = nominal_encoded_cols + ordinal_encoded_cols
    padj = []
    df_cat_cat = pd.DataFrame()
    for col in cat_cols:
        # Initialize the Dice columns
        df_tmp["{}_missing".format(col)] = df_encoded[col]
        # The column is false
        df_tmp["{}_missing".format(col)] = False
        # Replace where Height_missing with True where Height is missing
        df_tmp.loc[df_tmp[df_tmp[col].isnull()].index, "{}_missing".format(col)] = True

        # Chi2 to test for independence between categorical variables
        pvals = []
        cat_cat_cols = []
        for col_other in df_encoded[cat_cols].columns[(df_encoded[cat_cols].columns != col)]:
            true_counts = df_tmp[df_tmp["{}_missing".format(col)] == True].groupby(col_other)["{}_missing".format(col)].count()
            false_counts = df_tmp[df_tmp["{}_missing".format(col)] == False].groupby(col_other)["{}_missing".format(col)].count()
            # Chi2 square assumes at least a frequency 5 in a cell
            try:
                if not ((len(true_counts[true_counts < 1]) > 0) or (len(false_counts[false_counts < 1]) > 0)) and (
                        np.nansum(df_tmp[col_other].unique()) == sum(list(true_counts.index))):
                    table = [true_counts.astype(int).values, false_counts.astype(int).values]
                    chi2, p, dof, ex = chi2_contingency(table)

                    pvals.append(p)
                    cat_cat_cols.append(("{} vs. {}".format(col, col_other)))
            except ValueError:
                # print(col_other)
                # Contingency table (observed frequencies)
                contingency_table = [true_counts.astype(int).values.tolist(), false_counts.astype(int).values.tolist()]
                # Pad the rows with zeros to ensure uniform dimensions
                max_len = max(len(row) for row in contingency_table)
                contingency_table_padded = [row + [0] * (max_len - len(row)) for row in contingency_table]

                # Transpose the contingency table
                contingency_table_transposed = np.array(contingency_table_padded).T
                chi2, p, dof, ex = chi2_contingency(contingency_table_transposed)
                pvals.append(p)
                cat_cat_cols.append(("{} vs. {}".format(col, col_other)))

        if len(pvals) > 0:
            reject_list, corrected_p_vals = multipletests(pvals, method='fdr_bh')[:2]
            padj.extend(list(corrected_p_vals))

            df_cat_cat_tmp = pd.DataFrame.from_dict(
                {'comparison': cat_cat_cols, 'p-value': pvals, 'p-adj': corrected_p_vals})
            df_cat_cat = pd.concat([df_cat_cat, df_cat_cat_tmp])

    print("Relation between categories: ", np.asarray(padj)[np.asarray(padj) < 0.05])
    df_cat_cat.to_excel(os.path.join(save_folder, 'Chi2_cat_cat_dependency_test.xlsx'))

    # test for dependencies between continuous and categorical data
    padj_cat_con = []
    df_cat_con = pd.DataFrame()
    for col in cat_cols:
        pvals_cat_con = []
        cat_con_cols = []
        categorical_data = df_encoded[col]
        for con_col in continuous_cols:
            continuous_data = df_encoded[con_col]
            # Perform one-way ANOVA
            category_groups = []
            for category in categorical_data.unique()[pd.notna(categorical_data.unique())]:
                vals = continuous_data[categorical_data == category]
                mask_nan = vals != 'nan'
                vals = vals[mask_nan]
                vals = vals[~pd.isna(vals)]
                category_groups.append(vals.tolist())
            # Flatten the list of lists into a single list
            flattened_list = [item for sublist in category_groups for item in sublist]
            if len(np.unique(flattened_list)) > 1:
                f_statistic, p_value = f_oneway(*category_groups)
                pvals_cat_con.append(p_value)
                cat_con_cols.append(("{} vs. {}".format(col, con_col)))

        if len(pvals_cat_con) > 0:
            reject_list, corrected_p_vals = multipletests(pvals_cat_con, method='fdr_bh')[:2]
            padj_cat_con.extend(list(corrected_p_vals))
            df_cat_con_tmp = pd.DataFrame.from_dict(
                {'comparison': cat_con_cols, 'p-value': pvals_cat_con, 'p-adj': corrected_p_vals})
            df_cat_con = pd.concat([df_cat_con, df_cat_con_tmp])
    print("Relation between categories and continuous values: ", np.asarray(padj_cat_con)[np.asarray(padj_cat_con) < 0.05])

    # save as table and put in supplements
    df_cat_con.to_excel(os.path.join(save_folder, 'Anova_cat_con_dependency_test.xlsx'))

    variables_to_drop = list(df_percentage[df_percentage > 70].index)
    print('Highest missing data proportion is {:.2f}%'.format(df_percentage.max()))
    print('Attributes to drop due to missing values > 70% ', variables_to_drop)
    df_encoded = df_encoded.drop(variables_to_drop, axis=1)  # 73
    continuous_cols = [col for col in continuous_cols if col not in variables_to_drop]
    nominal_cols = [col for col in nominal_cols if col not in variables_to_drop]
    ordinal_cols = [col for col in ordinal_cols if col not in variables_to_drop]

    print('Number of variables after filtering: ', df_encoded.shape[1])

    # Normalise ordinal and continuous data before imputation as the distance
    # measure in KNN imputation suffers from large values
    # TODO check if I can impute categorical and continuous together
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


def cohens_d(d1, d2):
    """ function to calculate Cohen's d for independent samples

    Parameters
    ----------
    d1
    d2

    Returns
    -------

    """
    # calculate the size of samples
    n1, n2 = len(d1), len(d2)
    # calculate the variance of the samples
    s1, s2 = np.var(d1, ddof=1), np.var(d2, ddof=1)
    # calculate the pooled standard deviation
    s = np.sqrt(((n1 - 1) * s1 + (n2 - 1) * s2) / (n1 + n2 - 2))
    # calculate the means of the samples
    u1, u2 = np.mean(d1), np.mean(d2)
    # calculate the effect size
    return (u1 - u2) / s


def categorize_cohens_d(cohens_d):
    """
    Categorize Cohen's d into effect size categories.

    Parameters:
    cohens_d (float): The Cohen's d value to be categorized.

    Returns:
    str: The category of Cohen's d.
    """
    dict_sign = {'-1.0': 'less', '1.0': 'greater', '0.0': 'equal'}
    if abs(cohens_d) < 0.2:
        return dict_sign[str(np.sign(cohens_d))] + ", Negligible effect"
    elif 0.2 <= abs(cohens_d) < 0.5:
        return dict_sign[str(np.sign(cohens_d))] + ", Small effect"
    elif 0.5 <= abs(cohens_d) < 0.8:
        return dict_sign[str(np.sign(cohens_d))] + ", Medium effect"
    elif abs(cohens_d) >= 0.8:
        return dict_sign[str(np.sign(cohens_d))] + ", Large effect"
    else:
        return "Undefined effect size"


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
    dict_all_ordinal_attributes_effect_size = {
        "Anova": pd.DataFrame(index=list(df_continuous.columns), columns=names),
        "Kruskal-Wallis": pd.DataFrame(index=list(df_continuous.columns), columns=names)}
    df_kruskal = pd.DataFrame()
    df_anova = pd.DataFrame()
    df_effect_size = pd.DataFrame()
    df_effect_size_categorised = pd.DataFrame()
    df_mean_group_0 = pd.DataFrame()
    df_mean_group_rest = pd.DataFrame()
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
        pval_list = {'Anova': [], "Welch's t-test": [], 'Kruskal-Wallis': []}
        effect_size = []
        effect_size_categorised = []
        mean_group_0 = []
        mean_group_rest = []
        # Loop through variables
        for ind_var_normal in df_var_normal_test.index:
            variable_classes = [df_ind.loc[:, ind_var_normal].values for df_ind in data_temp]

            # plot class distributions
            res = dict(zip(np.arange(0, len(np.unique(encoded__labels))), variable_classes))

            # Calculate Cohen's f
            effect_size_tmp = cohens_d(variable_classes[0], variable_classes[1])
            categorised_effect_size_tmp = categorize_cohens_d(cohens_d=effect_size_tmp)

            effect_size.append((ind_var_normal, effect_size_tmp))
            effect_size_categorised.append((ind_var_normal, categorised_effect_size_tmp))
            mean_group_0.append((ind_var_normal, np.mean(variable_classes[0])))
            mean_group_rest.append((ind_var_normal, np.mean(variable_classes[1])))

            # Apply first layer Anova or Kruskal Wallis test
            if df_var_normal_test.loc[ind_var_normal, 'Normal_distribution']:
                # Apply Anova
                if df_var_normal_test.loc[ind_var_normal, 'var_homogeneity_pval'] > 0.05:
                    # Perform oneway-Anova
                    _, p_temp = f_oneway(variable_classes)
                    pval_list['Anova'].append((p_temp, ind_var_normal))
                else:
                    # print('Test for Normal + unequal Variance: Not yet implemented')
                    # Welch one-factor ANOVA with ad hoc Welch t comparisons
                    # Conduct Welch's t-Test and print the result
                    _, p_temp = ttest_ind(a=variable_classes[0], b=variable_classes[1], equal_var=False)
                    pval_list["Welch's t-test"].append((p_temp, ind_var_normal))
            else:
                # Apply Kruskal-Wallis
                # Perform non-parametric test eg Kruskal-Wallis
                _, p_temp = perform_test_nonparametric.perform_kruskal_wallis_test(
                    data=variable_classes, features=None, predictor=None)
                pval_list['Kruskal-Wallis'].append((p_temp, ind_var_normal))

        # Save Anova and Kruskal-Wallis test results to dataframe
        d_swap = {v: k for k, v in dict(pval_list['Anova']).items()}
        df_anova = pd.concat([
            df_anova, pd.DataFrame(d_swap, index=["{} vs {}".format(group[0], group[1])]).T], axis=1)
        d_swap = {v: k for k, v in dict(pval_list['Kruskal-Wallis']).items()}
        df_kruskal = pd.concat([
            df_kruskal, pd.DataFrame(d_swap, index=["{} vs {}".format(group[0], group[1])]).T], axis=1)

        # Save effect size
        df_effect_size = pd.concat([
            df_effect_size,
            pd.DataFrame(dict(effect_size), index=["{} vs {}".format(group[0], group[1])]).T], axis=1)
        df_effect_size_categorised = pd.concat([
            df_effect_size_categorised,
            pd.DataFrame(dict(effect_size_categorised), index=["{} vs {}".format(group[0], group[1])]).T], axis=1)

        df_mean_group_0 = pd.concat([
            df_mean_group_0, pd.DataFrame(dict(mean_group_0), index=[group[0]]).T], axis=1)
        df_mean_group_rest = pd.concat([
            df_mean_group_rest, pd.DataFrame(dict(mean_group_rest), index=["Rest ({})".format(group[0])]).T], axis=1)

        # Loop through applied tests
        for test_key in pval_list.keys():
            if len(pval_list[test_key]) > 0:
                for df_val, attr in pval_list[test_key]:
                    # Save padj values for each pairwise comparison in dataframe
                    dict_all_continuous_attributes[test_key].loc[attr, [group]] = df_val

            if len(pval_list[test_key]) > 0:
                for attr, effect in effect_size:
                    # Save padj values for each pairwise comparision in dataframe
                    dict_all_ordinal_attributes_effect_size[test_key].loc[attr, [group]] = effect

    df_anova.to_excel(os.path.join(save_folder_tmp, 'Stats_Continuous.xlsx'), sheet_name='Anova')
    # Save padj values to excel file
    with pd.ExcelWriter(os.path.join(save_folder_tmp, 'Stats_Continuous.xlsx'),
                        mode="a", engine="openpyxl") as writer:
        df_kruskal.to_excel(writer, sheet_name='Kruskal-Wallis', index=True)
        df_effect_size.to_excel(writer, sheet_name="Cohen's D", index=True)
        df_effect_size_categorised.to_excel(writer, sheet_name="Cohen's D categorised", index=True)
        df_mean_group_0.to_excel(writer, sheet_name='Mean group Endotypes', index=True)
        df_mean_group_rest.to_excel(writer, sheet_name='Mean group Rest', index=True)

    return dict_all_continuous_attributes, dict_all_ordinal_attributes_effect_size


# Calculate Cliff's delta
def cliffs_delta(lst1, lst2, **dull):
    # taken from https://github.com/neilernst/cliffsDelta

    """Returns delta and true if there are more than 'dull' differences"""
    if not dull:
        dull = {'small': 0.147, 'medium': 0.33, 'large': 0.474}  # effect sizes from (Hess and Kromrey, 2004)
    m, n = len(lst1), len(lst2)
    lst2 = sorted(lst2)
    j = more = less = 0
    for repeats, x in runs(sorted(lst1)):
        while j <= (n - 1) and lst2[j] < x:
            j += 1
        more += j*repeats
        while j <= (n - 1) and lst2[j] == x:
            j += 1
        less += (n - j)*repeats
    d = (more - less) / (m*n)
    size = lookup_size(d, dull)
    return d, size


def lookup_size(delta: float, dull: dict = None) -> str:
    """
    :type delta: float
    :type dull: dict, a dictionary of small, medium, large thresholds.
    """
    if dull is None:
        dull = {'small': 0.147, 'medium': 0.33, 'large': 0.474}

    dict_sign = {'-1.0': 'less', '1.0': 'greater', '0.0': 'equal'}
    if abs(delta) < dull['small']:
        return dict_sign[str(np.sign(delta))] + ', Negligible effect'
    if dull['small'] <= abs(delta) < dull['medium']:
        return dict_sign[str(np.sign(delta))] + ', Small effect'
    if dull['medium'] <= abs(delta) < dull['large']:
        return dict_sign[str(np.sign(delta))] + ', Medium effect'
    if abs(delta) >= dull['large']:
        return dict_sign[str(np.sign(delta))] + ', Large effect'


def runs(lst):
    """Iterator, chunks repeated values"""
    for j, two in enumerate(lst):
        if j == 0:
            one, i = two, 0
        if one != two:
            yield j - i, one
            i = j
        one = two
    yield j - i + 1, two


def perform_test_ordinal_category(df_scaled, ordinal_cols, encoded__labels, save_folder):
    save_folder_tmp = os.path.join(save_folder, 'Ordinal')
    os.makedirs(save_folder_tmp, exist_ok=True)

    unique_labels = np.unique(encoded__labels)
    names = [(str(s), 'Rest') for s in unique_labels]

    df_ordinal = df_scaled[ordinal_cols]
    df_temp_ordinal = df_ordinal.copy()

    df_all_ordinal_attributes = pd.DataFrame(index=list(df_ordinal.columns), columns=names)
    df_all_ordinal_attributes_effect_size = pd.DataFrame(index=list(df_ordinal.columns), columns=names)
    df_kruskal = pd.DataFrame()
    df_effect_size = pd.DataFrame()
    df_effect_size_categorised = pd.DataFrame()
    df_median_group_0 = pd.DataFrame()
    df_median_group_1 = pd.DataFrame()
    df_mean_group_0 = pd.DataFrame()
    df_mean_group_1 = pd.DataFrame()

    # Create 1 vs Rest groups
    for group in names:
        df_temp_ordinal['group'] = list(encoded__labels)
        df_temp_ordinal.loc[df_temp_ordinal['group'] != group[0], 'group'] = 'Rest'
        df_temp_ordinal['group'] = df_temp_ordinal['group'].astype('category')

        data_temp = [df_ordinal.loc[ids, :] for ids in df_temp_ordinal.groupby('group').groups.values()]

        # Compare in a one vs rest fashion
        pval_list = {'Kruskal': []}
        effect_size = []
        effect_size_categorised = []
        median_group_0 = []
        median_group_1 = []
        mean_group_0 = []
        mean_group_1 = []
        # Loop through variables
        for ordinal_attr in ordinal_cols:
            variable_classes = [df_ind.loc[:, ordinal_attr].values for df_ind in data_temp]

            # Apply Kruskal-Wallis
            # Perform non-parametric test eg Kruskal-Wallis
            _, p_temp = perform_test_nonparametric.perform_kruskal_wallis_test(
                data=variable_classes, features=None, predictor=None)
            pval_list['Kruskal'].append((p_temp, ordinal_attr))

            # Calculate Cliff's delta
            # Values closer to -1 or 1 indicating stronger association and
            # values closer to 0 indicating weaker association.
            effect_size_tmp, effect_size_categorised_tmp = cliffs_delta(
                lst1=variable_classes[0].astype(float), lst2=variable_classes[1].astype(float))
            effect_size.append((ordinal_attr, effect_size_tmp))
            effect_size_categorised.append((ordinal_attr, effect_size_categorised_tmp))
            median_group_0.append((ordinal_attr, np.median(variable_classes[0].astype(float))))
            median_group_1.append((ordinal_attr, np.median(variable_classes[1].astype(float))))
            mean_group_0.append((ordinal_attr, np.mean(variable_classes[0].astype(float))))
            mean_group_1.append((ordinal_attr, np.mean(variable_classes[1].astype(float))))

        # Save Kruskal-Wallis test results to dataframe
        d_swap = {v: k for k, v in dict(pval_list['Kruskal']).items()}
        df_kruskal = pd.concat([
            df_kruskal, pd.DataFrame(d_swap, index=["{} vs {}".format(group[0], group[1])]).T], axis=1)

        # Save effect size
        df_effect_size = pd.concat(
            [df_effect_size,
             pd.DataFrame(dict(effect_size), index=["{} vs {}".format(group[0], group[1])]).T], axis=1)
        df_effect_size_categorised = pd.concat(
            [df_effect_size_categorised,
             pd.DataFrame(dict(effect_size_categorised), index=["{} vs {}".format(group[0], group[1])]).T], axis=1)

        df_median_group_0 = pd.concat([
            df_median_group_0, pd.DataFrame(dict(median_group_0), index=[group[0]]).T], axis=1)
        df_median_group_1 = pd.concat([
            df_median_group_1, pd.DataFrame(dict(median_group_1), index=["Rest ({})".format(group[0])]).T], axis=1)
        df_mean_group_0 = pd.concat([
            df_mean_group_0, pd.DataFrame(dict(mean_group_0), index=[group[0]]).T], axis=1)
        df_mean_group_1 = pd.concat([
            df_mean_group_1, pd.DataFrame(dict(mean_group_1), index=["Rest ({})".format(group[0])]).T], axis=1)

        # Using the padj values -> find out in which class the attribute is significantly different
        test_key = "Kruskal"
        if len(pval_list[test_key]) > 0:
            for df_pval, attr in pval_list[test_key]:
                # Save padj values for each pairwise comparision in dataframe
                df_all_ordinal_attributes.loc[attr, [group]] = df_pval

        if len(pval_list[test_key]) > 0:
            for attr, effect in effect_size:
                # Save padj values for each pairwise comparision in dataframe
                df_all_ordinal_attributes_effect_size.loc[attr, [group]] = effect

    df_kruskal.to_excel(os.path.join(save_folder_tmp, 'Stats_Ordinal.xlsx'), sheet_name='Kruskal-Wallis')
    with pd.ExcelWriter(os.path.join(save_folder_tmp, 'Stats_Ordinal.xlsx'),
                        mode="a", engine="openpyxl") as writer:
        df_effect_size.to_excel(writer, sheet_name="Cliff's Delta", index=True)
        df_effect_size_categorised.to_excel(writer, sheet_name="Cliff's Delta categorised", index=True)
        df_median_group_0.to_excel(writer, sheet_name='Median group Endotypes', index=True)
        df_median_group_1.to_excel(writer, sheet_name='Median group Rest', index=True)
        df_mean_group_0.to_excel(writer, sheet_name='Mean group Endotypes', index=True)
        df_mean_group_1.to_excel(writer, sheet_name='Mean group Rest', index=True)

    return df_all_ordinal_attributes, df_all_ordinal_attributes_effect_size


def lookup_oddsratio(odds_ratio):
    if odds_ratio < 0.35:
        return 'less, Large effect'
    elif (odds_ratio < 0.70) & (odds_ratio >= 0.35):
        return 'less, Medium effect'
    elif (odds_ratio < 0.92) & (odds_ratio >= 0.70):
        return 'less, Small effect'
    elif (odds_ratio < 1.0) & (odds_ratio >= 0.92):
        return 'less, Negligible effect'
    elif odds_ratio == 1.0:
        return 'equal'
    elif (odds_ratio > 1.0) & (odds_ratio <= 1.05):
        return 'greater, Negligible effect'
    elif (odds_ratio > 1.05) & (odds_ratio <= 1.4):
        return 'greater, Small effect'
    elif (odds_ratio > 1.4) & (odds_ratio <= 2.9):
        return 'greater, Medium effect'
    elif odds_ratio > 2.9:
        return 'greater, Large effect'
    else:
        return 'not defined'


def perform_test_nominal_category(df, nominal_cols, encoded__labels, save_folder):
    save_folder_tmp = os.path.join(save_folder, 'Nominal')
    os.makedirs(save_folder_tmp, exist_ok=True)
    nominal_cols_encoded = df.columns[df.columns.str.contains("|".join(nominal_cols))]

    unique_labels = np.unique(encoded__labels)
    names = [(str(s), 'Rest') for s in unique_labels]

    df_nominal = df.loc[:, nominal_cols_encoded]
    df_temp_nominal = df_nominal.copy()

    df_all_nominal_attributes = pd.DataFrame(index=list(df_nominal.columns), columns=names)
    df_all_nominal_attributes_effect_size = pd.DataFrame(index=list(df_nominal.columns), columns=names)
    df_chisquare = pd.DataFrame()
    df_effect_size = pd.DataFrame()
    df_effect_size_categorised = pd.DataFrame()
    df_mode_group_0 = pd.DataFrame()
    df_mode_group_rest = pd.DataFrame()

    # Create 1 vs Rest groups
    for group in names:
        df_temp_nominal['group'] = list(encoded__labels)
        df_temp_nominal.loc[df_temp_nominal['group'] != group[0], 'group'] = 'Rest'
        df_temp_nominal['group'] = df_temp_nominal['group'].astype('category')

        data_temp = [df_nominal.loc[ids, :] for ids in df_temp_nominal.groupby('group').groups.values()]

        # Compare in a one vs rest fashion
        pval_list = {'Chi-square': []}
        effect_size = []
        effect_size_categorised = []
        mode_group_0 = []
        mode_group_rest = []
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
            cont_table = pd.crosstab(df_melted_tmp['value'], df_melted_tmp['variable'])

            # perform chi-square test of independence
            chi2, p_temp, dof, expected = chi2_contingency(cont_table)

            pval_list['Chi-square'].append((p_temp, nominal_attr))

            # Calculate odds ratio
            odds_ratio = (cont_table[0][0] * cont_table[1][1]) / (cont_table[1][0] * cont_table[0][1])
            effect_size.append((nominal_attr, odds_ratio))
            size_def_tmp = lookup_oddsratio(odds_ratio=odds_ratio)
            effect_size_categorised.append((nominal_attr, size_def_tmp))
            mode_group_0.append((nominal_attr, df_melted_tmp.loc[
                df_melted_tmp['variable'] == 0, 'value'].mode()[0]))
            mode_group_rest.append((nominal_attr, df_melted_tmp.loc[
                df_melted_tmp['variable'] == 1, 'value'].mode()[0]))

        # Save Chi-square test results to dataframe
        d_swap = {v: k for k, v in dict(pval_list['Chi-square']).items()}
        df_chisquare = pd.concat([
            df_chisquare, pd.DataFrame(d_swap, index=["{} vs {}".format(group[0], group[1])]).T], axis=1)

        # Save effect size
        df_effect_size = pd.concat(
            [df_effect_size,
             pd.DataFrame(dict(effect_size), index=["{} vs {}".format(group[0], group[1])]).T], axis=1)
        df_effect_size_categorised = pd.concat(
            [df_effect_size_categorised,
             pd.DataFrame(dict(effect_size_categorised), index=["{} vs {}".format(group[0], group[1])]).T], axis=1)

        df_mode_group_0 = pd.concat([
            df_mode_group_0, pd.DataFrame(dict(mode_group_0), index=[group[0]]).T], axis=1)
        df_mode_group_rest = pd.concat([
            df_mode_group_rest, pd.DataFrame(dict(mode_group_rest), index=["Rest ({})".format(group[0])]).T], axis=1)

        # Using the padj values -> find out in which class the attribute is significantly different
        test_key = "Chi-square"
        if len(pval_list[test_key]) > 0:
            # loop through padj list and read out df containing padj values and attribute
            for pvals, attr in pval_list[test_key]:
                # Save padj values for each pairwise comparison in dataframe
                df_all_nominal_attributes.loc[attr, [group]] = pvals

        if len(pval_list[test_key]) > 0:
            for attr, effect in effect_size:
                # Save padj values for each pairwise comparision in dataframe
                df_all_nominal_attributes_effect_size.loc[attr, [group]] = effect

    df_chisquare.to_excel(os.path.join(save_folder_tmp, 'Stats_Nominal.xlsx'), sheet_name='Chi-square')
    with pd.ExcelWriter(os.path.join(save_folder_tmp, 'Stats_Nominal.xlsx'),
                        mode="a", engine="openpyxl") as writer:
        df_effect_size.to_excel(writer, sheet_name="Cliff's Delta", index=True)
        df_effect_size_categorised.to_excel(writer, sheet_name="Cliff's Delta categorised", index=True)
        df_mode_group_0.to_excel(writer, sheet_name='Mode group Endotypes', index=True)
        df_mode_group_rest.to_excel(writer, sheet_name='Mode group Rest', index=True)

    return df_all_nominal_attributes, df_all_nominal_attributes_effect_size


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


# Function to sign padj values
def sign_padj(row):
    if 'greater' in row['Effect']:
        return -row['log10(padj)']
    elif 'less' in row['Effect']:
        return row['log10(padj)']
    else:
        return 0


def shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
    '''
    Function to offset the "center" of a colormap. Useful for
    data with a negative min and positive max and you want the
    middle of the colormap's dynamic range to be at zero.

    Input
    -----
      cmap : The matplotlib colormap to be altered
      start : Offset from lowest point in the colormap's range.
          Defaults to 0.0 (no lower offset). Should be between
          0.0 and `midpoint`.
      midpoint : The new center of the colormap. Defaults to
          0.5 (no shift). Should be between 0.0 and 1.0. In
          general, this should be  1 - vmax / (vmax + abs(vmin))
          For example if your data range from -15.0 to +5.0 and
          you want the center of the colormap at 0.0, `midpoint`
          should be set to  1 - 5/(5 + 15)) or 0.75
      stop : Offset from highest point in the colormap's range.
          Defaults to 1.0 (no upper offset). Should be between
          `midpoint` and 1.0.
    '''
    cdict = {
        'red': [],
        'green': [],
        'blue': [],
        'alpha': []
    }

    # regular index to compute the colors
    reg_index = np.linspace(start, stop, 257)

    # shifted index to match the data
    shift_index = np.hstack([
        np.linspace(0.0, midpoint, 128, endpoint=False),
        np.linspace(midpoint, 1.0, 129, endpoint=True)
    ])

    for ri, si in zip(reg_index, shift_index):
        r, g, b, a = cmap(ri)

        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))

    newcmap = matplotlib.colors.LinearSegmentedColormap(name, cdict)
    matplotlib.colormaps.register(cmap=newcmap)

    return newcmap


def main(input_dir, input_meta_dir, save_folder):
    adata = sc.read(os.path.join(
        input_dir,
        'Lesion_RNAseq_20210720_patient_meta_data_v04__CH_20220209_TherapyResponse_v04_v220715_TherapyResponseCorrected__Endotypes_230605.h5'))
    # # was wronly set to 0 -> delete zeros
    # adata.obs.loc[adata.obs.index.isin(['MUC3671', 'MUC20449']), 'Hty_exanth'] = np.nan
    df_corrected_metadata = pd.read_excel(os.path.join(
        '/Volumes', 'CH__data', 'Projects', 'Eyerich_AG_projects', 'BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer',
        'raw_data', 'clinical_data', '20240730_patient_meta_data_final_PS.xlsx'))
    # rename Sex.x, Helmholtz.identifyer', Diag_diag.x
    df_corrected_metadata = df_corrected_metadata.rename(
        columns={'Helmholtz.identifyer': 'Helmholtz_identifyer', 'Diag_diag.x': 'Diag_diag_x', 'Sex.x': 'Sex_x'})
    df_corrected_metadata.index = df_corrected_metadata['Helmholtz_identifyer']

    # common columns
    common_cols = np.intersect1d(adata.obs.columns, df_corrected_metadata.columns)
    assert np.all(df_corrected_metadata.loc[adata.obs.index, common_cols].index == adata.obs.index)
    # TODO might have to adjust some variables -> Peter
    adata.obs.loc[:, common_cols] = df_corrected_metadata.loc[adata.obs.index, common_cols]
    adata.obs = adata.obs.replace({"Pattern": {"undefined": "UD"}})

    obs_name = 'Endotypes'  # 'Molecular Subtype res0.9'

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
    dict_all_continuous_attributes, dict_all_continuous_attributes_effect_size = perform_test_continuous_variables(
        df=df, continuous_cols=continuous_cols, encoded__labels=encoded__labels, save_folder=save_folder)

    # 2. Ordinal (Dichotomus) categories
    df_all_ordinal_attributes, df_all_ordinal_attributes_effect_size = perform_test_ordinal_category(
        df_scaled=df, ordinal_cols=ordinal_cols, encoded__labels=encoded__labels, save_folder=save_folder)

    # 3. Nominal categories
    df_all_nominal_attributes, df_all_nominal_attributes_effect_size = perform_test_nominal_category(
        df=df, nominal_cols=nominal_cols, encoded__labels=encoded__labels, save_folder=save_folder)

    # Plot summary plot
    df_all_continuous_attributes = pd.DataFrame()
    for key_name in dict_all_continuous_attributes:
        df_all_continuous_attributes = pd.concat(
            [df_all_continuous_attributes, dict_all_continuous_attributes[key_name]])

        plot_summary(df=dict_all_continuous_attributes[key_name], save_folder=save_folder,
                     key='continuous_{}'.format(key_name))

    df_all_continuous_attributes_effect_size = pd.DataFrame()
    for key_name in dict_all_continuous_attributes_effect_size:
        df_all_continuous_attributes_effect_size = pd.concat(
            [df_all_continuous_attributes_effect_size, dict_all_continuous_attributes_effect_size[key_name]])
    # Drop na rows
    df_all_continuous_attributes_effect_size = df_all_continuous_attributes_effect_size.dropna(
        axis='index', how='all').astype(float)
    # Categorise effect size (Cohen's D) of continuous attributes
    df_all_continuous_attributes_effect_size_categorised = pd.DataFrame(
        index=df_all_continuous_attributes_effect_size.index, columns=df_all_continuous_attributes_effect_size.columns)
    for idx in df_all_continuous_attributes_effect_size_categorised.index:
        for col in df_all_continuous_attributes_effect_size_categorised.columns:
            df_all_continuous_attributes_effect_size_categorised.loc[idx][col] = categorize_cohens_d(
                cohens_d=df_all_continuous_attributes_effect_size.loc[idx][col])

    plot_summary(df=df_all_ordinal_attributes, save_folder=save_folder, key='ordinal')
    plot_summary(df=df_all_nominal_attributes, save_folder=save_folder, key='nominal')

    # Categorise effect size (Cliff's delta) of ordinal attributes
    df_all_ordinal_attributes_effect_size_categorised = pd.DataFrame(
        index=df_all_ordinal_attributes_effect_size.index, columns=df_all_ordinal_attributes_effect_size.columns)
    for idx in df_all_ordinal_attributes_effect_size_categorised.index:
        for col in df_all_ordinal_attributes_effect_size_categorised.columns:
            df_all_ordinal_attributes_effect_size_categorised.loc[idx][col] = lookup_size(
                delta=df_all_ordinal_attributes_effect_size.loc[idx][col])

    # Categorise effect size (Odds Ratio) of nominal attributes
    df_all_nominal_attributes_effect_size_categorised = pd.DataFrame(
        index=df_all_nominal_attributes_effect_size.index, columns=df_all_nominal_attributes_effect_size.columns)
    for idx in df_all_nominal_attributes_effect_size_categorised.index:
        for col in df_all_nominal_attributes_effect_size_categorised.columns:
            df_all_nominal_attributes_effect_size_categorised.loc[idx][col] = lookup_oddsratio(
                odds_ratio=df_all_nominal_attributes_effect_size.loc[idx][col])

    # concatenate categorised effect sizes
    df_effectsize_categorised = pd.concat([
        df_all_continuous_attributes_effect_size_categorised, df_all_ordinal_attributes_effect_size_categorised,
        df_all_nominal_attributes_effect_size_categorised])

    # concatenate dfs
    df_merged = pd.concat([df_all_continuous_attributes, df_all_ordinal_attributes, df_all_nominal_attributes])
    # reorder columns like sorted in dendrogram
    adata.obs['Endotypes'] = adata.obs['Endotypes'].cat.reorder_categories(
        ['E5', 'E6', 'E7', 'E8', 'E9', 'E10', 'E11', 'E12', 'E13', 'E1', 'E2', 'E3', 'E4'])
    sorted_colnames = [(str(s), 'Rest') for s in list(adata.obs[obs_name].cat.categories)]
    df_merged = df_merged.loc[:, sorted_colnames]

    # Drop na rows
    df_merged = df_merged.dropna(axis='index', how='all').astype(float)

    # All p-values at once
    p_values = df_merged.values.flatten()
    _, pvals_fdr, _, _ = multipletests(p_values, method='fdr_bh')
    # Reshape the adjusted p-values back to the original DataFrame shape
    df_padj = pd.DataFrame(pvals_fdr, index=pd.MultiIndex.from_product([df_merged.index, df_merged.columns]),
                           columns=['Adjusted_p_value']).reset_index()
    # Convert MultiIndex to DataFrame
    df_padj = df_padj.pivot(index='level_0', columns='level_1', values='Adjusted_p_value')
    # Reorder columns to match original DataFrame
    df_padj = df_padj[df_merged.columns]
    df_padj = df_padj.loc[df_merged.index, :]

    df_merged_effect_size = pd.concat([
        df_all_continuous_attributes_effect_size, df_all_ordinal_attributes_effect_size,
        df_all_nominal_attributes_effect_size])
    df_merged_effect_size = df_merged_effect_size.loc[:, sorted_colnames]
    # Drop na rows
    df_merged_effect_size = df_merged_effect_size.dropna(axis='index', how='all').astype(float)

    # Save all in one excel sheet
    df_merged.to_excel(os.path.join(save_folder, 'Stats_all_enriched_clinical_attributes_Endotypes.xlsx'),
                       sheet_name='p-value', index=True)
    with pd.ExcelWriter(os.path.join(save_folder, 'Stats_all_enriched_clinical_attributes_Endotypes.xlsx'),
                        mode="a", engine="openpyxl") as writer:
        df_padj.to_excel(writer, sheet_name='padj. value', index=True)
        df_merged_effect_size.to_excel(writer, sheet_name="effect size", index=True)
        df_effectsize_categorised.to_excel(writer, sheet_name="effect size categorised", index=True)

    # save data
    df_merged_effect_size.to_excel(os.path.join(
        save_folder, 'Effectsize_all_enriched_clinical_attributes_Endotypes.xlsx'))

    # Waterfallplot
    # Unpivot dataframes to long format
    df_padj_long = df_padj.reset_index().melt(id_vars='index', var_name='Comparison', value_name='padj')
    df_effect_long = df_effectsize_categorised.reset_index().melt(
        id_vars='index', var_name='Comparison', value_name='Effect')

    # Combine the DataFrames
    df_combined = pd.merge(df_padj_long, df_effect_long, on=['index', 'Comparison'])
    # Add log10 padj value
    df_combined['log10(padj)'] = df_combined['padj'].apply(np.log10, axis=1)

    # Apply the function
    df_combined['signed_padj'] = df_combined.apply(sign_padj, axis=1)

    # Mapping effect categories to colors
    effect_color_map = {
        'less, Negligible effect': '#d3d3d3',   # light gray
        'less, Small effect': '#87ceeb',        # light blue
        'less, Medium effect': '#4682b4',       # steel blue
        'less, Large effect': '#000080',        # navy
        'greater, Negligible effect': '#ffcccb',# light red
        'greater, Small effect': '#ff7f7f',     # salmon
        'greater, Medium effect': '#ff4500',    # orange red
        'greater, Large effect': '#8b0000',     # dark red
        'equal, Negligible effect': '#a9a9a9',                     # dark gray
        'not defined': '#808080'                # gray
    }

    # Plotting per comparison
    for endotype in df_combined['Comparison'].unique():
        mask = df_combined['Comparison'].isin([endotype])
        df_tmp = df_combined.loc[mask, :]
        df_tmp['Effect'] = df_tmp['Effect'].astype(str)

        df_tmp.loc[:, 'color'] = df_tmp.loc[:, 'Effect'].map(effect_color_map)

        # Sort by signed padj values for the waterfall effect
        df_tmp = df_tmp.sort_values(by='signed_padj').reset_index(drop=True)

        # Drop all attributes with padj value > 0.2
        df_tmp = df_tmp[df_tmp['padj'] < 0.2]

        # Filter the effect_colors dictionary to include only the effects present in the data
        used_effect_colors = {
            effect: color for effect, color in effect_color_map.items() if effect in df_tmp['Effect'].unique()}

        # Sort the dataframe by (Effect), then by signed_padj (negative values in descending order)
        # df_tmp = df_tmp.sort_values(by=['Effect', 'signed_padj']).reset_index(drop=True)
        df_tmp = df_tmp.sort_values(by=['signed_padj'], ascending=False).reset_index(drop=True)
        # df_tmp = df_tmp.sort_values(
        #     by=['Effect', 'signed_padj'],
        #     key=lambda col: custom_sort(col) if col.name == 'signed_padj' else col
        # ).reset_index(drop=True)

        fig, ax = plt.subplots(figsize=(12, 4))
        # Plot bars
        ax.bar(df_tmp.index, df_tmp['signed_padj'], color=df_tmp['color'])

        # Draw horizontal lines between effect categories
        # prev_effect = None
        # for i, effect in enumerate(df_tmp['Effect']):
        #     if effect != prev_effect:
        #         if prev_effect is not None:
        #             ax.axvline(x=i-0.5, color='black', linestyle='-')
        #         prev_effect = effect

        # Add legend
        handles = [plt.Rectangle((0, 0), 1, 1, color=color) for color in used_effect_colors.values()]
        labels = used_effect_colors.keys()
        ax.legend(handles, labels, title="Effect Categories", loc='upper left', bbox_to_anchor=(1, 1))

        # Set x-axis labels to show the original categories
        ax.set_xticks(df_tmp.index)
        ax.set_xticklabels(df_tmp['index'], rotation=90, fontsize=8)
        # Adding labels
        ax.set_xlabel('Clinical attributes')
        ax.set_ylabel('Signed log10(p-adj. values)')
        sns.despine(fig=fig)

        fig.savefig(
            os.path.join(save_folder, '{}_Waterfallplot_clinicalattributes_padj_cutoff_02.pdf'.format(
                "_vs_".join(endotype))),
            bbox_inches='tight')
        plt.close(fig=fig)

    # # Melt the DataFrame
    # df_combined_tmp = df_combined[['index', 'Comparison', 'Effect', 'signed_padj']].copy()
    # df_melted = pd.melt(df_combined_tmp, id_vars=['Comparison', 'index'], value_vars=['signed_padj'])
    #
    # df_pivoted_padj = df_melted.pivot_table(columns=['Comparison'], values='value', index=['index'])
    #
    # circle_sizes = {'Large effect': 50, 'Medium effect': 25, 'Small effect': 10,
    #                 'Negligible effect': 1, 'not defined': 0}
    # df_sizes = df_effectsize_categorised.copy()
    # for col in df_sizes.columns:
    #     df_sizes[col] = df_sizes[col].str.split(', ').str[1]
    #     # transform effect size in size values
    #     df_sizes[col] = df_sizes[col].replace(circle_sizes, regex=True)
    #
    # # Extract values and labels
    # values = df_pivoted_padj.values
    # xlabels = df_pivoted_padj.columns
    # ylabels = df_pivoted_padj.index
    #
    # N, M = values.shape  # Columns, Rows
    #
    # # Values in heatmap
    # x, y = np.meshgrid(np.arange(M), np.arange(N))
    # # Size of circles
    # sizes = df_sizes.values
    # # s = np.random.randint(0, 180, size=(N, M))
    # # Colors in heatmap
    # c = values
    #
    # R = sizes/sizes.max()/2
    #
    # pval_cut = 0.05
    #
    # # Create colormap
    # divnorm = matplotlib.colors.TwoSlopeNorm(vmin=values.min(), vcenter=0, vmax=values.max())
    #
    # cmap_reds = plt.get_cmap('Reds')
    # cmap_blues = plt.get_cmap('Blues')
    # num_colors = 80
    # colors = [cmap_blues(i / num_colors) for i in range(8, num_colors + 4)][::-1] + [
    #     cmap_reds(i / num_colors) for i in range(8, num_colors + 4)]
    # cmap = LinearSegmentedColormap.from_list('', colors, len(colors))
    # cmap.set_bad("lightgrey")
    #
    # fig, ax = plt.subplots(figsize=(6, 6))  # (w, h)
    # circles = [plt.Circle((j, i), radius=r) for r, j, i in zip(R.flat, x.flat, y.flat)]
    # # Colorcode circles matplotlib.colors.CenteredNorm()
    # col = PatchCollection(circles, array=c.flatten(), norm=divnorm, cmap=matplotlib.colormaps.get_cmap('RdBu_r'))
    # ax.add_collection(col)
    #
    # ax.set(xticks=np.arange(M), yticks=np.arange(N),
    #        xticklabels=xlabels, yticklabels=ylabels)
    # ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    # ax.set_xticks(np.arange(M+1)-0.5, minor=True)
    # ax.set_yticks(np.arange(N+1)-0.5, minor=True)
    # ax.grid(which='minor')
    #
    # # Add colorbar
    # # fig.colorbar(col)
    # axcolor = fig.add_axes([0.2, 0.5, 0.65, 0.02])  # (left, bottom, width, height)
    # cbar = plt.colorbar(
    #     mappable=ax.collections[0], cax=axcolor, orientation="horizontal",
    #     ticks=[np.amin(values)] + [-10] + [
    #         np.log10(pval_cut)] + [0] + [-np.log10(pval_cut)] + [12] + [np.amax(values)])
    # # cbar.ax.get_xaxis().labelpad = 15
    # cbar.ax.set_xlabel('-log10(p-values)', rotation=0, fontsize=12)
    #
    # # With seaborn
    # df_combined['Effect size plot'] = df_combined['Effect']
    # df_combined['Effect size plot'] = df_combined['Effect'].str.split(', ').str[1]
    # # transform effect size in size values
    # df_combined['Effect size plot'] = df_combined['Effect size plot'].replace(circle_sizes, regex=True)
    # df_combined['Comparison'] = df_combined['Comparison'].apply(lambda x: " vs ".join(x))
    #
    # fig, ax = plt.subplots(figsize=(8, 12))
    # sns.scatterplot(x="Comparison", y="index", hue="signed_padj", size="Effect size plot",
    #                 sizes=(1, 250), data=df_combined, palette='RdBu_r', linewidth=0, legend=True, ax=ax)
    # ax.invert_yaxis()
    # ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    # sm = plt.cm.ScalarMappable(cmap="RdBu_r", norm=divnorm)
    # sm.set_array([])
    # # Remove the legend and add a colorbar
    # ax.get_legend().remove()
    # cbar = ax.figure.colorbar(sm)
    #
    # # cbar.ax.get_xaxis().labelpad = 15
    # cbar.ax.set_ylabel('-log10(p-values)', rotation=90, fontsize=12)
    #
    #
    # keep = ['Hist_ID_quant', 'Hist_Nr_Keratoses', 'Hist_Granu', 'Hist_Para_quant', 'Hist_Aka_quant',
    #         'Hist_Neutro_quant', 'Hist_Lympho', 'Morph_pustules', 'Morph_scaling', 'Hist_LCV', 'Com_muc_2.0']
    # attr_names = ['Histology - Interface dermatitis (quantitative)', 'Histology - Dyskeratoses (quantitative)',
    #               'Histology - Abnormal Str. granulosum', 'Histology - Parakeratosis (quantitative)',
    #               'Histology - Akanthosis (quantitative)', 'Histology - Neutrophils (quantitative)',
    #               'Histology - Dermal lymphocytes (quantitative)', 'Morphology - Pustules', 'Morphology - Scaling',
    #               'History - Leukocytoclastic vasculitis', 'Comorbidity - Mucosa involved']
    # rename_attr = dict(zip(keep, attr_names))
    #
    # # Select attributes - see email BRAIN - Figure 3D Heatmap 1 vs Rest clinical attributes from 14.09.2023
    # df_combined_preselected = df_combined.loc[df_combined['index'].isin(keep), :]
    # df_combined_preselected['index'] = df_combined_preselected['index'].replace(rename_attr, regex=True)
    #
    # pws = df_combined_preselected['Effect size plot'].unique()
    # fig, ax = plt.subplots(figsize=(8, 12))
    # ax.scatter(df_combined_preselected.Comparison, df_combined_preselected['index'],
    #            c=df_combined_preselected.signed_padj, s=(df_combined_preselected['Effect size plot']),
    #            cmap="RdBu_r")
    # ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    # # Add colorbar
    # axcolor = fig.add_axes([0.2, 0.02, 0.65, 0.02])  # (left, bottom, width, height)
    # cbar = plt.colorbar(
    #     mappable=ax.collections[0], cax=axcolor, orientation="horizontal")
    # cbar.ax.set_xlabel('-log10(p-values)', rotation=90, fontsize=12)
    # # make a legend:
    # for pw in pws:
    #     ax.scatter([], [], s=pw, c="k", label=str(pw))
    # h, l = ax.get_legend_handles_labels()
    # fig.legend(h[1:], l[1:], labelspacing=1.2, title="Effect size", borderpad=1,
    #            frameon=True, framealpha=0.6, edgecolor="k", facecolor="w")
    # sns.despine(fig=fig, ax=ax)


    pval_cut = 0.05
    circle_sizes = {'Large effect': 50, 'Medium effect': 25, 'Small effect': 10,
                    'Negligible effect': 1, 'not defined': 0}
    keep = ['Hist_ID_quant', 'Hist_Nr_Keratoses', 'Hist_Granu', 'Hist_Para_quant', 'Hist_Aka_quant',
            'Hist_Neutro_quant', 'Hist_Lympho', 'Morph_pustules', 'Morph_scaling', 'Hist_LCV', 'Com_muc_2.0']
    attr_names = ['Histology - Interface dermatitis (quantitative)', 'Histology - Dyskeratoses (quantitative)',
                  'Histology - Abnormal Str. granulosum', 'Histology - Parakeratosis (quantitative)',
                  'Histology - Akanthosis (quantitative)', 'Histology - Neutrophils (quantitative)',
                  'Histology - Dermal lymphocytes (quantitative)', 'Morphology - Pustules', 'Morphology - Scaling',
                  'History - Leukocytoclastic vasculitis', 'Comorbidity - Mucosa involved']
    rename_attr = dict(zip(keep, attr_names))

    # Select attributes - see email BRAIN - Figure 3D Heatmap 1 vs Rest clinical attributes from 14.09.2023
    df_combined_preselected = df_combined.loc[df_combined['index'].isin(keep), :]
    df_combined_preselected['index'] = df_combined_preselected['index'].replace(rename_attr, regex=True)

    df_combined_preselected['Effect size plot'] = df_combined_preselected['Effect']
    df_combined_preselected['Effect size plot'] = df_combined_preselected['Effect'].str.split(', ').str[1]
    # transform effect size in size values
    df_combined_preselected['Effect size plot'] = df_combined_preselected['Effect size plot'].replace(circle_sizes, regex=True)
    df_combined_preselected['Comparison'] = df_combined_preselected['Comparison'].apply(lambda x: " vs ".join(x))
    df_combined_preselected.to_excel(os.path.join(save_folder, 'Selected_Attributes_plot_information.xlsx'))

    pws = df_combined_preselected['Effect size plot'].unique()
    reverse_circle_sizes = {y: x for x, y in circle_sizes.items()}
    df_colors = pd.DataFrame(
        data=np.asarray([["#ff6db6"], ["#490092"],
                         ["#b66dff"], ["#000000"], ["#920000"], ["#E69F00"], ["#D55E00"], ["#8B4513"],
                         ["#999999"], ["#006ddb"], ["#b6dbff"], ["#004949"], ["#009292"]]).T,
        columns=['E1', 'E2', 'E3', 'E4', 'E5', 'E6', 'E7', 'E8', 'E9', 'E10', 'E11', 'E12', 'E13'])
    vmin = -np.ceil(np.abs(df_combined_preselected.signed_padj.min()))
    vmax = np.ceil(df_combined_preselected.signed_padj.max())
    norm = matplotlib.colors.TwoSlopeNorm(vcenter=0, vmin=vmin, vmax=vmax)

    # Create colormaps for blue and red
    cmap_blues = plt.get_cmap('Blues')
    cmap_reds = plt.get_cmap('Reds')
    num_colors = 40  # defines color intensity in formula below

    colors = [cmap_blues(i / num_colors) for i in range(4, 24 + 4)][::-1] + [
        'lightgrey'] * 9 + [cmap_reds(i / num_colors) for i in range(4, 28 + 4)]
    custom_cmap = LinearSegmentedColormap.from_list('', colors, len(colors))

    fig = plt.figure(figsize=(6, 4))  # (w, h)
    # Add heatmap showing the colors of the clusters - TODO split into three columns
    ax_colors = fig.add_axes([0.2, 0.19, 0.65, 0.02])  # (left, bottom, width, height)
    cmap_colors = LinearSegmentedColormap.from_list('', df_colors.iloc[0, :].values, len(df_colors.iloc[0, :].values))
    gradient = np.linspace(0, 1, len(df_colors.columns))
    gradient = np.vstack((gradient, gradient))
    ax_colors.imshow(gradient, aspect='auto', cmap=cmap_colors)
    ax_colors.get_xticklabels([])
    ax_colors.get_yticklabels([])
    ax_colors.get_xaxis().set_visible(False)
    ax_colors.get_yaxis().set_visible(False)

    # Compute and plot second dendrogram. (left)
    ax2 = fig.add_axes([0.16, 0.2, 0.03, 0.72])  # (left, bottom, width, height) fig.add_axes([0.095, 0.1, 0.1, 0.7])
    ax2.set_yticks([])
    ax2.set_xticks([])
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.spines['bottom'].set_visible(False)
    ax2.spines['left'].set_visible(False)

    # Split heatmap into Pattern 2a-like, 3-like, and 1-like columns
    # Pattern 2a-like
    df_pattern_2a = df_combined_preselected.loc[
                    df_combined_preselected['Comparison'].isin(
                        ['E5 vs Rest',  'E6 vs Rest',  'E7 vs Rest', 'E8 vs Rest', 'E9 vs Rest', 'E10 vs Rest']), :]
    axmatrix = fig.add_axes([0.2, 0.2, 0.28, 0.72])  # (left, bottom, width, height)
    axmatrix.scatter(df_pattern_2a.Comparison, df_pattern_2a['index'],
                     c=df_pattern_2a.signed_padj, s=(df_pattern_2a['Effect size plot']),
                     cmap=custom_cmap, norm=norm)
    axmatrix.invert_yaxis()
    axmatrix.set_yticklabels(axmatrix.get_yticklabels(), rotation=0)
    axmatrix.set_xticklabels(axmatrix.get_xticklabels(), rotation=90)
    axmatrix.set_yticklabels(axmatrix.get_yticklabels(), rotation=0)
    # axmatrix.set_yticklabels([])
    # axmatrix.set_yticks([])
    axmatrix.spines['top'].set_visible(True)
    axmatrix.spines['right'].set_visible(True)
    axmatrix.spines['bottom'].set_visible(True)
    axmatrix.spines['left'].set_visible(True)

    # Pattern 3-like
    df_pattern_3 = df_combined_preselected.loc[
                    df_combined_preselected['Comparison'].isin(
                        ['E11 vs Rest', 'E12 vs Rest', 'E13 vs Rest']), :]
    axmatrix = fig.add_axes([0.50, 0.2, 0.14, 0.72])  # (left, bottom, width, height)
    axmatrix.scatter(df_pattern_3.Comparison, df_pattern_3['index'],
                     c=df_pattern_3.signed_padj, s=(df_pattern_3['Effect size plot']),
                     cmap=custom_cmap, norm=norm)
    axmatrix.invert_yaxis()
    axmatrix.set_yticklabels(axmatrix.get_yticklabels(), rotation=0)
    axmatrix.set_xticklabels(axmatrix.get_xticklabels(), rotation=90)
    axmatrix.set_yticklabels([])
    axmatrix.set_yticks([])
    axmatrix.spines['top'].set_visible(True)
    axmatrix.spines['right'].set_visible(True)
    axmatrix.spines['bottom'].set_visible(True)
    axmatrix.spines['left'].set_visible(True)

    # Pattern 1-like
    df_pattern_1 = df_combined_preselected.loc[
                    df_combined_preselected['Comparison'].isin(
                        ['E1 vs Rest',  'E2 vs Rest',  'E3 vs Rest', 'E4 vs Rest']), :]
    axmatrix = fig.add_axes([0.66, 0.2, 0.19, 0.72])  # (left, bottom, width, height)
    axmatrix.scatter(df_pattern_1.Comparison, df_pattern_1['index'],
                     c=df_pattern_1.signed_padj, s=(df_pattern_1['Effect size plot']),
                     cmap=custom_cmap, norm=norm)
    axmatrix.invert_yaxis()
    axmatrix.set_yticks([])
    # axmatrix.set_yticklabels(axmatrix.get_yticklabels(), rotation=0)
    axmatrix.set_xticklabels(axmatrix.get_xticklabels(), rotation=90)
    # axmatrix.yaxis.set_label_position("right")
    # axmatrix.yaxis.tick_right()
    axmatrix.spines['top'].set_visible(True)
    axmatrix.spines['right'].set_visible(True)
    axmatrix.spines['bottom'].set_visible(True)
    axmatrix.spines['left'].set_visible(True)

    # Add colorbar
    axcolor = fig.add_axes([0.2, 0.02, 0.65, 0.02])  # (left, bottom, width, height)
    cbar = plt.colorbar(
        mappable=axmatrix.collections[0], cax=axcolor, orientation="horizontal",
        ticks=[vmin, -4.5, -3] + [np.log10(pval_cut)] + [0] + [-np.log10(pval_cut)] + [5, 10] + [vmax])
    cbar.ax.set_xlabel('signed -log10(p-adj. values)', rotation=0, fontsize=12)

    # make a legend:
    for pw in pws:
        axmatrix.scatter([], [], s=pw, c="k", label=reverse_circle_sizes[pw])
    h, l = axmatrix.get_legend_handles_labels()
    fig.legend(h, l, labelspacing=1.2, title="Effect size", borderpad=1, loc='center left', bbox_to_anchor=(0.8, 0.5),
               frameon=True, framealpha=0.6, edgecolor="k", facecolor="w")

    plt.savefig(os.path.join(save_folder, "Selected_attr_Heatmap_pval_{}_ok.pdf".format(pval_cut)),
                bbox_inches='tight')
    plt.close()


    df_combined['Effect size plot'] = df_combined['Effect']
    df_combined['Effect size plot'] = df_combined['Effect'].str.split(', ').str[1]
    # transform effect size in size values
    df_combined['Effect size plot'] = df_combined['Effect size plot'].replace(circle_sizes, regex=True)
    df_combined['Comparison'] = df_combined['Comparison'].apply(lambda x: " vs ".join(x))
    pws = df_combined['Effect size plot'].unique()
    reverse_circle_sizes = {y: x for x, y in circle_sizes.items()}
    df_colors = pd.DataFrame(
        data=np.asarray([["#ff6db6"], ["#490092"],
                         ["#b66dff"], ["#000000"], ["#920000"], ["#E69F00"], ["#D55E00"], ["#8B4513"],
                         ["#999999"], ["#006ddb"], ["#b6dbff"], ["#004949"], ["#009292"]]).T,
        columns=['E1', 'E2', 'E3', 'E4', 'E5', 'E6', 'E7', 'E8', 'E9', 'E10', 'E11', 'E12', 'E13'])
    vmin = -np.ceil(np.abs(df_combined.signed_padj.min()))
    vmax = np.ceil(df_combined.signed_padj.max())
    norm = matplotlib.colors.TwoSlopeNorm(vcenter=0, vmin=vmin, vmax=vmax)

    # Create colormaps for blue and red
    cmap_blues = plt.get_cmap('Blues')
    cmap_reds = plt.get_cmap('Reds')
    num_colors = 40  # defines color intensity in formula below

    colors = [cmap_blues(i / num_colors) for i in range(4, 24 + 4)][::-1] + [
        'lightgrey'] * 9 + [cmap_reds(i / num_colors) for i in range(4, 28 + 4)]
    custom_cmap = LinearSegmentedColormap.from_list('', colors, len(colors))

    fig = plt.figure(figsize=(6, 12))  # (w, h)
    # Add heatmap showing the colors of the clusters - TODO split into three columns
    ax_colors = fig.add_axes([0.2, 0.19, 0.65, 0.02])  # (left, bottom, width, height)
    cmap_colors = LinearSegmentedColormap.from_list('', df_colors.iloc[0, :].values, len(df_colors.iloc[0, :].values))
    gradient = np.linspace(0, 1, len(df_colors.columns))
    gradient = np.vstack((gradient, gradient))
    ax_colors.imshow(gradient, aspect='auto', cmap=cmap_colors)
    ax_colors.get_xticklabels([])
    ax_colors.get_yticklabels([])
    ax_colors.get_xaxis().set_visible(False)
    ax_colors.get_yaxis().set_visible(False)

    # Compute and plot second dendrogram. (left)
    ax2 = fig.add_axes([0.16, 0.2, 0.03, 0.72])  # (left, bottom, width, height) fig.add_axes([0.095, 0.1, 0.1, 0.7])
    ax2.set_yticks([])
    ax2.set_xticks([])
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.spines['bottom'].set_visible(False)
    ax2.spines['left'].set_visible(False)

    # Split heatmap into Pattern 2a-like, 3-like, and 1-like columns
    # Pattern 2a-like
    df_pattern_2a = df_combined.loc[
                    df_combined['Comparison'].isin(
                        ['E5 vs Rest',  'E6 vs Rest',  'E7 vs Rest', 'E8 vs Rest', 'E9 vs Rest', 'E10 vs Rest']), :]
    axmatrix = fig.add_axes([0.2, 0.2, 0.28, 0.72])  # (left, bottom, width, height)
    axmatrix.scatter(df_pattern_2a.Comparison, df_pattern_2a['index'],
                     c=df_pattern_2a.signed_padj, s=(df_pattern_2a['Effect size plot']),
                     cmap=custom_cmap, norm=norm)
    axmatrix.invert_yaxis()
    axmatrix.set_yticklabels(axmatrix.get_yticklabels(), rotation=0)
    axmatrix.set_xticklabels(axmatrix.get_xticklabels(), rotation=90)
    axmatrix.set_yticklabels(axmatrix.get_yticklabels(), rotation=0)
    # axmatrix.set_yticklabels([])
    # axmatrix.set_yticks([])
    axmatrix.spines['top'].set_visible(True)
    axmatrix.spines['right'].set_visible(True)
    axmatrix.spines['bottom'].set_visible(True)
    axmatrix.spines['left'].set_visible(True)

    # Pattern 3-like
    df_pattern_3 = df_combined.loc[
                    df_combined['Comparison'].isin(
                        ['E11 vs Rest', 'E12 vs Rest', 'E13 vs Rest']), :]
    axmatrix = fig.add_axes([0.50, 0.2, 0.14, 0.72])  # (left, bottom, width, height)
    axmatrix.scatter(df_pattern_3.Comparison, df_pattern_3['index'],
                     c=df_pattern_3.signed_padj, s=(df_pattern_3['Effect size plot']),
                     cmap=custom_cmap, norm=norm)
    axmatrix.invert_yaxis()
    axmatrix.set_yticklabels(axmatrix.get_yticklabels(), rotation=0)
    axmatrix.set_xticklabels(axmatrix.get_xticklabels(), rotation=90)
    axmatrix.set_yticklabels([])
    axmatrix.set_yticks([])
    axmatrix.spines['top'].set_visible(True)
    axmatrix.spines['right'].set_visible(True)
    axmatrix.spines['bottom'].set_visible(True)
    axmatrix.spines['left'].set_visible(True)

    # Pattern 1-like
    df_pattern_1 = df_combined.loc[
                    df_combined['Comparison'].isin(
                        ['E1 vs Rest',  'E2 vs Rest',  'E3 vs Rest', 'E4 vs Rest']), :]
    axmatrix = fig.add_axes([0.66, 0.2, 0.19, 0.72])  # (left, bottom, width, height)
    axmatrix.scatter(df_pattern_1.Comparison, df_pattern_1['index'],
                     c=df_pattern_1.signed_padj, s=(df_pattern_1['Effect size plot']),
                     cmap=custom_cmap, norm=norm)
    axmatrix.invert_yaxis()
    axmatrix.set_yticks([])
    # axmatrix.set_yticklabels(axmatrix.get_yticklabels(), rotation=0)
    axmatrix.set_xticklabels(axmatrix.get_xticklabels(), rotation=90)
    # axmatrix.yaxis.set_label_position("right")
    # axmatrix.yaxis.tick_right()
    axmatrix.spines['top'].set_visible(True)
    axmatrix.spines['right'].set_visible(True)
    axmatrix.spines['bottom'].set_visible(True)
    axmatrix.spines['left'].set_visible(True)

    # Add colorbar
    axcolor = fig.add_axes([0.2, 0.02, 0.65, 0.02])  # (left, bottom, width, height)
    cbar = plt.colorbar(
        mappable=axmatrix.collections[0], cax=axcolor, orientation="horizontal",
        ticks=[vmin, -4.5, -3] + [np.log10(pval_cut)] + [0] + [-np.log10(pval_cut)] + [5, 10] + [vmax])
    cbar.ax.set_xlabel('signed -log10(p-adj. values)', rotation=0, fontsize=12)

    # make a legend:
    for pw in pws:
        axmatrix.scatter([], [], s=pw, c="k", label=reverse_circle_sizes[pw])
    h, l = axmatrix.get_legend_handles_labels()
    fig.legend(h, l, labelspacing=1.2, title="Effect size", borderpad=1, loc='center left', bbox_to_anchor=(0.8, 0.5),
               frameon=True, framealpha=0.6, edgecolor="k", facecolor="w")

    plt.savefig(os.path.join(save_folder, "Scatterplot_padj_{}_ok.pdf".format(pval_cut)),
                bbox_inches='tight')
    plt.close()


if __name__ == '__main__':
    general_dir = os.path.join(
        '/Volumes', 'CH__data', 'Projects', 'Eyerich_AG_projects', 'BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer')
    output_dir = os.path.join(general_dir, 'analysis', 'Molecular_subtypes', 'output',
                              'Figure_S3_Waterfallplot', str(date.today()))
    os.makedirs(output_dir, exist_ok=True)

    data_root = os.path.join(general_dir, 'analysis', 'Molecular_subtypes', 'input', 'h5_files', 'LESION')

    meta_dir = os.path.join(general_dir, 'raw_data', 'clinical_data')

    main(save_folder=output_dir, input_dir=data_root, input_meta_dir=meta_dir)
