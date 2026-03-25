# Create spiderplot of responder and non responder using the identified 15 genes
# A spiderplot shows the mean of a genes expression count, scaled between 0 and 1
# add on top the Michigan profile cohort
from scripts.GEx_classifier.config import config_tnf

from sklearn.metrics.pairwise import cosine_similarity
from scipy.stats import pearsonr
from datetime import date
import copy
import scanpy as sc
import pandas as pd
import numpy as np
from math import pi
import pickle
import os

import matplotlib.pyplot as plt
import seaborn as sns

# read in DEGs of L vs NL Michigan and Eyerich cohort
# for Eyerich and Michigan cohort:
# subset to classifier genes
# plot spiderplot
# optional: scale mean counts between 0 and 1


def NormalizeData(data):
    return (data - np.min(data)) / (np.max(data) - np.min(data))


def plot_radar_plot(values, cols, ax, color, alpha=1.):
    N = len(cols)

    # What will be the angle of each axis in the plot? (we divide the plot / number of variable)
    angles = [n / float(N) * 2 * pi for n in range(N)]
    angles += angles[:1]

    # Draw one axe per variable + add labels
    plt.xticks(angles[:-1], cols, color='k', size=16)
    # Draw ylabels
    ax.set_rlabel_position(0)
    plt.yticks([0.2, 0.4, 0.6, 0.8], ["0.2", "0.4", "0.6", "0.8"], color="black", size=12)
    ax.set_ylim(0, 1)
    # Plot data
    ax.plot(angles, values, linewidth=1, linestyle='solid', color=color)
    # Fill area
    ax.fill(angles, values, color, alpha=alpha)

    return ax


def plot_profile_plots(df_group_1, df_group_2, cols, color, save_folder, key):
    group_1_values = df_group_1['log2FoldChange'].values.flatten().tolist()
    # scale between 0 and 1
    group_1_values = list(NormalizeData(data=group_1_values))
    normed_group_1_values = copy.copy(group_1_values)
    group_1_values += group_1_values[:1]

    group_2_values = df_group_2['log2FoldChange'].values.flatten().tolist()
    # scale between 0 and 1
    group_2_values = list(NormalizeData(data=group_2_values))
    normed_group_2_values = copy.copy(group_2_values)
    group_2_values += group_2_values[:1]

    # # Calculate similarity
    # intersection = len(set(normed_group_1_values) & set(normed_group_2_values))
    # union = len(set(normed_group_1_values) | set(normed_group_2_values))
    # jaccard_similarity = intersection / union
    # print(f"Jaccard Similarity: {jaccard_similarity}")
    #
    # # Convert lists to vectors
    # normed_group_1_vector = np.array(normed_group_1_values).reshape(1, -1)
    # normed_group_2_vector = np.array(normed_group_2_values).reshape(1, -1)
    # print('Cosine Similarity', cosine_similarity(normed_group_1_vector, normed_group_2_vector)[0, 0])
    #
    # hamming_distance = sum(el1 != el2 for el1, el2 in zip(normed_group_1_values, normed_group_2_values))
    # hamming_similarity = 1 - (hamming_distance / len(normed_group_2_values))
    # print(f"Hamming Similarity: {hamming_similarity}")

    correlation_coefficient, _ = pearsonr(normed_group_1_values, normed_group_2_values)
    print(f"Correlation Coefficient: {correlation_coefficient}")

    # df_plot = pd.DataFrame(data=np.asarray([normed_group_1_values, normed_group_2_values]).T,
    #                        columns=['Michigan', 'Eyerich'])
    # sns.scatterplot(data=df_plot, x='Michigan', y='Eyerich')

    # Initialise the spider plot
    ax = plt.subplot(111, polar=True)
    ax = plot_radar_plot(values=group_1_values, cols=cols, ax=ax, color=color[0])
    ax = plot_radar_plot(values=group_2_values, cols=cols, ax=ax, color=color[1], alpha=0.4)
    plt.text(s="{}={:.2f}".format(r'$r_{pearson}$', correlation_coefficient), y=1, x=1, fontsize=16)
    plt.savefig(os.path.join(save_folder, '{}.pdf'.format(key)))
    plt.close()


def plot_boxplot(df, obs, palette, save_folder, key):
    df__normed_boxplot = df.copy()
    df__normed_boxplot['Endotypes'] = obs
    df_boxplot = pd.melt(df__normed_boxplot, id_vars="Endotypes")
    df_boxplot["variable"] = df_boxplot["variable"].astype("category")

    df_boxplot["Endotypes"] = df_boxplot["Endotypes"].astype("category")
    df_boxplot['Endotypes'] = df_boxplot['Endotypes'].cat.reorder_categories(['responder', 'non-responder'])

    fig, ax = plt.subplots(figsize=(8, 6))
    sns.boxplot(data=df_boxplot, x="variable", y="value", hue='Endotypes',
                ax=ax, color=".8", palette=palette)
    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False)
    sns.despine(fig=fig, ax=ax)
    ax.set_ylabel("scaled counts", fontsize=12)
    ax.set_xlabel("")
    plt.setp(ax.get_xticklabels(), ha="center", rotation=90, fontsize=12)
    plt.savefig(
        os.path.join(save_folder, "Boxplot_{}_cohort.pdf".format(key)),
        bbox_inches="tight",
    )
    plt.close(fig=fig)


def load_model(data_root):
    # Load configurations -> see config_TNF.py file
    cfg = config_tnf.Config()

    # Specify date - will load .pkl from folder e.g.
    # ../Molecular_subtypes/output/Classifier/TNF_IL23_HC_classifier_E13_vs_E12/2023-09-20
    dates = {cfg.dataset.name: '2024-01-08'}  #'2023-11-07'}

    # Load best model
    model = BestClassifier(model_name=cfg.classifier.name)
    model.load_model(date_func=dates[cfg.dataset.name], dataset_name=cfg.dataset.name,
                     fs_method=cfg.feature_selection.method, vif=cfg.feature_selection.vif.apply, data_root=data_root)

    return model



class BestClassifier:
    def __init__(self, model_name):
        self.cfg = None
        self.palette = sc.pl.palettes.default_20

        self.dict_classifier_res = None
        self.estimator_information = None
        self.best_model = None
        self.best_model_parameters = None
        self.fs = None
        self.parts = None

        self.label_encoder = None
        self.norm = None

        self.model_name = model_name

    def load_model(self: str, date_func: str, fs_method: str, dataset_name: str, vif: str, data_root: str):
        # Load the best classifier object
        filename = os.path.join(data_root, 'Molecular_subtypes', 'output', 'NestedCrossValidation_fs_classifier',
                                'Classifier_N1000', dataset_name, self.model_name,
                                'VIF_{}_FSmethod_{}'.format(vif, fs_method), date_func, 'Final_classifier.pkl')
        with open(filename, 'rb') as ff:
            self.estimator_information = pickle.load(ff)

        self.best_model = self.estimator_information['model_object'].model
        self.dict_classifier_res = self.estimator_information['gridsearchCV']

        self.fs = self.estimator_information['model_object'].fs
        self.cfg = self.estimator_information['model_object'].cfg
        self.parts = self.estimator_information['model_object'].parts

        self.label_encoder = self.estimator_information['model_object'].label_encoder
        self.norm = self.estimator_information['model_object'].norm


def main(save_folder):
    data_root = '/Volumes/CH__data/Projects/Eyerich_AG_projects/BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer/analysis'
    # Load model and read out features
    model = load_model(data_root=data_root)

    for responder_type in ['Responder', 'NONResponder']:
        if responder_type == 'Responder':
            endotype = 'E13'
        else:
            endotype = 'E12'
        df_degs_eyerich = pd.read_excel(os.path.join(
            data_root, 'Molecular_subtypes', 'output', 'DGE_analysis', 'TNF_responder_E13_TNF_nonresponder_E12_LNL_',
            '2024-01-24', 'DEGs_{}_LvsNL___Age_Sex_sampleType.xlsx'.format(endotype)))

        # Load DEGs Michigan
        df_degs_michigan = pd.read_excel(os.path.join(
            data_root, 'Molecular_subtypes', 'output', 'DGE_analysis',
            'Michigan_Responder_NONResponder_L_vs_NL', '2024-01-24',
            'DEGs_{}_LvsNL___Age_Sex_SampleType.xlsx'.format(responder_type)))

        # subset to classifier genes
        features = model.fs.get_feature_names_out()
        df_degs_michigan = df_degs_michigan.loc[df_degs_michigan['hgnc_symbol'].isin(features), :]
        df_degs_eyerich = df_degs_eyerich.loc[df_degs_eyerich['hgnc_symbol'].isin(features), :]

        # get intersection of both datasets
        shared_genes = np.intersect1d(df_degs_eyerich['hgnc_symbol'], df_degs_michigan['hgnc_symbol'])
        df_degs_michigan = df_degs_michigan.loc[df_degs_michigan['hgnc_symbol'].isin(shared_genes), :]
        df_degs_eyerich = df_degs_eyerich.loc[df_degs_eyerich['hgnc_symbol'].isin(shared_genes), :]

        # resort dataframes
        df_degs_michigan.index = df_degs_michigan['hgnc_symbol']
        df_degs_eyerich.index = df_degs_eyerich['hgnc_symbol']
        df_degs_michigan = df_degs_michigan.loc[df_degs_eyerich.index, :]
        assert np.all(df_degs_michigan.index == df_degs_eyerich.index)

        # Spiderplots
        plot_profile_plots(df_group_1=df_degs_michigan, df_group_2=df_degs_eyerich,
                           cols=list(df_degs_eyerich.index), color=['darkred', 'darkorange'], save_folder=save_folder,
                           key='{}_DEGs_Eyerich_vs_Michigan'.format(responder_type))

    # Load DEGs Eyerich
    df_degs_eyerich_responder = pd.read_excel(os.path.join(
        data_root, 'Molecular_subtypes', 'output', 'DGE_analysis', 'TNF_responder_E13_TNF_nonresponder_E12_LNL_',
        '2024-01-24', 'DEGs_{}_LvsNL___Age_Sex_sampleType.xlsx'.format('E13')))
    df_degs_eyerich_nonresponder = pd.read_excel(os.path.join(
        data_root, 'Molecular_subtypes', 'output', 'DGE_analysis', 'TNF_responder_E13_TNF_nonresponder_E12_LNL_',
        '2024-01-24', 'DEGs_{}_LvsNL___Age_Sex_sampleType.xlsx'.format('E12')))

    # Load DEGs Michigan
    df_degs_michigan_responder = pd.read_excel(os.path.join(
        data_root, 'Molecular_subtypes', 'output', 'DGE_analysis',
        'Michigan_Responder_NONResponder_L_vs_NL', '2024-01-24',
        'DEGs_{}_LvsNL___Age_Sex_SampleType.xlsx'.format('Responder')))
    df_degs_michigan_nonresponder = pd.read_excel(os.path.join(
        data_root, 'Molecular_subtypes', 'output', 'DGE_analysis',
        'Michigan_Responder_NONResponder_L_vs_NL', '2024-01-24',
        'DEGs_{}_LvsNL___Age_Sex_SampleType.xlsx'.format('NONResponder')))

    # subset to classifier genes
    features = model.fs.get_feature_names_out()
    df_degs_eyerich_responder = df_degs_eyerich_responder.loc[df_degs_eyerich_responder['hgnc_symbol'].isin(features), :]
    df_degs_eyerich_nonresponder = df_degs_eyerich_nonresponder.loc[df_degs_eyerich_nonresponder['hgnc_symbol'].isin(features), :]
    df_degs_michigan_responder = df_degs_michigan_responder.loc[df_degs_michigan_responder['hgnc_symbol'].isin(features), :]
    df_degs_michigan_nonresponder = df_degs_michigan_nonresponder.loc[df_degs_michigan_nonresponder['hgnc_symbol'].isin(features), :]

    # get intersection of both datasets
    shared_genes_eyerich = np.intersect1d(df_degs_eyerich_responder['hgnc_symbol'], df_degs_eyerich_nonresponder['hgnc_symbol'])
    shared_genes_michigan = np.intersect1d(df_degs_michigan_responder['hgnc_symbol'], df_degs_michigan_nonresponder['hgnc_symbol'])
    shared_genes = np.intersect1d(shared_genes_eyerich, shared_genes_michigan)
    df_degs_eyerich_responder = df_degs_eyerich_responder.loc[df_degs_eyerich_responder['hgnc_symbol'].isin(shared_genes), :]
    df_degs_eyerich_nonresponder = df_degs_eyerich_nonresponder.loc[df_degs_eyerich_nonresponder['hgnc_symbol'].isin(shared_genes), :]
    df_degs_michigan_responder = df_degs_michigan_responder.loc[df_degs_michigan_responder['hgnc_symbol'].isin(shared_genes), :]
    df_degs_michigan_nonresponder = df_degs_michigan_nonresponder.loc[df_degs_michigan_nonresponder['hgnc_symbol'].isin(shared_genes), :]

    # resort dataframes
    df_degs_eyerich_responder.index = df_degs_eyerich_responder['hgnc_symbol']
    df_degs_eyerich_nonresponder.index = df_degs_eyerich_nonresponder['hgnc_symbol']
    df_degs_michigan_responder.index = df_degs_michigan_responder['hgnc_symbol']
    df_degs_michigan_nonresponder.index = df_degs_michigan_nonresponder['hgnc_symbol']
    df_degs_eyerich_nonresponder = df_degs_eyerich_nonresponder.loc[df_degs_eyerich_responder.index, :]
    df_degs_michigan_responder = df_degs_michigan_responder.loc[df_degs_eyerich_responder.index, :]
    df_degs_michigan_nonresponder = df_degs_michigan_nonresponder.loc[df_degs_eyerich_responder.index, :]
    assert np.all(df_degs_eyerich_responder.index == df_degs_eyerich_nonresponder.index)
    assert np.all(df_degs_michigan_responder.index == df_degs_michigan_nonresponder.index)

    # Spiderplots
    plot_profile_plots(df_group_1=df_degs_eyerich_responder, df_group_2=df_degs_eyerich_nonresponder,
                       cols=list(df_degs_eyerich_responder.index), color=['darkgreen', 'darkblue'], save_folder=save_folder,
                       key='Responder_vs_NONResponder_DEGs_Eyerich')
    plot_profile_plots(df_group_1=df_degs_michigan_responder, df_group_2=df_degs_michigan_nonresponder,
                       cols=list(df_degs_eyerich_responder.index), color=['darkgreen', 'darkblue'], save_folder=save_folder,
                       key='Responder_vs_NONResponder_DEGs_Michigan')


if __name__ == '__main__':
    bc = True
    imputing = False
    output_dir = os.path.join(
        '/Volumes', 'CH__data', 'Projects', 'Eyerich_AG_projects',
        'BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer', 'analysis', 'Molecular_subtypes',
        'output', 'Figure_5_Michigan', "DEGs_spiderplots".format(
            str(date.today()), bc, imputing))
    os.makedirs(output_dir, exist_ok=True)

    main(save_folder=output_dir)



