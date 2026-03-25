import numpy as np
import scanpy as sc
from datetime import date
import pickle
import platform
import pandas as pd

import igraph
from sklearn.metrics.cluster import adjusted_mutual_info_score

import itertools

import os
from matplotlib import pyplot as plt
import seaborn as sns


plt.ion()

if platform.system() == 'Linux':
    data_root = '/home/christina/2022_BRAIN_Eyerich'
else:
    data_root = '/Users/christina.hillig/PycharmProjects/Eyerich_ERCGrant'

fontsize_labels = 16
fontsize_ticks = 14


def compare_clustering(label_1, label_2, score_dict):
    # The lower the score the more similar are the cluserting results, score of 0 means they are identical
    # Desired output would be a low mean and std of the score across the seeds for a given resolution
    vi = igraph.compare_communities(comm1=label_1, comm2=label_2)
    # Adjusted rand index
    # The adjusted Rand index is thus ensured to have a value close to 0.0 for random labeling independently of
    # the number of clusters and samples and exactly 1.0 when the clusterings are identical (up to a permutation).
    # The adjusted Rand index is bounded below by -0.5 for especially discordant clusterings.
    ari = igraph.compare_communities(comm1=label_1, comm2=label_2, method='adjusted_rand')
    # adjusted mutual information score
    # The AMI returns a value of 1 when the two partitions are identical (ie perfectly matched).
    # Random partitions (independent labellings) have an expected AMI around 0 on average hence can be negative.
    # The value is in adjusted nats (based on the natural logarithm).
    amis = adjusted_mutual_info_score(label_1, label_2)

    score_dict['variation of information'].append(vi)
    score_dict['ARI'].append(ari)
    score_dict['AMI'].append(amis)

    return score_dict


def plot_compare_cluster_seeds(
        mean_score, std_scores, xvalues, optpercentage, save_folder, metric, key, mark_xvalue=None):
    fig, ax = plt.subplots()
    ax.grid(linewidth=0.5)
    for res in xvalues:
        ax.errorbar(res, mean_score[res][metric], yerr=std_scores[res][metric], fmt='-o', alpha=1, color='black',
                    lw=0.8, capsize=4)
    # Mark highest overall acc for each method with vertical line
    ax.axhline(y=0, xmin=-0.02, xmax=1.02, color='red', ls='--', lw=1)
    ax.axvline(x=xvalues[mark_xvalue], ymin=-1, ymax=1, color='grey', ls='--', lw=1)
    ax.set_xlabel(key, fontsize=fontsize_labels)
    ax.set_ylabel(metric, fontsize=fontsize_labels)
    if key == 'n clusters':
        ax.set_xlim([np.min(xvalues) - 0.5, np.max(xvalues) + 0.5])
        ax.set_xticks(np.arange(np.min(xvalues),  np.max(xvalues), 2))
    else:
        ax.set_xlim([0, np.max(xvalues) + 0.1])
    # Increase ticks label size
    ax.tick_params(axis='both', which='major', labelsize=fontsize_ticks)
    sns.despine(fig=fig, ax=ax)
    plt.tight_layout()

    fig.savefig(os.path.join(
        save_folder, 'GridSearch_compare_seeds_optPercentage{}_metric{}_{}.png'.format(optpercentage, metric, key)))
    plt.close(fig=fig)


def plot_metric(mean_score, std_scores, xvalues, optpercentage, save_folder, metric, key, mark_xvalue=None):

    fig, ax = plt.subplots()
    ax.grid(linewidth=0.5)
    ax.errorbar(xvalues, mean_score, yerr=std_scores, fmt='-o', alpha=1, color='black', lw=0.8, capsize=4)
    # Mark highest overall acc for each method with vertical line
    if mark_xvalue is not None:
        ax.axvline(x=xvalues[mark_xvalue], ymin=0, ymax=5, color='darkblue', ls='--', lw=1)
    ax.set_xlabel(key, fontsize=fontsize_labels)
    ax.set_ylabel(metric, fontsize=fontsize_labels)
    if key == 'n clusters':
        ax.set_xlim([np.min(xvalues) - 0.5, np.max(xvalues) + 0.5])
        ax.set_xticks(np.arange(np.min(xvalues),  np.max(xvalues), 2))
    else:
        ax.set_xlim([0, np.max(xvalues) + 0.1])
    # Increase ticks label size
    ax.tick_params(axis='both', which='major', labelsize=fontsize_ticks)
    sns.despine(fig=fig, ax=ax)
    plt.tight_layout()

    fig.savefig(os.path.join(
        save_folder, 'GridSearch_optPercentage{}_metric{}_{}.png'.format(optpercentage, metric, key)))
    plt.close(fig=fig)


def plot_resolution_vs_nclusters(mean_nclusters, std_nclusters, xvalues, save_folder):
    fig, ax = plt.subplots()
    ax.grid(linewidth=0.5)
    ax.errorbar(xvalues, mean_nclusters, yerr=std_nclusters,
                fmt='-o', alpha=1, color='black', lw=0.8, capsize=4)
    ax.set_xlabel('resolution', fontsize=fontsize_labels)
    ax.set_ylabel('mean number of clusters', fontsize=fontsize_labels)
    ax.set_xlim([0, np.max(xvalues) + 0.1])
    # Increase ticks label size
    ax.tick_params(axis='both', which='major', labelsize=fontsize_ticks)
    sns.despine(fig=fig, ax=ax)
    plt.tight_layout()
    fig.savefig(os.path.join(
        save_folder, 'GridSearch_optPercentage{}_meannumberofclusters.png'.format(7)))
    plt.close(fig=fig)


def main():
    path_input = "/Users/christina.hillig/PycharmProjects/Eyerich_ERCGrant/Molecular_subtypes/output/Leiden"
    output_dir = os.path.join(path_input, 'compare_resolutions', str(date.today()))
    os.makedirs(output_dir, exist_ok=True)
    file_hig = os.path.join(
        path_input, '2023-03-30', 'resolution_[0.5 0.6 0.7 0.8 0.9]__featureselection_HIG',
        'Scores_GridSearch_9_minPercentage7_maxPercentage7_minResolution0.5_maxResolution0.9_.pkl')
    with open(file_hig, 'rb') as ff:
        gs_hig = pickle.load(ff)

    df_res = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in gs_hig.items()]))

    seeds = df_res['seed'].unique()

    # TODO run application with optimal values
    adata = sc.read(os.path.join(
        data_root, 'Molecular_subtypes', 'input',
        'data_freeze', 'Lesion_RNAseq_20210720_patient_meta_data_v04__CH_20220209_TherapyResponse_v04_v220715.h5'))

    # TODO check out resolution 0.8 is approx. 10 clusters
    # Select opt resolution and save to dataframe
    for opt_res in [0.5, 0.6, 0.7, 0.8, 0.9]:
        df_res_opt = df_res.loc[df_res['resolution'] == opt_res, :]
        df_res_opt['MUC IDs'] = None
        df_res_opt['Gene names'] = None
        for seed in seeds:
            ind_seed = df_res_opt[df_res_opt['seed'] == seed].index[0]
            df_res_opt.at[ind_seed, 'MUC IDs'] = list(adata.obs.index)
            df_res_opt.at[ind_seed, 'Gene names'] = list(adata.var_names[df_res_opt.at[ind_seed, 'index selected genes']])

        with open(os.path.join(output_dir, 'Optimalres_{}__OptimalGPTK_7.pkl'.format(opt_res)), 'wb') as ff:
            pickle.dump(df_res_opt, ff)

        df_res_opt.to_excel(os.path.join(output_dir, 'Optimalres_{}__OptimalGPTK_7.xlsx'.format(opt_res)))


if __name__ == '__main__':
    main()
