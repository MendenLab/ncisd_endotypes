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
    data_root = '/Volumes/CH__data/Projects/Eyerich_AG_projects/BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer/analysis'

fontsize_labels = 18
fontsize_ticks = 16
fileformat = ".pdf"


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
        mean_score, std_scores, xvalues, optpercentage, save_folder, metric, key, xticks, mark_xvalue=None):
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.grid(linewidth=0.5)
    for res in xvalues:
        ax.errorbar(res, mean_score[res][metric], yerr=std_scores[res][metric], fmt='-o', alpha=1, color='black',
                    lw=0.8, capsize=4)
    # ax.axhline(y=0, xmin=-0.02, xmax=1.02, color='red', ls='--', lw=1)
    # Mark highest overall acc for each method with vertical line
    # ax.axvline(x=xvalues[mark_xvalue], ymin=-1, ymax=1, color='grey', ls='--', lw=1)
    ax.set_xlabel(key, fontsize=fontsize_labels)
    ax.set_ylabel(metric, fontsize=fontsize_labels)
    if key == 'n clusters':
        ax.set_xlim([np.min(xvalues) - 0.5, np.max(xvalues) + 0.5])
        ax.set_xticks(np.arange(np.min(xvalues),  np.max(xvalues), 2))
    else:
        ax.set_xlim([0.55, np.max(xvalues) + 0.05])
        ax.set_xticks(xticks)
        # ax.set_xticks(xticks[::2])
    if metric == 'AMI':
        ax.set_ylim([0.55, 0.75])
    # Increase ticks label size
    ax.tick_params(axis='both', which='major', labelsize=fontsize_ticks)
    sns.despine(fig=fig, ax=ax)
    plt.tight_layout()

    fig.savefig(os.path.join(
        save_folder, 'GridSearch_compare_seeds_optPercentage{}_metric{}_{}{}'.format(
            optpercentage, metric, key, fileformat)))
    plt.close(fig=fig)


def plot_metric(mean_score, std_scores, xvalues, optpercentage, save_folder, metric, key, xticks, mark_xvalue=None):

    fig, ax = plt.subplots(figsize=(6, 6))
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
        # ax.set_xlim([0, np.max(xvalues) + 0.1])
        ax.set_xlim([0.55, np.max(xvalues) + 0.1])
        ax.set_xticks(xticks)
    # Increase ticks label size
    ax.tick_params(axis='both', which='major', labelsize=fontsize_ticks)
    sns.despine(fig=fig, ax=ax)
    plt.tight_layout()

    fig.savefig(os.path.join(
        save_folder, 'GridSearch_optPercentage{}_metric{}_{}{}'.format(optpercentage, metric, key, fileformat)))
    plt.close(fig=fig)


def plot_resolution_vs_nclusters(mean_nclusters, std_nclusters, xvalues, save_folder, opt_percentage=7):
    fig, ax = plt.subplots(figsize=(6, 6))
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
        save_folder, 'GridSearch_optPercentage{}_meannumberofclusters{}'.format(opt_percentage, fileformat)))
    plt.close(fig=fig)


def main():
    path_input = "/Volumes/CH__data/Projects/Eyerich_AG_projects/BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer/analysis/Molecular_subtypes/output/Leiden"
    output_dir = os.path.join(path_input, 'compare_resolutions', str(date.today()))
    os.makedirs(output_dir, exist_ok=True)
    file_hig = os.path.join(
        path_input, '2023-02-20', 'resolution_[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.  1.1 1.2 1.3 1.4 1.5]__featureselection_HIG',
        'Scores_GridSearch_9_minPercentage7_maxPercentage7_minResolution0.1_maxResolution1.5_.pkl')
    with open(file_hig, 'rb') as ff:
        gs_hig = pickle.load(ff)

    opt_percentage = 7
    resolutions = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.5]

    df_res = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in gs_hig.items()]))

    resolution = df_res['resolution'].unique()
    seeds = df_res['seed'].unique()
    seed_pair_comb = list(itertools.combinations(seeds, 2))
    scores = {k: [] for k in resolution}
    for res in resolution:
        # ARI: Adjusted rand index; AMI: adjusted mutual information
        scores[res] = {'ARI': [], 'AMI': [], 'variation of information': []}
        for seed_1, seed_2 in seed_pair_comb:
            df_res_01 = df_res[(df_res['resolution'] == res) & (df_res['seed'] == seed_1)]
            df_res_02 = df_res[(df_res['resolution'] == res) & (df_res['seed'] == seed_2)]

            # Calculate scores comparig the different clustering results for varying seeds
            scores[res] = compare_clustering(
                label_1=df_res_01['Molecular subtypes'].values[0], label_2=df_res_02['Molecular subtypes'].values[0],
                score_dict=scores[res])

    # Find most stable (scores across seeds) resolution and well separated clusters (silhouette score)
    # 1. Find most stable resolution
    mean = {k: [] for k in resolution}
    std = {k: [] for k in resolution}
    for res in resolution:
        mean[res] = {'ARI': None, 'AMI': None, 'variation of information': None}
        std[res] = {'ARI': None, 'AMI': None, 'variation of information': None}
        for score_metrics in scores[res].keys():
            mean[res][score_metrics] = np.nanmean(scores[res][score_metrics])
            std[res][score_metrics] = np.nanstd(scores[res][score_metrics])

    tmp_scores = np.asarray([(v['variation of information'], v['ARI'], v['AMI']) for k, v in mean.items()])
    vi_score, ari_score, amis_score = tmp_scores.T[0], tmp_scores.T[1], tmp_scores.T[2]

    # ARI
    plot_compare_cluster_seeds(
        mean_score=mean, std_scores=std, xvalues=resolution, optpercentage=opt_percentage, save_folder=output_dir,
        metric='ARI', mark_xvalue=np.nanargmax(ari_score), key='resolution', xticks=resolutions)
    # AMI score
    plot_compare_cluster_seeds(
        mean_score=mean, std_scores=std, xvalues=resolution, optpercentage=opt_percentage, save_folder=output_dir,
        metric='AMI', mark_xvalue=np.nanargmax(amis_score), key='resolution', xticks=resolutions[5:])
    # VI score
    plot_compare_cluster_seeds(
        mean_score=mean, std_scores=std, xvalues=resolution, optpercentage=opt_percentage, save_folder=output_dir,
        metric='variation of information', mark_xvalue=np.nanargmin(vi_score), key='resolution', xticks=resolutions)

    # 2. Find best resolution
    metrics_mean = {k: [] for k in resolution}
    metrics_std = {k: [] for k in resolution}
    for res in resolution:
        metrics_mean[res] = {'Davies Bouldin Score': None, 'Calinski-Harabasz Index': None, 'Silhouette Score': None,
                             'Density Based clustering validation': None}
        metrics_std[res] = {'Davies Bouldin Score': None, 'Calinski-Harabasz Index': None, 'Silhouette Score': None,
                            'Density Based clustering validation': None}
        df_res_01 = df_res[(df_res['resolution'] == res)]
        for metric_name in metrics_mean[res].keys():
            metrics_mean[res][metric_name] = np.nanmean(df_res_01[metric_name].values)
            metrics_std[res][metric_name] = np.nanstd(df_res_01[metric_name].values)

    tmp_metrics = np.asarray(
        [(v['Davies Bouldin Score'], v['Calinski-Harabasz Index'],
          v['Silhouette Score'], v['Density Based clustering validation']) for k, v in metrics_mean.items()])
    mean_db_score, mean_ch_score, mean_s_score, mean_dbcv_score = tmp_metrics.T[0], tmp_metrics.T[1], \
                                                                  tmp_metrics.T[2], tmp_metrics.T[3]
    tmp_metrics = np.asarray(
        [(v['Davies Bouldin Score'], v['Calinski-Harabasz Index'],
          v['Silhouette Score'], v['Density Based clustering validation']) for k, v in metrics_std.items()])
    std_db_score, std_ch_score, std_s_score, std_dbcv_score = tmp_metrics.T[0], tmp_metrics.T[1], \
                                                              tmp_metrics.T[2], tmp_metrics.T[3]

    # Average change of DB for gamma > 0.5
    steps = []
    for ind, val in enumerate(mean_db_score[5:-1]):
        steps.append(mean_db_score[5:][ind+1] - mean_db_score[5:][ind])
    np.mean(steps)
    np.std(np.asarray(steps)[np.asarray(steps) > 0])
    np.std(np.asarray(steps)[np.asarray(steps) < 0])

    plot_metric(mean_score=mean_db_score[5:], std_scores=std_db_score[5:], xvalues=resolution[5:],
                optpercentage=opt_percentage, save_folder=output_dir, metric='Davies Bouldin index',
                mark_xvalue=np.nanargmin(mean_db_score), xticks=resolution[5:], key='resolution')
    plot_metric(mean_score=mean_ch_score, std_scores=std_ch_score, xvalues=resolution,
                optpercentage=opt_percentage, save_folder=output_dir, metric='Calinski-Harabasz index',
                mark_xvalue=np.nanargmax(mean_ch_score), xticks=resolution, key='resolution')
    plot_metric(mean_score=mean_s_score, std_scores=std_s_score, xvalues=resolution,
                optpercentage=opt_percentage, save_folder=output_dir, metric='Silhouette coefficient',
                mark_xvalue=None, xticks=resolution, key='resolution')
    # index produces values between −1 and +1, with greater values of the measure indicating
    # better density-based clustering solutions.
    plot_metric(mean_score=mean_dbcv_score, std_scores=std_dbcv_score, xvalues=resolution,
                optpercentage=opt_percentage, save_folder=output_dir, metric='DBCV',
                mark_xvalue=np.nanargmax(mean_dbcv_score), xticks=resolution, key='resolution')

    # Find resolutions with same number of clusters adn compare these

    res_mean_number_clusters = []
    res_std_number_clusters = []
    num_clusters = {k: [] for k in resolution}
    df_res['n clusters'] = 0
    for res in resolution:
        df_res_01 = df_res[df_res['resolution'] == res]
        lables = df_res_01['Molecular subtypes'].values
        num_clusters[res] = []
        for label in lables:
            num_clusters[res].append(len(np.unique(label)))
            # print("{}".format(res), len(np.unique(label)))
        res_mean_number_clusters.append(np.nanmean(num_clusters[res]))
        res_std_number_clusters.append(np.nanstd(num_clusters[res]))
        df_res.loc[df_res['resolution'] == res, 'n clusters'] = num_clusters[res]

    # Plot resolution vs number of clusters
    plot_resolution_vs_nclusters(mean_nclusters=res_mean_number_clusters, std_nclusters=res_std_number_clusters,
                                 xvalues=resolution, save_folder=output_dir, opt_percentage=opt_percentage)

    # Apply same metrics again but using number of clusters instead of resolution
    nclusters = df_res['n clusters'].unique()
    nclusters.sort()
    nclusters = list(nclusters)
    # remove 1 from list
    nclusters.remove(1)
    scores_ncluster = {k: [] for k in nclusters}
    for ncluster in nclusters:
        scores_ncluster[ncluster] = {'ARI': [], 'AMI': [],
                                     'variation of information': []}
        for seed_1, seed_2 in seed_pair_comb:
            df_res_01 = df_res[(df_res['n clusters'] == ncluster) & (df_res['seed'] == seed_1)]
            df_res_02 = df_res[(df_res['n clusters'] == ncluster) & (df_res['seed'] == seed_2)]

            # Calculate scores comparig the different clustering results for varying seeds
            if not df_res_01.empty and not df_res_02.empty:
                scores_ncluster[ncluster] = compare_clustering(
                    label_1=df_res_01['Molecular subtypes'].values[0],
                    label_2=df_res_02['Molecular subtypes'].values[0],
                    score_dict=scores_ncluster[ncluster])
            else:
                scores_ncluster[ncluster]['ARI'].append(np.nan)
                scores_ncluster[ncluster]['AMI'].append(np.nan)
                scores_ncluster[ncluster]['variation of information'].append(np.nan)

    mean_ncluster_scores = {k: [] for k in nclusters}
    std_ncluster_scores = {k: [] for k in nclusters}
    for ncluster in nclusters:
        mean_ncluster_scores[ncluster] = {'ARI': None, 'AMI': None, 'variation of information': None}
        std_ncluster_scores[ncluster] = {'ARI': None, 'AMI': None, 'variation of information': None}
        for score_metrics in scores_ncluster[ncluster].keys():
            mean_ncluster_scores[ncluster][score_metrics] = np.nanmean(scores_ncluster[ncluster][score_metrics])
            std_ncluster_scores[ncluster][score_metrics] = np.nanstd(scores_ncluster[ncluster][score_metrics])

    tmp_scores = np.asarray([(v['variation of information'], v['ARI'],
                              v['AMI']) for k, v in mean_ncluster_scores.items()])
    vi_ncluster_score, ari_ncluster_score, amis_ncluster_score = tmp_scores.T[0], tmp_scores.T[1], tmp_scores.T[2]

    # ARI
    plot_compare_cluster_seeds(
        mean_score=mean_ncluster_scores, std_scores=std_ncluster_scores, xvalues=nclusters, optpercentage=7,
        save_folder=output_dir,
        metric='ARI', mark_xvalue=np.nanargmax(ari_ncluster_score), key='n clusters', xticks=resolutions)
    # AMI score
    plot_compare_cluster_seeds(
        mean_score=mean_ncluster_scores, std_scores=std_ncluster_scores, xvalues=nclusters, optpercentage=7,
        save_folder=output_dir,
        metric='AMI', mark_xvalue=np.nanargmax(amis_ncluster_score), key='n clusters', xticks=resolutions)
    # VI score
    plot_compare_cluster_seeds(
        mean_score=mean_ncluster_scores, std_scores=std_ncluster_scores, xvalues=nclusters, optpercentage=7,
        save_folder=output_dir,
        metric='variation of information', mark_xvalue=np.nanargmin(vi_ncluster_score),
        key='n clusters', xticks=resolutions)

    # Find optimal number of clusters
    metrics_mean_nclusters = {k: [] for k in nclusters}
    metrics_std_nclusters = {k: [] for k in nclusters}
    for ncluster in nclusters:
        metrics_mean_nclusters[ncluster] = {'Davies Bouldin Score': None, 'Calinski-Harabasz Index': None,
                                            'Silhouette Score': None, 'Density Based clustering validation': None}
        metrics_std_nclusters[ncluster] = {'Davies Bouldin Score': None, 'Calinski-Harabasz Index': None,
                                           'Silhouette Score': None, 'Density Based clustering validation': None}
        df_res_01 = df_res[(df_res['n clusters'] == ncluster)]
        for metric_name in metrics_mean_nclusters[ncluster].keys():
            metrics_mean_nclusters[ncluster][metric_name] = np.nanmean(df_res_01[metric_name].values)
            metrics_std_nclusters[ncluster][metric_name] = np.nanstd(df_res_01[metric_name].values)

    tmp_metrics = np.asarray(
        [(v['Davies Bouldin Score'], v['Calinski-Harabasz Index'],
          v['Silhouette Score'], v['Density Based clustering validation']) for k, v in metrics_mean_nclusters.items()])
    mean_db_score_nclusters, mean_ch_score_nclusters, \
    mean_s_score_nclusters, mean_dbcv_score_nclusters = tmp_metrics.T[0], tmp_metrics.T[1], \
                                                        tmp_metrics.T[2], tmp_metrics.T[3]
    tmp_metrics = np.asarray(
        [(v['Davies Bouldin Score'], v['Calinski-Harabasz Index'],
          v['Silhouette Score'], v['Density Based clustering validation']) for k, v in metrics_std_nclusters.items()])
    std_db_score_nclusters, std_ch_score_nclusters, \
    std_s_score_nclusters, std_dbcv_score_nclusters = tmp_metrics.T[0], tmp_metrics.T[1], \
                                                      tmp_metrics.T[2], tmp_metrics.T[3]

    plot_metric(mean_score=mean_db_score_nclusters, std_scores=std_db_score_nclusters, xvalues=nclusters,
                optpercentage=opt_percentage, save_folder=output_dir, metric='Davies Bouldin index',
                mark_xvalue=np.nanargmin(mean_db_score_nclusters), xticks=nclusters, key='n clusters')
    plot_metric(mean_score=mean_ch_score_nclusters, std_scores=std_ch_score_nclusters, xvalues=nclusters,
                optpercentage=opt_percentage, save_folder=output_dir, metric='Calinski-Harabasz index',
                mark_xvalue=np.nanargmax(mean_ch_score_nclusters), xticks=nclusters, key='n clusters')
    plot_metric(mean_score=mean_s_score_nclusters, std_scores=std_s_score_nclusters, xvalues=nclusters,
                optpercentage=opt_percentage, save_folder=output_dir, metric='Silhouette coefficient',
                mark_xvalue=None, xticks=nclusters, key='n clusters')
    plot_metric(mean_score=mean_dbcv_score_nclusters, std_scores=std_dbcv_score_nclusters, xvalues=nclusters,
                optpercentage=opt_percentage, save_folder=output_dir, metric='DBCV',
                mark_xvalue=np.nanargmax(mean_dbcv_score_nclusters), xticks=nclusters, key='n clusters')

    # TODO run application with optimal values
    adata = sc.read(os.path.join(
        data_root, 'Molecular_subtypes', 'input', 'h5_files', 'LESION',
        'Lesion_RNAseq_20210720_patient_meta_data_v04__CH_20220209_TherapyResponse_v04_v220715.h5'))

    # TODO check out resolution 0.8 is approx. 10 clusters
    # Select opt resolution and save to dataframe
    for opt_res in [0.5, 0.6, 0.7, 0.8, 0.9]:
        df_res_opt = df_res.loc[df_res['resolution'] == opt_res, :].copy()
        df_res_opt['MUC IDs'] = None
        df_res_opt['Gene names'] = None
        for seed in seeds:
            ind_seed = df_res_opt[df_res_opt['seed'] == seed].index[0]
            df_res_opt.at[ind_seed, 'MUC IDs'] = list(adata.obs.index)
            df_res_opt.at[ind_seed, 'Gene names'] = list(adata.var_names[df_res_opt.at[ind_seed, 'index genes']])

        with open(os.path.join(output_dir, 'Optimalres_{}__OptimalGPTK_{}.pkl'.format(
                opt_res, opt_percentage)), 'wb') as ff:
            pickle.dump(df_res_opt, ff)

        df_res_opt.to_excel(os.path.join(output_dir, 'Optimalres_{}__OptimalGPTK_{}.xlsx'.format(
            opt_res, opt_percentage)))


if __name__ == '__main__':
    main()
