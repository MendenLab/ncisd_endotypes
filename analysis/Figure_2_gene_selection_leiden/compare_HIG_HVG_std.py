import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import os
from datetime import date
import pickle
import seaborn as sns
import collections
import pandas as pd
from collections import Counter
from scipy.stats import mannwhitneyu
import platform

from sklearn.decomposition import PCA
# PCA.explained_variance_ratio_()

# import the necessary functions from rpy2
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
pandas2ri.activate()

# import the BayesFactor package
BayesFactor = importr('BayesFactor')
base = importr("base")

fontsize_title = 18
fontsize_labels = 18
fontsize_ticks = 15
fileformat = 'pdf'


if platform.system() == 'Linux':
    data_root = '/home/christina/2022_BRAIN_Eyerich'
else:
    data_root = '/Volumes/CH__data/Projects/Eyerich_AG_projects/BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer/analysis'


def pandasToR(df: pd.DataFrame):
    """
    Transforms a pandas counts to R
    Parameters
    ----------
    df: Pandas counts

    Returns
    -------
    R counts corresponding to given pandas Dataframe

    """
    with localconverter(ro.default_converter + pandas2ri.converter):
        return ro.conversion.py2rpy(df)


def plot_errorbars(dict_res, save_folder, default_title):
    unique_percentage = np.unique(dict_res['gene percentage to keep'])
    num_unique_percentage = len(unique_percentage)
    num_unique_resolutions = np.unique(dict_res['resolution'])
    # per cluster
    max_mean_acc_cluster = dict.fromkeys(list(num_unique_resolutions))
    counter = 0
    for ind in range(0, len(num_unique_resolutions)):
        fig, ax = plt.subplots()
        ax.grid(linewidth=0.5)
        ax.errorbar(np.unique(dict_res['gene percentage to keep']),
                    dict_res['mean_overall_acc'][counter:counter + num_unique_percentage],
                    yerr=dict_res['std_overall_acc'][counter:counter + num_unique_percentage], fmt='-o')
        ax.set_ylabel('mean accuracy', fontsize=fontsize_labels)
        ax.set_xlabel('gene percentage to keep', fontsize=fontsize_labels)
        ax.set_title('Leiden resolution {}'.format(num_unique_resolutions[ind]), fontsize=fontsize_title)
        # Increase ticks label size
        ax.tick_params(axis='both', which='major', labelsize=fontsize_ticks)
        sns.despine()
        fig.savefig(os.path.join(
            save_folder,
            'MeanAccOverallSeeds_{}_Resolution_{}.{}'.format(default_title, num_unique_resolutions[ind], fileformat)))
        plt.close(fig=fig)

        ind_max_mean_acc = np.nanargmax(dict_res['mean_overall_acc'][counter:counter + num_unique_percentage])
        max_mean_acc_cluster[num_unique_resolutions[ind]] = (
            dict_res['mean_overall_acc'][counter:counter + num_unique_percentage][ind_max_mean_acc],
            unique_percentage[ind_max_mean_acc])

        counter += num_unique_percentage

    return max_mean_acc_cluster


def plot_errorbars_HIG_vs_HVG(dict_res_hvg, dict_res_hig, save_folder, default_title):
    """
    Mean acc = acc overall seeds
    Overall Mean acc = Mean acc over all resolutions

    Parameters
    ----------
    dict_res_hvg
    dict_res_hig
    save_folder
    default_title

    Returns
    -------

    """
    unique_percentage = np.unique(dict_res_hvg['gene percentage to keep'])
    num_unique_percentage = len(unique_percentage)
    num_unique_resolutions = np.unique(dict_res_hvg['resolution'])

    # per cluster
    counter = 0
    for ind in range(0, len(num_unique_resolutions)):
        ind_max_mean_acc_hvg = np.nanargmax(dict_res_hvg['mean_overall_acc'][counter:counter + num_unique_percentage])
        ind_max_mean_acc_hig = np.nanargmax(dict_res_hig['mean_overall_acc'][counter:counter + num_unique_percentage])

        stat, p_value = mannwhitneyu(dict_res_hig['mean_overall_acc'][counter:counter + num_unique_percentage],
                                     dict_res_hvg['mean_overall_acc'][counter:counter + num_unique_percentage],
                                     alternative='greater')
        fig, ax = plt.subplots()
        ax.grid(linewidth=0.5)
        ax.errorbar(np.unique(dict_res_hvg['gene percentage to keep']),
                    dict_res_hvg['mean_overall_acc'][counter:counter + num_unique_percentage],
                    yerr=dict_res_hvg['std_overall_acc'][counter:counter + num_unique_percentage], fmt='-o',
                    color='darkblue', alpha=0.7, label='HVG')
        ax.errorbar(np.unique(dict_res_hig['gene percentage to keep']),
                    dict_res_hig['mean_overall_acc'][counter:counter + num_unique_percentage],
                    yerr=dict_res_hig['std_overall_acc'][counter:counter + num_unique_percentage], fmt='-o',
                    color='darkorange', alpha=0.7, label='HIG')
        # Mark highest overall acc for each method with vertical line
        ax.axvline(x=unique_percentage[ind_max_mean_acc_hvg], ymin=0, ymax=1, color='darkblue', ls='--', lw=0.5)
        ax.text(x=unique_percentage[ind_max_mean_acc_hvg], y=0.85,
                s="{}".format(unique_percentage[ind_max_mean_acc_hvg]), bbox=dict(facecolor='white', alpha=1),
                verticalalignment='center')
        ax.axvline(x=unique_percentage[ind_max_mean_acc_hig], ymin=0, ymax=1, color='darkorange', ls='--', lw=0.5)
        ax.text(x=unique_percentage[ind_max_mean_acc_hig], y=0.85,
                s="{}".format(unique_percentage[ind_max_mean_acc_hig]), bbox=dict(facecolor='white', alpha=1),
                verticalalignment='center')
        # Add Statistic result
        ax.text(x=70, y=0.75, s="statistics = {:.2f}\np-value = {:.2e}".format(stat, p_value))
        # y-limits
        ax.set_ylim([0.3, 0.9])
        ax.set_ylabel('mean accuracy', fontsize=fontsize_labels)
        ax.set_xlabel('gene percentage to keep', fontsize=fontsize_labels)
        ax.set_title('Leiden resolution {}'.format(num_unique_resolutions[ind]), fontsize=fontsize_title)
        # Increase ticks label size
        ax.tick_params(axis='both', which='major', labelsize=fontsize_ticks)
        ax.legend(frameon=False)
        sns.despine()
        plt.tight_layout()

        fig.savefig(os.path.join(
            save_folder,
            'MeanAccOverallSeeds_{}_Resolution_{}.{}'.format(default_title, num_unique_resolutions[ind], fileformat)))
        plt.close(fig=fig)

        counter += num_unique_percentage


def plot_genepercentagetokeep_vs_overallaccoverallresolutions_method(
        dict_res_method, method, save_folder, default_title):
    """
    Mean acc = acc overall seeds
    Overall Mean acc = Mean acc over all resolutions

    Parameters
    ----------
    dict_res_method
    method
    save_folder
    default_title

    Returns
    -------

    """
    if method == 'HVG':
        color = 'darkblue'
    else:
        color = 'darkorange'

    x_values = np.unique(dict_res_method['gene percentage to keep'])
    # Get max overall_acc_overall_resolutions value
    ind_max_mean_acc_method = np.nanargmax(dict_res_method['mean_overall_acc_overall_resolutions'])

    fig, ax = plt.subplots()
    ax.grid(linewidth=0.5)
    ax.errorbar(x_values,
                dict_res_method['mean_overall_acc_overall_resolutions'],
                yerr=dict_res_method['std_overall_acc_overall_resolutions'],
                label=method, fmt='-o', color=color, alpha=0.7)
    # Mark highest overall acc for each method with vertical line
    ax.axvline(x=x_values[ind_max_mean_acc_method], ymin=0, ymax=1, color='darkblue', ls='--', lw=0.5)
    plt.xticks(list(plt.xticks()[0]) + [x_values[ind_max_mean_acc_method], x_values[ind_max_mean_acc_method]])
    ax.set_ylim([0.45, 0.9])
    ax.set_xlim([0, 100.5])
    ax.set_xlabel('gene percentage to keep', fontsize=fontsize_labels)
    ax.set_ylabel('overall mean accuracy', fontsize=fontsize_labels)
    # Increase ticks label size
    ax.tick_params(axis='both', which='major', labelsize=fontsize_ticks)
    ax.legend()
    sns.despine(fig=fig, ax=ax)
    plt.tight_layout()

    fig.savefig(os.path.join(
        save_folder,
        'OverallMeanAccOverallResolutionsOverallseeds_{}.{}'.format(default_title, fileformat)))
    plt.close(fig=fig)


def plot_davies_bouldin_score(dict_res, save_folder, default_title):
    num_unique_resolutuions = np.unique(dict_res['resolution'])
    fig, ax = plt.subplots()
    ax.errorbar(num_unique_resolutuions, dict_res['mean_Davies Bouldin Score'],
                yerr=dict_res['std_Davies Bouldin Score'], fmt='-o')
    ax.set_xlabel('Leiden resolution', fontsize=fontsize_labels)
    ax.set_ylabel('Davies Bouldin index', fontsize=fontsize_labels)
    # Increase ticks label size
    ax.tick_params(axis='both', which='major', labelsize=fontsize_ticks)
    sns.despine()
    plt.tight_layout()

    fig.savefig(os.path.join(
        save_folder,
        'DaviesBouldinScore_{}.{}'.format(default_title, fileformat)))
    plt.close(fig=fig)


def plot_calinski_harabasz_index(dict_res, save_folder, default_title):
    """
    Higher value of CH index means the clusters are dense and well separated, although there is no
    “acceptable” cut-off value. We need to choose that solution which gives a peak or at least an
    abrupt elbow on the line plot of CH indices. On the other hand, if the line is smooth
    (horizontal or ascending or descending) then there is no such reason to prefer one solution over others.

    Parameters
    ----------
    dict_res
    save_folder
    default_title

    Returns
    -------

    """
    num_unique_resolutuions = np.unique(dict_res['resolution'])
    fig, ax = plt.subplots()
    ax.errorbar(num_unique_resolutuions, dict_res['mean_calinski_haravasz_index'],
                yerr=dict_res['std_calinski_haravasz_index'], fmt='-o')
    ax.set_xlabel('Leiden resolution', fontsize=fontsize_labels)
    ax.set_ylabel('Calinski-Harabasz Index', fontsize=fontsize_labels)
    # Increase ticks label size
    ax.tick_params(axis='both', which='major', labelsize=fontsize_ticks)
    sns.despine()
    plt.tight_layout()

    fig.savefig(os.path.join(
        save_folder,
        'CalinskiHarabaszIndex_{}.{}'.format(default_title, fileformat)))
    plt.close(fig=fig)


def plot_genepercentagetokeep_vs_overallaccoverallresolutions(dict_res_hvg, dict_res_hig, save_folder, default_title):
    """
    Mean acc = acc overall seeds
    Overall Mean acc = Mean acc over all resolutions

    Parameters
    ----------
    dict_res_hvg
    dict_res_hig
    save_folder
    default_title

    Returns
    -------

    """
    x_values = np.unique(dict_res_hvg['gene percentage to keep'])
    # Get max overall_acc_overall_resolutions value
    ind_max_mean_acc_hvg = np.nanargmax(dict_res_hvg['mean_overall_acc_overall_resolutions'])
    ind_max_mean_acc_hig = np.nanargmax(dict_res_hig['mean_overall_acc_overall_resolutions'])

    # Apply stats test to test whether the distribution differ significantly
    stat, p_value = mannwhitneyu(dict_res_hig['mean_overall_acc_overall_resolutions'],
                                 dict_res_hvg['mean_overall_acc_overall_resolutions'],
                                 alternative='greater')

    fig, ax = plt.subplots()
    ax.grid(linewidth=0.5)
    ax.errorbar(x_values,
                dict_res_hvg['mean_overall_acc_overall_resolutions'],
                yerr=dict_res_hvg['std_overall_acc_overall_resolutions'],
                label='HVG', fmt='-o', color='darkblue', alpha=0.7)
    ax.errorbar(x_values,
                dict_res_hig['mean_overall_acc_overall_resolutions'],
                yerr=dict_res_hig['std_overall_acc_overall_resolutions'],
                label='HIG', fmt='-o', color='darkorange', alpha=0.7)
    # Mark highest overall acc for each method with vertical line
    ax.axvline(x=x_values[ind_max_mean_acc_hvg], ymin=0, ymax=1, color='darkblue', ls='--', lw=0.5)
    ax.axvline(x=x_values[ind_max_mean_acc_hig], ymin=0, ymax=1, color='darkorange', ls='--', lw=0.5)
    plt.xticks(list(plt.xticks()[0]) + [x_values[ind_max_mean_acc_hvg], x_values[ind_max_mean_acc_hig]])
    # Add Statistic result
    ax.text(x=70, y=0.80, s="statistics = {:.2f}\np-value = {:.2e}".format(stat, p_value), fontsize=10)
    ax.set_ylim([0.45, 0.9])
    ax.set_xlim([0, 100.5])
    ax.set_xlabel('gene percentage to keep', fontsize=fontsize_labels)
    ax.set_ylabel('overall mean accuracy', fontsize=fontsize_labels)
    # Increase ticks label size
    ax.tick_params(axis='both', which='major', labelsize=fontsize_ticks)
    ax.legend()
    sns.despine(fig=fig, ax=ax)
    plt.tight_layout()

    fig.savefig(os.path.join(
        save_folder,
        'OverallMeanAccOverallResolutionsOverallseeds_{}.{}'.format(default_title, fileformat)))
    plt.close(fig=fig)

    return x_values[ind_max_mean_acc_hvg], x_values[ind_max_mean_acc_hig]


def plot_resolution_vs_overallacc_optimalgenekeep(dict_res_hvg, dict_res_hig, save_folder, default_title):
    """
    Mean acc = acc overall seeds
    Overall Mean acc = Mean acc over all resolutions

    Parameters
    ----------
    dict_res_hvg
    dict_res_hig
    save_folder
    default_title

    Returns
    -------

    """
    x_values = np.unique(dict_res_hvg['resolution'])
    # Get max overall_acc_overall_resolutions value
    ind_max_mean_acc_hvg = np.nanargmax(dict_res_hvg['mean_acc_overall_seeds_optimal_geneperkeep'])
    ind_max_mean_acc_hig = np.nanargmax(dict_res_hig['mean_acc_overall_seeds_optimal_geneperkeep'])

    # Apply stats test to test whether the distribution differ significantly
    stat, p_value = mannwhitneyu(dict_res_hig['mean_acc_overall_seeds_optimal_geneperkeep'],
                                 dict_res_hvg['mean_acc_overall_seeds_optimal_geneperkeep'],
                                 alternative='greater')

    fig, ax = plt.subplots()
    ax.grid(linewidth=0.5)
    ax.errorbar(x_values,
                dict_res_hvg['mean_acc_overall_seeds_optimal_geneperkeep'],
                yerr=dict_res_hvg['std_acc_overall_seeds_optimal_geneperkeep'],
                label='HVG', fmt='-o', color='darkblue', alpha=0.7)
    ax.errorbar(x_values,
                dict_res_hig['mean_acc_overall_seeds_optimal_geneperkeep'],
                yerr=dict_res_hig['std_acc_overall_seeds_optimal_geneperkeep'],
                label='HIG', fmt='-o', color='darkorange', alpha=0.7)
    # Mark highest overall acc for each method with vertical line
    ax.axvline(x=x_values[ind_max_mean_acc_hvg], ymin=0, ymax=1, color='darkblue', ls='--', lw=0.5)
    ax.axvline(x=x_values[ind_max_mean_acc_hig], ymin=0, ymax=1, color='darkorange', ls='--', lw=0.5)
    plt.xticks(list(plt.xticks()[0]) + [x_values[ind_max_mean_acc_hvg], x_values[ind_max_mean_acc_hig]])
    # Add Statistic result
    ax.text(x=0.35, y=0.45, s="statistics = {:.2f}\np-value = {:.2e}".format(stat, p_value), fontsize=10)
    ax.set_ylim([0.40, 0.85])
    ax.set_xlim([0, 1.1])
    ax.set_xlabel('resolution', fontsize=fontsize_labels)
    ax.set_ylabel('mean accuracy', fontsize=fontsize_labels)
    # Increase ticks label size
    ax.tick_params(axis='both', which='major', labelsize=fontsize_ticks)
    ax.legend()
    sns.despine(fig=fig, ax=ax)
    plt.tight_layout()

    fig.savefig(os.path.join(
        save_folder,
        'MeanAcc_OptimalGenePercKeep_{}.{}'.format(default_title, fileformat)))
    plt.close(fig=fig)


def plot_resolution_vs_overallmetric_optimalgenekeep(
        dict_res_hvg, dict_res_hig, save_folder, default_title, metric='acc', key=''):
    """
    Mean acc = acc overall seeds
    Overall Mean acc = Mean acc over all resolutions

    Parameters
    ----------
    dict_res_hvg
    dict_res_hig
    save_folder
    default_title
    metric

    Returns
    -------

    """
    obs_name = 'mean_{}_overall_seeds_optimal_geneperkeep_{}'.format(metric, key)
    std_obs_name = 'std_{}_overall_seeds_optimal_geneperkeep_{}'.format(metric, key)

    x_values = np.unique(dict_res_hvg['resolution'])
    # Get max overall_acc_overall_resolutions value
    ind_max_mean_acc_hvg = np.nanargmax(dict_res_hvg[obs_name])
    ind_max_mean_acc_hig = np.nanargmax(dict_res_hig[obs_name])

    # Apply stats test to test whether the distribution differ significantly
    stat, p_value = mannwhitneyu(dict_res_hig[obs_name],
                                 dict_res_hvg[obs_name],
                                 alternative='greater')

    fig, ax = plt.subplots()
    ax.grid(linewidth=0.5)
    ax.errorbar(x_values,
                dict_res_hvg[obs_name],
                yerr=dict_res_hvg[std_obs_name],
                label='HVG', fmt='-o', color='darkblue', alpha=0.7)
    ax.errorbar(x_values,
                dict_res_hig[obs_name],
                yerr=dict_res_hig[std_obs_name],
                label='HIG', fmt='-o', color='darkorange', alpha=0.7)
    # Mark highest overall acc for each method with vertical line
    ax.axvline(x=x_values[ind_max_mean_acc_hvg], ymin=0, ymax=1, color='darkblue', ls='--', lw=0.5)
    ax.axvline(x=x_values[ind_max_mean_acc_hig], ymin=0, ymax=1, color='darkorange', ls='--', lw=0.5)
    plt.xticks(list(plt.xticks()[0]) + [x_values[ind_max_mean_acc_hvg], x_values[ind_max_mean_acc_hig]])
    # Add Statistic result
    ax.text(x=0.35, y=0.45, s="statistics = {:.2f}\np-value = {:.2e}".format(stat, p_value), fontsize=10)
    # ax.set_ylim([0.40, 0.85])
    ax.set_xlim([0, 1.1])
    ax.set_xlabel('resolution', fontsize=fontsize_labels)
    ax.set_ylabel('mean {}'.format(metric), fontsize=fontsize_labels)
    # Increase ticks label size
    ax.tick_params(axis='both', which='major', labelsize=fontsize_ticks)
    ax.legend()
    sns.despine(fig=fig, ax=ax)
    plt.tight_layout()

    fig.savefig(os.path.join(
        save_folder,
        'Mean{}_OptimalGenePercKeep_{}.{}'.format(metric, default_title, fileformat)))
    plt.close(fig=fig)


def plot_genepercentagetokeep_vs_overallmetricoverallresolutions(
        dict_res_hvg, dict_res_hig, save_folder, default_title, metric='acc', alternative='greater', y_text=0.80):
    """
    Mean acc = acc overall seeds
    Overall Mean acc = Mean acc over all resolutions

    Parameters
    ----------
    dict_res_hvg
    dict_res_hig
    save_folder
    default_title
    metric
    alternative
    y_text

    Returns
    -------

    """
    x_values = np.unique(dict_res_hvg['gene percentage to keep'])
    obs_name = 'mean_overall_{}_overall_resolutions'.format(metric)
    std_obs_name = 'std_overall_{}_overall_resolutions'.format(metric)
    # Get optimal overall_acc_overall_resolutions value
    if alternative == 'greater':
        ind_mean_acc_hvg = np.nanargmax(dict_res_hvg[obs_name])
        ind_mean_acc_hig = np.nanargmax(dict_res_hig[obs_name])
    else:
        ind_mean_acc_hvg = np.nanargmin(dict_res_hvg[obs_name])
        ind_mean_acc_hig = np.nanargmin(dict_res_hig[obs_name])

    # Apply stats test to test whether the distribution differ significantly
    stat, p_value = mannwhitneyu(dict_res_hig[obs_name], dict_res_hvg[obs_name], alternative=alternative)

    fig, ax = plt.subplots()
    ax.grid(linewidth=0.5)
    ax.errorbar(x_values, dict_res_hvg[obs_name],
                yerr=dict_res_hvg[std_obs_name],
                label='HVG', fmt='-o', color='darkblue', alpha=0.7)
    ax.errorbar(x_values, dict_res_hig[obs_name],
                yerr=dict_res_hig[std_obs_name],
                label='HIG', fmt='-o', color='darkorange', alpha=0.7)
    # Mark highest overall acc for each method with vertical line
    ax.axvline(x=x_values[ind_mean_acc_hvg], ymin=0, ymax=1, color='darkblue', ls='--', lw=0.5)
    ax.axvline(x=x_values[ind_mean_acc_hig], ymin=0, ymax=1, color='darkorange', ls='--', lw=0.5)
    plt.xticks(list(plt.xticks()[0]) + [x_values[ind_mean_acc_hvg], x_values[ind_mean_acc_hig]])
    # Add Statistic result
    ax.text(x=70, y=y_text, s="statistics = {:.2f}\np-value = {:.2e}".format(stat, p_value), fontsize=10)
    ax.set_xlim([0, 100.5])
    ax.set_xlabel('gene percentage to keep', fontsize=fontsize_labels)
    # 'Mean {} overall resolutions overall seeds'.format(metric)
    ax.set_ylabel('overall mean {}'.format(metric), fontsize=fontsize_labels)
    # Increase ticks label size
    ax.tick_params(axis='both', which='major', labelsize=fontsize_ticks)
    ax.legend()
    sns.despine(fig=fig, ax=ax)
    plt.tight_layout()

    fig.savefig(os.path.join(save_folder, 'OverallMean{}OverallResolutionsOverallSeeds_{}.{}'.format(
        metric, default_title, fileformat)))
    plt.close(fig=fig)

    return x_values[ind_mean_acc_hvg], x_values[ind_mean_acc_hig]


def plot_resolution_vs_overallmetric_overallseeds_optimalgenekeep(
        dict_res_hvg, dict_res_hig, save_folder, default_title, metric='acc', key='',
        alternative='greater', y_text=0.45):
    """
    Mean acc = acc overall seeds
    Overall Mean acc = Mean acc over all resolutions

    Parameters
    ----------
    dict_res_hvg
    dict_res_hig
    save_folder
    default_title
    metric
    key
    alternative
    y_text

    Returns
    -------

    """
    x_values = np.unique(dict_res_hvg['resolution'])
    obs_name = 'mean_{}_overall_seeds_optimal_geneperkeep_{}'.format(metric, key)
    std_obs_name = 'std_{}_overall_seeds_optimal_geneperkeep_{}'.format(metric, key)

    # Get optimal overall_acc_overall_resolutions value
    if alternative == 'greater':
        ind_mean_acc_hvg = np.nanargmax(dict_res_hvg[obs_name])
        ind_mean_acc_hig = np.nanargmax(dict_res_hig[obs_name])
    else:
        ind_mean_acc_hvg = np.nanargmin(dict_res_hvg[obs_name])
        ind_mean_acc_hig = np.nanargmin(dict_res_hig[obs_name])

    # Apply stats test to test whether the distribution differ significantly
    stat, p_value = mannwhitneyu(dict_res_hig[obs_name], dict_res_hvg[obs_name], alternative=alternative)

    fig, ax = plt.subplots()
    ax.grid(linewidth=0.5)
    ax.errorbar(x_values, dict_res_hvg[obs_name], yerr=dict_res_hvg[std_obs_name],
                label='HVG', fmt='-o', color='darkblue', alpha=0.7)
    ax.errorbar(x_values, dict_res_hig[obs_name], yerr=dict_res_hig[std_obs_name],
                label='HIG', fmt='-o', color='darkorange', alpha=0.7)
    # Mark highest overall acc for each method with vertical line
    ax.axvline(x=x_values[ind_mean_acc_hvg], ymin=0, ymax=5, color='darkblue', ls='--', lw=0.5)
    ax.axvline(x=x_values[ind_mean_acc_hig], ymin=0, ymax=5, color='darkorange', ls='--', lw=0.5)
    plt.xticks(list(plt.xticks()[0]) + [x_values[ind_mean_acc_hvg], x_values[ind_mean_acc_hig]])
    # Add Statistic result
    ax.text(x=0.8, y=y_text, s="statistics = {:.2f}\np-value = {:.2e}".format(stat, p_value), fontsize=10)
    ax.set_xlim([0, 1.1])
    ax.set_xlabel('resolution', fontsize=fontsize_labels)
    ax.set_ylabel('overall mean {}'.format(metric), fontsize=fontsize_labels)
    # Increase ticks label size
    ax.tick_params(axis='both', which='major', labelsize=fontsize_ticks)
    ax.legend()
    sns.despine(fig=fig, ax=ax)
    plt.tight_layout()

    fig.savefig(os.path.join(
        save_folder,
        'OverallMean{}_OptimalGenePercKeep_{}.{}'.format(metric, default_title, fileformat)))
    plt.close(fig=fig)


def plot_probability_gene_to_be_selected(probability, genes, method, save_folder, num_ticks=1000):
    if method == 'HVG':
        color = 'darkblue'
    else:
        color = 'darkorange'

    xticks = np.arange(0, len(list(genes)), num_ticks)
    ind_xticks = list(xticks)

    fig, ax = plt.subplots()
    ax.grid(linewidth=0.5)
    ax.plot(np.arange(0, len(list(genes))), probability, color=color)
    ax.set_xticks(xticks)
    ax.set_xticklabels(np.asarray(list(genes))[ind_xticks])
    ax.set_ylim(0, 1.02)
    ax.set_ylabel('Probability', fontsize=fontsize_labels)
    ax.set_xlabel('Genes selected by {}'.format(method), fontsize=fontsize_labels)
    # Increase ticks label size
    ax.tick_params(axis='both', which='major', labelsize=fontsize_ticks)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    sns.despine(ax=ax)
    plt.tight_layout()

    plt.savefig(os.path.join(save_folder, 'Genes selected by {}.{}'.format(method, fileformat)))
    plt.close()


def calculate_bayesfactor(dict_res_hvg, dict_res_hig):
    unique_percentage = np.unique(dict_res_hvg['gene percentage to keep'])
    num_unique_percentage = len(unique_percentage)
    num_unique_resolutions = np.unique(dict_res_hvg['resolution'])

    dollar = base.__dict__["$"]
    at = base.__dict__["@"]

    # per cluster
    bf = dict.fromkeys(list(num_unique_resolutions))
    counter = 0
    for ind in range(0, len(num_unique_resolutions)):
        # df = pd.DataFrame.from_dict({'y': [0.80, 0.60, 0.75, 0.72, 0.77], 'x': [0.70, 0.80, 0.60, 0.74, 0.73]})
        df = pd.DataFrame.from_dict(
            {'HVG': dict_res_hvg['mean_overall_acc'][counter:counter + num_unique_percentage],
             'HIG': dict_res_hig['mean_overall_acc'][counter:counter + num_unique_percentage]})
        df_r = pandasToR(df=df)

        # compute the Bayes factor
        bf_temp = BayesFactor.ttestBF(y=dollar(df_r, 'HVG'), x=dollar(df_r, 'HIG'),
                                      nullInterval=np.asarray([0, 100]))
        bf_temp_values = dollar(BayesFactor.extractBF(bf_temp), 'bf')
        # Calculate BF score
        bf_score = bf_temp_values[0] / bf_temp_values[1]

        counter += num_unique_percentage

        bf[num_unique_resolutions[ind]] = bf_score

    return bf


def calculate_mannwhitneyu(dict_res_hvg, dict_res_hig):
    """
    Values smaller than 0.05 reject Null Hypothesis -> distributions are the "same"

    Parameters
    ----------
    dict_res_hvg
    dict_res_hig

    Returns
    -------

    """
    unique_percentage = np.unique(dict_res_hvg['gene percentage to keep'])
    num_unique_percentage = len(unique_percentage)
    num_unique_resolutions = np.unique(dict_res_hvg['resolution'])

    # per cluster
    score = dict.fromkeys(list(num_unique_resolutions))
    counter = 0
    for ind in range(0, len(num_unique_resolutions)):
        stat, p_value = mannwhitneyu(dict_res_hig['mean_overall_acc'][counter:counter + num_unique_percentage],
                                     dict_res_hvg['mean_overall_acc'][counter:counter + num_unique_percentage],
                                     alternative='greater')

        counter += num_unique_percentage

        score[num_unique_resolutions[ind]] = p_value

    return score


def get_mean_std_acc_ch_index(dict_res):
    mean_ch_index = []
    st_ch_index = []

    num_unique_resolutuions = np.unique(dict_res['resolution'])
    for resolution_val in num_unique_resolutuions:
        ind_clusters = np.where(np.asarray(dict_res['resolution']) == resolution_val)[0]
        # Get Mean and std  of Davies Boulding score over all gene percentage to keep
        nclusters_db_score = list(map(dict_res['Calinski-Harabasz Index'].__getitem__, ind_clusters))
        mean_ch_index.append(np.nanmean(nclusters_db_score))
        st_ch_index.append(np.nanstd(nclusters_db_score))

    dict_res['mean_calinski_haravasz_index'] = mean_ch_index
    dict_res['std_calinski_haravasz_index'] = st_ch_index

    return dict_res


def get_overall_acc_overall_resolutions(dict_res):
    mean_overall_acc_overall_resolutions = []
    st_overall_acc_overall_resolutions = []

    df_res = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in dict_res.items()]))

    df_res = df_res[['seed', 'resolution', 'gene percentage to keep', 'acc']]

    # for res in df_res['resolution'].unique():
    for gene_to_kepp_per in df_res['gene percentage to keep'].unique():
        # Calculate acc of overall_resolutions per gene percentage to keep
        mean_overall_acc_overall_resolutions.append(
            df_res.loc[(df_res['gene percentage to keep'] == gene_to_kepp_per), 'acc'].mean())
        st_overall_acc_overall_resolutions.append(
            df_res.loc[(df_res['gene percentage to keep'] == gene_to_kepp_per), 'acc'].std())

    dict_res['mean_overall_acc_overall_resolutions'] = mean_overall_acc_overall_resolutions
    dict_res['std_overall_acc_overall_resolutions'] = st_overall_acc_overall_resolutions

    return dict_res


def get_overall_metric_overall_resolutions(dict_res, metric='acc'):
    """

    Parameters
    ----------
    dict_res
    metric

    Returns
    -------

    """
    mean_overall_acc_overall_resolutions = []
    st_overall_acc_overall_resolutions = []

    df_res = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in dict_res.items()]))

    df_res = df_res[['seed', 'resolution', 'gene percentage to keep', metric]]

    # for res in df_res['resolution'].unique():
    for gene_to_kepp_per in df_res['gene percentage to keep'].unique():
        # Calculate acc of overall_resolutions per gene percentage to keep
        mean_overall_acc_overall_resolutions.append(
            df_res.loc[(df_res['gene percentage to keep'] == gene_to_kepp_per), metric].mean())
        st_overall_acc_overall_resolutions.append(
            df_res.loc[(df_res['gene percentage to keep'] == gene_to_kepp_per), metric].std())

    dict_res['mean_overall_{}_overall_resolutions'.format(metric)] = mean_overall_acc_overall_resolutions
    dict_res['std_overall_{}_overall_resolutions'.format(metric)] = st_overall_acc_overall_resolutions

    return dict_res


def get_mean_acc_overall_seeds_for_optimal_genepercentagetokeep_value(dict_res, opt_geneperckeep):
    mean_acc_overall_seeds = []
    std_acc_overall_seeds = []

    df_res = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in dict_res.items()]))

    df_res = df_res[['seed', 'resolution', 'gene percentage to keep', 'acc']]
    df_res = df_res.loc[df_res['gene percentage to keep'] == opt_geneperckeep, :]

    for res in df_res['resolution'].unique():
        mean_acc_overall_seeds.append(df_res.loc[df_res['resolution'] == res, 'acc'].mean())
        std_acc_overall_seeds.append(df_res.loc[df_res['resolution'] == res, 'acc'].std())

    dict_res['mean_acc_overall_seeds_optimal_geneperkeep'] = mean_acc_overall_seeds
    dict_res['std_acc_overall_seeds_optimal_geneperkeep'] = std_acc_overall_seeds

    return dict_res


def get_mean_calinskiharabszscore_overall_seeds_for_optimal_genepercentagetokeep(dict_res, opt_geneperckeep):
    mean_acc_overall_seeds = []
    std_acc_overall_seeds = []

    df_res = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in dict_res.items()]))

    df_res = df_res[['seed', 'resolution', 'gene percentage to keep', 'Calinski-Harabasz Index']]
    df_res = df_res.loc[df_res['gene percentage to keep'] == opt_geneperckeep, :]

    for res in df_res['resolution'].unique():
        mean_acc_overall_seeds.append(df_res.loc[df_res['resolution'] == res, 'Calinski-Harabasz Index'].mean())
        std_acc_overall_seeds.append(df_res.loc[df_res['resolution'] == res, 'Calinski-Harabasz Index'].std())

    dict_res['mean_CHI_overall_seeds_optimal_geneperkeep'] = mean_acc_overall_seeds
    dict_res['std_CHI_overall_seeds_optimal_geneperkeep'] = std_acc_overall_seeds

    return dict_res


def get_mean_metric_overall_seeds_for_optimal_genepercentagetokeep(dict_res, opt_geneperckeep, metric='acc', key=''):
    mean_acc_overall_seeds = []
    std_acc_overall_seeds = []

    df_res = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in dict_res.items()]))

    df_res = df_res[['seed', 'resolution', 'gene percentage to keep', metric]]
    df_res = df_res.loc[df_res['gene percentage to keep'] == opt_geneperckeep, :]

    for res in df_res['resolution'].unique():
        mean_acc_overall_seeds.append(df_res.loc[df_res['resolution'] == res, metric].mean())
        std_acc_overall_seeds.append(df_res.loc[df_res['resolution'] == res, metric].std())

    dict_res['mean_{}_overall_seeds_optimal_geneperkeep_{}'.format(metric, key)] = mean_acc_overall_seeds
    dict_res['std_{}_overall_seeds_optimal_geneperkeep_{}'.format(metric, key)] = std_acc_overall_seeds

    return dict_res


def get_probability_gene_topnpercent(dict_res, genes, metric):
    # Calculate probability that a genes falls into top 10% of feature selection method (HVG & HIG)
    df_res = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in dict_res.items()]))
    df_res = df_res[['gene percentage to keep', 'seed', 'resolution', 'index genes']]
    # overall seeds and resolutions for optimal gene percentage to keep value
    df_res = df_res.loc[df_res['gene percentage to keep'] == 10, :]

    gene_names = np.asarray(genes)
    genenames = []
    for res in df_res['resolution'].unique():
        for seed in df_res['seed'].unique():
            if metric == 'HVG':
                genenames.extend(
                    list(gene_names[list(
                        df_res.loc[(df_res['resolution'] == res) & (
                                df_res['seed'] == seed), 'index genes'].values[0].values)]))
            else:
                genenames.extend(
                    list(gene_names[list(
                        df_res.loc[(df_res['resolution'] == res) & (
                                df_res['seed'] == seed), 'index genes'].values[0])]))

    # max occurence
    max_occurence = df_res.shape[0]
    # Count frequencies
    count = Counter(genenames)
    # sort frequencies
    count.most_common()
    # norm frequencies by number of parameter settings
    for val in count.keys():
        count[val] = count[val] / max_occurence
    count.most_common()
    unzipped_freq = list(zip(*count.most_common()))
    sorted_freq = list(count.values())
    sorted_freq.sort(reverse=True)

    return unzipped_freq, sorted_freq


if __name__ == '__main__':
    path_input = "/Volumes/CH__data/Projects/Eyerich_AG_projects/BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer/analysis/Molecular_subtypes/output/Leiden"
    file_hig = os.path.join(
        path_input, '2023-02-12', 'resolution_[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1. ]__featureselection_HIG',
        'Scores_GridSearch_9_minPercentage0.5_maxPercentage100.0_minResolution0.1_maxResolution1.0_.pkl')
    with open(file_hig, 'rb') as ff:
        gs_hig = pickle.load(ff)

    file_hvg = os.path.join(
        path_input, '2023-02-13', 'resolution_[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1. ]__featureselection_HVG',
        'Scores_GridSearch_9_minPercentage0.5_maxPercentage100.0_minResolution0.1_maxResolution1.0_.pkl')
    with open(file_hvg, 'rb') as ff:
        gs_hvg = pickle.load(ff)

    # file_hvg = os.path.join(
    #     path_input, '2023-11-25', 'resolution_[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1. ]__featureselection_std',
    #     'Scores_GridSearch_9_minPercentage0.5_maxPercentage100.0_minResolution0.1_maxResolution1.0_.pkl')
    # with open(file_hvg, 'rb') as ff:
    #     gs_std = pickle.load(ff)

    # Load gene names
    adata = sc.read(os.path.join(
        data_root, 'Molecular_subtypes', 'input', 'h5_files', 'LESION',
        'Lesion_RNAseq_20210720_patient_meta_data_v04__CH_20220209_TherapyResponse_v04_v220715.h5'))
    var_names = list(adata.var_names)
    del adata

    output_dir = os.path.join(path_input, 'comparing_edgeR_HIG_HVG', str(date.today()))
    os.makedirs(output_dir, exist_ok=True)
    highest_mean_resolution_acc_hvg = plot_errorbars(dict_res=gs_hvg, save_folder=output_dir, default_title='HVG')
    plot_davies_bouldin_score(dict_res=gs_hvg, save_folder=output_dir, default_title='HVG')

    gs_hvg = get_mean_std_acc_ch_index(dict_res=gs_hvg)
    plot_calinski_harabasz_index(dict_res=gs_hvg, save_folder=output_dir, default_title='HVG')

    # Identify highest mean acc (mean by seeds) across gene percentage to keep
    highest_mean_resolution_acc_hig = plot_errorbars(dict_res=gs_hig, save_folder=output_dir, default_title='HIG')
    plot_davies_bouldin_score(dict_res=gs_hig, save_folder=output_dir, default_title='HIG')

    gs_hig = get_mean_std_acc_ch_index(dict_res=gs_hig)
    plot_calinski_harabasz_index(dict_res=gs_hig, save_folder=output_dir, default_title='HIG')

    # Read out gene percentage to keep
    gene_per_keep_hvg = []
    for dict_item in highest_mean_resolution_acc_hvg.items():
        gene_per_keep_hvg.append(dict_item[1][1])
    counts_hvg = collections.Counter(gene_per_keep_hvg)

    gene_per_keep_hig = []
    for dict_item in highest_mean_resolution_acc_hig.items():
        gene_per_keep_hig.append(dict_item[1][1])
    counts_hig = collections.Counter(gene_per_keep_hig)

    # Take gene percentage to keep which occurred the most in counter
    # Draw plot gene percentage to keep vs mean overall acc comparing HVG and HIG per resolution
    plot_errorbars_HIG_vs_HVG(
        dict_res_hvg=gs_hvg, dict_res_hig=gs_hig, save_folder=output_dir, default_title='HVG_vs_HIG')
    # Apply Mann-WhitneyU test
    mannwhitneyu_pvals = calculate_mannwhitneyu(dict_res_hvg=gs_hvg, dict_res_hig=gs_hig)

    # Get mean acc overall seeds and resolutions
    gs_hig = get_overall_acc_overall_resolutions(dict_res=gs_hig)
    gs_hvg = get_overall_acc_overall_resolutions(dict_res=gs_hvg)

    optimal_genekeep_hvg, optimal_genekeep_hig = plot_genepercentagetokeep_vs_overallaccoverallresolutions(
        dict_res_hvg=gs_hvg, dict_res_hig=gs_hig, save_folder=output_dir, default_title='HVG_vs_HIG')
    plot_genepercentagetokeep_vs_overallaccoverallresolutions_method(
        dict_res_method=gs_hig, method='HIG', save_folder=output_dir, default_title='HIG')

    """ 
    Note: 
    - Optimal gene percentage to keep parameter is 6 for HVG and 7 for HIG
    - Mann WhitneyU test rejects null hypothesis -> Methods HVG and HIG differ significantly
    - HIG is significantly better than HVG -> use HIG for feature selection
    """

    # Find highest resolution for gene percentage to keep parameter 7:
    # x = resolution y = mean acc over all seeds
    gs_hvg = get_mean_acc_overall_seeds_for_optimal_genepercentagetokeep_value(
        dict_res=gs_hvg, opt_geneperckeep=optimal_genekeep_hvg)
    gs_hig = get_mean_acc_overall_seeds_for_optimal_genepercentagetokeep_value(
        dict_res=gs_hig, opt_geneperckeep=optimal_genekeep_hig)
    plot_resolution_vs_overallacc_optimalgenekeep(
        dict_res_hvg=gs_hvg, dict_res_hig=gs_hig, save_folder=output_dir, default_title='HVG_vs_HIG')

    # ================ Calinski Harabasz Index
    gs_hvg = get_overall_metric_overall_resolutions(dict_res=gs_hvg, metric='Calinski-Harabasz Index')
    gs_hig = get_overall_metric_overall_resolutions(dict_res=gs_hig, metric='Calinski-Harabasz Index')
    optimal_genekeep_hvg_chi, optimal_genekeep_hig_chi = plot_genepercentagetokeep_vs_overallmetricoverallresolutions(
        dict_res_hvg=gs_hvg, dict_res_hig=gs_hig, save_folder=output_dir, default_title='HVG_vs_HIG',
        metric='Calinski-Harabasz Index', alternative='greater', y_text=2)

    # Using optimal Gene percentage to keep derived from Calinski-Harabasz Index
    gs_hvg = get_mean_metric_overall_seeds_for_optimal_genepercentagetokeep(
        dict_res=gs_hvg, opt_geneperckeep=optimal_genekeep_hvg_chi, metric='Calinski-Harabasz Index', key='best_CHI')
    gs_hig = get_mean_metric_overall_seeds_for_optimal_genepercentagetokeep(
        dict_res=gs_hig, opt_geneperckeep=optimal_genekeep_hig_chi, metric='Calinski-Harabasz Index', key='best_CHI')

    plot_resolution_vs_overallmetric_optimalgenekeep(
        dict_res_hvg=gs_hvg, dict_res_hig=gs_hig, save_folder=output_dir, default_title='HVG_vs_HIG',
        metric='Calinski-Harabasz Index', key='best_CHI')

    plot_resolution_vs_overallmetric_overallseeds_optimalgenekeep(
        dict_res_hvg=gs_hvg, dict_res_hig=gs_hig, save_folder=output_dir, default_title='HVG_vs_HIG',
        metric='Calinski-Harabasz Index', alternative='greater', y_text=2, key='best_CHI')

    # Using optimal Gene percentage to keep derived from acc
    gs_hvg = get_mean_metric_overall_seeds_for_optimal_genepercentagetokeep(
        dict_res=gs_hvg, opt_geneperckeep=optimal_genekeep_hvg, metric='Calinski-Harabasz Index', key='best_acc')
    gs_hig = get_mean_metric_overall_seeds_for_optimal_genepercentagetokeep(
        dict_res=gs_hig, opt_geneperckeep=optimal_genekeep_hig, metric='Calinski-Harabasz Index', key='best_acc')

    plot_resolution_vs_overallmetric_optimalgenekeep(
        dict_res_hvg=gs_hvg, dict_res_hig=gs_hig, save_folder=output_dir, default_title='HVG_vs_HIG',
        metric='Calinski-Harabasz Index', key='best_acc')

    plot_resolution_vs_overallmetric_overallseeds_optimalgenekeep(
        dict_res_hvg=gs_hvg, dict_res_hig=gs_hig, save_folder=output_dir, default_title='HVG_vs_HIG',
        metric='Calinski-Harabasz Index', alternative='greater', y_text=2, key='best_acc')

    # ================ Davies Bouldin Score
    gs_hvg = get_overall_metric_overall_resolutions(dict_res=gs_hvg, metric='Davies Bouldin Score')
    gs_hig = get_overall_metric_overall_resolutions(dict_res=gs_hig, metric='Davies Bouldin Score')

    optimal_genekeep_hvg_dbs, optimal_genekeep_hig_dbs = plot_genepercentagetokeep_vs_overallmetricoverallresolutions(
        dict_res_hvg=gs_hvg, dict_res_hig=gs_hig, save_folder=output_dir, default_title='HVG_vs_HIG',
        metric='Davies Bouldin Score', alternative='less', y_text=3.4)

    gs_hvg = get_mean_metric_overall_seeds_for_optimal_genepercentagetokeep(
        dict_res=gs_hvg, opt_geneperckeep=optimal_genekeep_hvg_dbs, metric='Davies Bouldin Score', key='best_DBS')
    gs_hig = get_mean_metric_overall_seeds_for_optimal_genepercentagetokeep(
        dict_res=gs_hig, opt_geneperckeep=optimal_genekeep_hig_dbs, metric='Davies Bouldin Score', key='best_DBS')

    plot_resolution_vs_overallmetric_optimalgenekeep(
        dict_res_hvg=gs_hvg, dict_res_hig=gs_hig, save_folder=output_dir, default_title='HVG_vs_HIG',
        metric='Davies Bouldin Score', key='best_DBS')

    plot_resolution_vs_overallmetric_overallseeds_optimalgenekeep(
        dict_res_hvg=gs_hvg, dict_res_hig=gs_hig, save_folder=output_dir, default_title='HVG_vs_HIG',
        metric='Davies Bouldin Score', alternative='less', y_text=2.5, key='best_DBS')

    # Using optimal Genen percentage to keep of acc
    gs_hvg = get_mean_metric_overall_seeds_for_optimal_genepercentagetokeep(
        dict_res=gs_hvg, opt_geneperckeep=optimal_genekeep_hvg, metric='Davies Bouldin Score', key='best_DBS')
    gs_hig = get_mean_metric_overall_seeds_for_optimal_genepercentagetokeep(
        dict_res=gs_hig, opt_geneperckeep=optimal_genekeep_hig, metric='Davies Bouldin Score', key='best_DBS')

    plot_resolution_vs_overallmetric_optimalgenekeep(
        dict_res_hvg=gs_hvg, dict_res_hig=gs_hig, save_folder=output_dir, default_title='HVG_vs_HIG',
        metric='Davies Bouldin Score', key='best_DBS')

    plot_resolution_vs_overallmetric_overallseeds_optimalgenekeep(
        dict_res_hvg=gs_hvg, dict_res_hig=gs_hig, save_folder=output_dir, default_title='HVG_vs_HIG',
        metric='Davies Bouldin Score', alternative='less', y_text=2.5, key='best_DBS')

    # ================ Comparing Feature selection methods HVG vs HIG using likelihood a gene falls into 10%
    # Calculate probability that a genes falls into top 10% of feature selection method (HVG & HIG)
    unzipped_freq_hvg, sorted_freq_hvg = get_probability_gene_topnpercent(
        dict_res=gs_hvg, genes=var_names, metric='HVG')
    unzipped_freq_hig, sorted_freq_hig = get_probability_gene_topnpercent(
        dict_res=gs_hig, genes=var_names, metric='HIG')

    plot_probability_gene_to_be_selected(
        probability=sorted_freq_hvg, genes=unzipped_freq_hvg[0], method='HVG', save_folder=output_dir, num_ticks=1000)

    plot_probability_gene_to_be_selected(
        probability=sorted_freq_hig, genes=unzipped_freq_hig[0], method='HIG', save_folder=output_dir, num_ticks=500)

    # save dataframe
    dict_selected_genes = {'HVG': unzipped_freq_hvg[0], 'HIG': unzipped_freq_hig[0]}
    unique_hvg = list(set(dict_selected_genes['HVG']) - set(dict_selected_genes['HIG']))
    unique_hig = list(set(dict_selected_genes['HIG']) - set(dict_selected_genes['HVG']))
    dict_selected_genes['unique_HIG'] = unique_hig
    dict_selected_genes['unique_HVG'] = unique_hvg
    shared = np.intersect1d(dict_selected_genes['HVG'], dict_selected_genes['HIG'])
    dict_selected_genes['Shared_genes'] = shared

    df_selected_genes = pd.DataFrame.from_dict(dict([(k, pd.Series(v)) for k, v in dict_selected_genes.items()]))

    df_selected_genes.to_excel(os.path.join(output_dir, 'Top10per_Selected_Genes.xlsx'))
