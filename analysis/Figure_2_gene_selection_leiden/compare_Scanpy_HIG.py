import numpy as np
import scanpy as sc
import os
from datetime import date
import pickle
import collections
import pandas as pd
import baycomp as bt
import itertools
from tabulate import tabulate
from collections import Counter
from scipy.stats import mannwhitneyu, wilcoxon
import platform

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import seaborn as sns

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


def run_baycomp(combinations: dict, name_models_dict: dict):
    # Calculate probability that classifier 1 is better than 2 and vice versa
    cct = dict()

    for comb in combinations.keys():
        cct[comb] = bt.two_on_single(
            x=np.asarray(combinations[comb][0]),
            y=np.asarray(combinations[comb][1]),
            rope=0,
            runs=1,  # same number as n_repeats used in RepeatedStratifiedKFold
            plot=False,
        )

    df_cct = pd.DataFrame.from_dict(cct)
    df_cct.index = ["BF_x", "BF_y"]
    # Get Bayes factor: calculate the difference between BF_x and BF_y
    df_cct.loc['bayes_factor', :] = df_cct.loc['BF_x', :] / df_cct.loc['BF_y', :]
    print(tabulate(df_cct.T, headers="keys", tablefmt="psql"))

    # Save as 2D matrix
    num_models = len(name_models_dict.keys())
    comparison_matrix = np.zeros((num_models, num_models))
    names_models = list(name_models_dict.keys())
    for ind, comb in enumerate(names_models):
        for j in range(ind + 1, num_models):
            bf = bt.two_on_single(
                x=np.asarray(name_models_dict[comb]),
                y=np.asarray(name_models_dict[names_models[j]]),
                rope=0,
                runs=10,  # same number as n_repeats used in RepeatedStratifiedKFold
                plot=False,
            )

            # Store Bayes factor in the comparison matrix (symmetrically)
            bayes_factor = bf[0] / bf[1]
            comparison_matrix[ind, j] = bayes_factor
            comparison_matrix[j, ind] = 1 / bayes_factor  # Inverse as it's symmetric

    df_bf = pd.DataFrame(data=comparison_matrix, columns=names_models, index=names_models)

    return df_cct, df_bf


def _pairwise(dict_comb: dict):
    # combination of keys
    key_comb = list(itertools.combinations(dict_comb.keys(), 2))

    combinations = {}
    for pair in key_comb:
        combinations["{} vs {}".format(pair[0], pair[1])] = (
            dict_comb[pair[0]],
            dict_comb[pair[1]],
        )

    return combinations


def plot_errorbars_HIG_vs_HVG_std(dict_res_hvg, dict_res_hig, dict_res_std, save_folder, default_title):
    """
    Mean acc = acc overall seeds
    Overall Mean acc = Mean acc over all resolutions

    Parameters
    ----------
    dict_res_hvg
    dict_res_hig
    dict_res_std
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
        ind_max_mean_acc_std = np.nanargmax(dict_res_std['mean_overall_acc'][counter:counter + num_unique_percentage])

        print('Highest acc for HRG at resolution {}: {:.2f}%'.format(
            num_unique_resolutions[ind],
            dict_res_hig['mean_overall_acc'][counter:counter + num_unique_percentage][ind_max_mean_acc_hig] * 100))
        print('Highest acc for SCANPY at resolution {}: {:.2f}%'.format(
            num_unique_resolutions[ind],
            dict_res_hvg['mean_overall_acc'][counter:counter + num_unique_percentage][ind_max_mean_acc_hvg] * 100))
        print('Highest acc for SD at resolution {}: {:.2f}%'.format(
            num_unique_resolutions[ind],
            dict_res_std['mean_overall_acc'][counter:counter + num_unique_percentage][ind_max_mean_acc_std] * 100))

        stat_p_value_hig_vs_hvg = wilcoxon(
            dict_res_hig['mean_overall_acc'][counter:counter + num_unique_percentage],
            dict_res_hvg['mean_overall_acc'][counter:counter + num_unique_percentage], alternative='greater')
        stat_p_value_hig_vs_std = wilcoxon(
            dict_res_hig['mean_overall_acc'][counter:counter + num_unique_percentage],
            dict_res_std['mean_overall_acc'][counter:counter + num_unique_percentage], alternative='greater')

        fig, ax = plt.subplots(figsize=(6, 6))
        ax.grid(linewidth=0.5)
        # HVG
        ax.errorbar(np.unique(dict_res_hvg['gene percentage to keep']),
                    dict_res_hvg['mean_overall_acc'][counter:counter + num_unique_percentage],
                    yerr=dict_res_hvg['std_overall_acc'][counter:counter + num_unique_percentage], fmt='-o',
                    color='darkblue', alpha=0.7, label='HVG')
        # HIG
        ax.errorbar(np.unique(dict_res_hig['gene percentage to keep']),
                    dict_res_hig['mean_overall_acc'][counter:counter + num_unique_percentage],
                    yerr=dict_res_hig['std_overall_acc'][counter:counter + num_unique_percentage], fmt='-o',
                    color='darkorange', alpha=0.7, label='HIG')
        # Std
        ax.errorbar(np.unique(dict_res_std['gene percentage to keep']),
                    dict_res_std['mean_overall_acc'][counter:counter + num_unique_percentage],
                    yerr=dict_res_std['std_overall_acc'][counter:counter + num_unique_percentage], fmt='-o',
                    color='darkred', alpha=0.7, label='std')

        # Mark highest overall acc for each method with vertical line
        # HVG
        ax.axvline(x=unique_percentage[ind_max_mean_acc_hvg], ymin=0, ymax=1, color='darkblue', ls='--', lw=0.5)
        ax.text(x=unique_percentage[ind_max_mean_acc_hvg], y=0.85,
                s="{}".format(unique_percentage[ind_max_mean_acc_hvg]), bbox=dict(facecolor='white', alpha=1),
                verticalalignment='center')
        # HIG
        ax.axvline(x=unique_percentage[ind_max_mean_acc_hig], ymin=0, ymax=1, color='darkorange', ls='--', lw=0.5)
        ax.text(x=unique_percentage[ind_max_mean_acc_hig], y=0.9,
                s="{}".format(unique_percentage[ind_max_mean_acc_hig]), bbox=dict(facecolor='white', alpha=1),
                verticalalignment='center')
        # Std
        ax.axvline(x=unique_percentage[ind_max_mean_acc_std], ymin=0, ymax=1, color='darkred', ls='--', lw=0.5)
        ax.text(x=unique_percentage[ind_max_mean_acc_std], y=0.85,
                s="{}".format(unique_percentage[ind_max_mean_acc_std]), bbox=dict(facecolor='white', alpha=1),
                verticalalignment='center')

        # Add Statistic result
        ax.text(x=70, y=0.75, s="p{} = {:.2e}\np{} = {:.2e}".format(
            r'$_{HIG\; vs\; HVG}$', stat_p_value_hig_vs_hvg[1],
            r'$_{HIG\; vs\; std}$', stat_p_value_hig_vs_std[1]),
                fontsize=10, horizontalalignment='center', verticalalignment='center')
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
            'MeanAccOverallSeeds_{}_Resolution_{}.{}'.format(
                default_title, num_unique_resolutions[ind], fileformat)))
        plt.close(fig=fig)

        counter += num_unique_percentage


def plot_genepercentagetokeep_vs_overallaccoverallresolutions(
        dict_res_hvg, dict_res_hig, dict_res_std, save_folder, default_title):
    """
    Mean acc = acc overall seeds
    Overall Mean acc = Mean acc over all resolutions

    Parameters
    ----------
    dict_res_hvg
    dict_res_hig
    dict_res_std
    save_folder
    default_title

    Returns
    -------

    """
    metric = 'acc'

    x_values = np.unique(dict_res_hvg['gene percentage to keep'])
    # Get max overall_acc_overall_resolutions value
    ind_max_mean_acc_hvg = np.nanargmax(dict_res_hvg['mean_overall_acc_overall_resolutions'])
    ind_max_mean_acc_hig = np.nanargmax(dict_res_hig['mean_overall_acc_overall_resolutions'])
    ind_max_mean_acc_std = np.nanargmax(dict_res_std['mean_overall_acc_overall_resolutions'])

    # Apply stats test to test whether the distribution differ significantly;
    # Wilcoxon signed-rank test is used for paired or related samples where each observation in one group
    # corresponds to a specific observation in the other group.
    result_hig_vs_hvg = wilcoxon(dict_res_hig['mean_overall_{}_overall_resolutions'.format(metric)],
                                 dict_res_hvg['mean_overall_{}_overall_resolutions'.format(metric)],
                                 alternative="greater")
    result_hig_vs_std = wilcoxon(dict_res_hig['mean_overall_{}_overall_resolutions'.format(metric)],
                                 dict_res_std['mean_overall_{}_overall_resolutions'.format(metric)],
                                 alternative="greater")

    print('Highest {} for HRG over all resolutions: {:.2f}% {} {:.2f}'.format(
        metric, dict_res_hig['mean_overall_acc_overall_resolutions'][ind_max_mean_acc_hig],
        r'$\pm$', dict_res_hig['std_overall_acc_overall_resolutions'][ind_max_mean_acc_hig]))
    print('Highest {} for SCANPY over all resolutions: {:.2f}% {} {:.2f}'.format(
        metric, dict_res_hvg['mean_overall_acc_overall_resolutions'][ind_max_mean_acc_hvg],
        r'$\pm$', dict_res_hvg['std_overall_acc_overall_resolutions'][ind_max_mean_acc_hvg]))
    print('Highest {} for SD over all resolutions: {:.2f}% {} {:.2f}'.format(
        metric, dict_res_std['mean_overall_acc_overall_resolutions'][ind_max_mean_acc_std],
        r'$\pm$', dict_res_std['std_overall_acc_overall_resolutions'][ind_max_mean_acc_std]))

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
    # Std
    ax.errorbar(x_values,
                dict_res_std['mean_overall_acc_overall_resolutions'],
                yerr=dict_res_std['std_overall_acc_overall_resolutions'],
                label='std', fmt='-o', color='darkred', alpha=0.7)
    # Mark highest overall acc for each method with vertical line
    ax.axvline(x=x_values[ind_max_mean_acc_hvg], ymin=0, ymax=1, color='darkblue', ls='--', lw=0.5)
    ax.axvline(x=x_values[ind_max_mean_acc_std], ymin=0, ymax=1, color='darkred', ls='--', lw=0.5)
    ax.axvline(x=x_values[ind_max_mean_acc_hig], ymin=0, ymax=1, color='darkorange', ls='--', lw=0.5)
    plt.xticks(list(plt.xticks()[0]) + [x_values[ind_max_mean_acc_hvg], x_values[ind_max_mean_acc_std],
                                        x_values[ind_max_mean_acc_hig]])
    # Add Statistic result
    ax.text(x=25, y=48, s="p{} = {:.2e}\np{} = {:.2e}".format(
        r'$_{HIG\; vs\; HVG}$', result_hig_vs_hvg[1], r'$_{HIG\; vs\; std}$', result_hig_vs_std[1]),
            fontsize=10, horizontalalignment='center', verticalalignment='center')
    ax.set_ylim([45, 90])
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


def plot_resolution_vs_overallacc_optimalgenekeep(
        dict_res_hvg, dict_res_hig, dict_res_std, save_folder, default_title, metric='acc', key=''):
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
    obs = 'mean_{}_overall_seeds_optimal_geneperkeep_{}'.format(metric, key)
    obs_std = 'std_{}_overall_seeds_optimal_geneperkeep_{}'.format(metric, key)
    x_values = np.unique(dict_res_hvg['resolution'])
    # Get max overall_acc_overall_resolutions value
    ind_max_mean_acc_hvg = np.nanargmax(dict_res_hvg[obs])
    ind_max_mean_acc_hig = np.nanargmax(dict_res_hig[obs])
    ind_max_mean_acc_std = np.nanargmax(dict_res_std[obs])

    # Apply stats test to test whether the distribution differ significantly
    stat_p_value_hrg_scanpy = wilcoxon(dict_res_hig[obs], dict_res_hvg[obs], alternative='greater')
    stat_p_value_hrg_sd = wilcoxon(dict_res_hig[obs], dict_res_std[obs], alternative='greater')

    fig, ax = plt.subplots()
    ax.grid(linewidth=0.5)
    ax.errorbar(x_values,
                dict_res_hvg[obs],
                yerr=dict_res_hvg[obs_std],
                label='HVG', fmt='-o', color='darkblue', alpha=0.7)
    ax.errorbar(x_values,
                dict_res_hig[obs],
                yerr=dict_res_hig[obs_std],
                label='HIG', fmt='-o', color='darkorange', alpha=0.7)
    # std
    ax.errorbar(x_values,
                dict_res_std[obs],
                yerr=dict_res_std[obs_std],
                label='std', fmt='-o', color='darkred', alpha=0.7)
    # Mark highest overall acc for each method with vertical line
    ax.axvline(x=x_values[ind_max_mean_acc_hvg], ymin=0, ymax=1, color='darkblue', ls='--', lw=0.5)
    ax.axvline(x=x_values[ind_max_mean_acc_std], ymin=0, ymax=1, color='darkred', ls='--', lw=0.5)
    ax.axvline(x=x_values[ind_max_mean_acc_hig], ymin=0, ymax=1, color='darkorange', ls='--', lw=0.5)
    plt.xticks(list(plt.xticks()[0]) + [x_values[ind_max_mean_acc_hvg], x_values[ind_max_mean_acc_std],
                                        x_values[ind_max_mean_acc_hig]])
    # Add Statistic result
    # Add Statistic result
    ax.text(x=0.35, y=0.45, s="p{} = {:.2e}\np{} = {:.2e}".format(
        r'$_{HIG\; vs\; HVG}$', stat_p_value_hrg_scanpy[1], r'$_{HIG\; vs\; std}$', stat_p_value_hrg_sd[1]),
            fontsize=10, horizontalalignment='center', verticalalignment='center')
    if metric == 'acc':
        ax.set_ylim([0, 100])
        ax.set_xlim([0, 1.5])
    else:
        ax.set_ylim([0.40, 0.85])
        ax.set_xlim([0, 1.5])
    ax.set_xlabel('resolution', fontsize=fontsize_labels)
    ax.set_ylabel('mean accuracy', fontsize=fontsize_labels)
    # Increase ticks label size
    ax.tick_params(axis='both', which='major', labelsize=fontsize_ticks)
    ax.legend()
    sns.despine(fig=fig, ax=ax)
    plt.tight_layout()

    fig.savefig(os.path.join(
        save_folder,
        'Mean{}_OptimalGenePercKeep_of_{}_{}.{}'.format(metric, key, default_title, fileformat)))
    plt.close(fig=fig)


def plot_resolution_vs_overallmetric_optimalgenekeep(
        dict_res_hvg, dict_res_hig, dict_res_std, save_folder, default_title, metric='acc', key=''):
    """
    Mean acc = acc overall seeds
    Overall Mean acc = Mean acc over all resolutions

    Parameters
    ----------
    dict_res_hvg
    dict_res_hig
    dict_res_std
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
    ind_max_mean_acc_std = np.nanargmax(dict_res_std[obs_name])

    # Apply stats test to test whether the distribution differ significantly
    stat_p_value_hrg_scanpy = wilcoxon(dict_res_hig[obs_name], dict_res_hvg[obs_name], alternative='greater')
    stat_p_value_hrg_sd = wilcoxon(dict_res_hig[obs_name], dict_res_std[obs_name], alternative='greater')

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
    ax.errorbar(x_values,
                dict_res_std[obs_name],
                yerr=dict_res_std[std_obs_name],
                label='std', fmt='-o', color='darkred', alpha=0.7)
    # Mark highest overall acc for each method with vertical line
    ax.axvline(x=x_values[ind_max_mean_acc_hvg], ymin=0, ymax=1, color='darkblue', ls='--', lw=0.5)
    ax.axvline(x=x_values[ind_max_mean_acc_std], ymin=0, ymax=1, color='darkred', ls='--', lw=0.5)
    ax.axvline(x=x_values[ind_max_mean_acc_hig], ymin=0, ymax=1, color='darkorange', ls='--', lw=0.5)
    plt.xticks(list(plt.xticks()[0]) + [x_values[ind_max_mean_acc_hvg], x_values[ind_max_mean_acc_std],
                                        x_values[ind_max_mean_acc_hig]])
    # Add Statistic result
    ax.text(x=0.35, y=0.45, s="p{} = {:.2e}\np{} = {:.2e}".format(
        r'$_{HIG\; vs\; HVG}$', stat_p_value_hrg_scanpy[1], r'$_{HIG\; vs\; std}$', stat_p_value_hrg_sd[1]),
            fontsize=10, horizontalalignment='center', verticalalignment='center')
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
        'Mean{}_OptimalGenePercKeep_of_{}_{}.{}'.format(metric, key, default_title, fileformat)))
    plt.close(fig=fig)


def plot_genepercentagetokeep_vs_overallmetricoverallresolutions(
        dict_res_hvg, dict_res_hig, dict_res_std, save_folder, default_title, metric='acc',
        alternative='greater', y_text=0.80):
    """
    Mean acc = acc overall seeds
    Overall Mean acc = Mean acc over all resolutions

    Parameters
    ----------
    dict_res_hvg
    dict_res_hig
    dict_res_std
    save_folder
    default_title
    metric
    alternative
    y_text

    Returns
    -------

    """
    # Apply stats test to test whether the distribution differ significantly;
    # Wilcoxon signed-rank test is used for paired or related samples where each observation in one group
    # corresponds to a specific observation in the other group.
    result_hig_vs_hvg = wilcoxon(dict_res_hig['mean_overall_{}_overall_resolutions'.format(metric)],
                                 dict_res_hvg['mean_overall_{}_overall_resolutions'.format(metric)],
                                 alternative=alternative)
    result_hig_vs_std = wilcoxon(dict_res_hig['mean_overall_{}_overall_resolutions'.format(metric)],
                                 dict_res_std['mean_overall_{}_overall_resolutions'.format(metric)],
                                 alternative=alternative)

    x_values = np.unique(dict_res_hvg['gene percentage to keep'])
    obs_name = 'mean_overall_{}_overall_resolutions'.format(metric)
    std_obs_name = 'std_overall_{}_overall_resolutions'.format(metric)

    # Get optimal overall_acc_overall_resolutions value
    if alternative == 'greater':
        ind_mean_acc_hvg = np.nanargmax(dict_res_hvg[obs_name])
        ind_mean_acc_hig = np.nanargmax(dict_res_hig[obs_name])
        ind_mean_acc_std = np.nanargmax(dict_res_std[obs_name])
    else:
        ind_mean_acc_hvg = np.nanargmin(dict_res_hvg[obs_name])
        ind_mean_acc_hig = np.nanargmin(dict_res_hig[obs_name])
        ind_mean_acc_std = np.nanargmin(dict_res_std[obs_name])

    print('Highest {} for HRG over all resolutions: {:.2f}% {} {:.2f}'.format(
        metric, dict_res_hig[obs_name][ind_mean_acc_hig],
        r'$\pm$', dict_res_hig[std_obs_name][ind_mean_acc_hig]))
    print('Highest {} for SCANPY over all resolutions: {:.2f}% {} {:.2f}'.format(
        metric, dict_res_hvg[obs_name][ind_mean_acc_hvg],
        r'$\pm$', dict_res_hvg[std_obs_name][ind_mean_acc_hvg]))
    print('Highest {} for SD over all resolutions: {:.2f}% {} {:.2f}'.format(
        metric, dict_res_std[obs_name][ind_mean_acc_std],
        r'$\pm$', dict_res_std[std_obs_name][ind_mean_acc_std]))

    fig, ax = plt.subplots()
    ax.grid(linewidth=0.5)
    ax.errorbar(x_values, dict_res_hvg[obs_name],
                yerr=dict_res_hvg[std_obs_name],
                label='HVG', fmt='-o', color='darkblue', alpha=0.7)
    ax.errorbar(x_values, dict_res_hig[obs_name],
                yerr=dict_res_hig[std_obs_name],
                label='HIG', fmt='-o', color='darkorange', alpha=0.7)
    # std
    ax.errorbar(x_values, dict_res_std[obs_name],
                yerr=dict_res_std[std_obs_name],
                label='std', fmt='-o', color='darkred', alpha=0.7)
    # Mark highest overall acc for each method with vertical line
    ax.axvline(x=x_values[ind_mean_acc_hvg], ymin=0, ymax=1, color='darkblue', ls='--', lw=0.5)
    ax.axvline(x=x_values[ind_mean_acc_std], ymin=0, ymax=1, color='darkred', ls='--', lw=0.5)
    ax.axvline(x=x_values[ind_mean_acc_hig], ymin=0, ymax=1, color='darkorange', ls='--', lw=0.5)
    # Add Statistic result
    ax.text(x=20, y=3, s="p{} = {:.2e}\np{} = {:.2e}".format(
        r'$_{HIG\; vs\; HVG}$', result_hig_vs_hvg[1], r'$_{HIG\; vs\; std}$', result_hig_vs_std[1]),
            fontsize=10, horizontalalignment='center', verticalalignment='center')
    plt.xticks(list(plt.xticks()[0]) + [x_values[ind_mean_acc_hvg], x_values[ind_mean_acc_std],
                                        x_values[ind_mean_acc_hig]])
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


def plot_resolution_vs_overallmetric_overallseeds_optimalgenekeep(
        dict_res_hvg, dict_res_hig, dict_res_std, save_folder, default_title, metric='acc', key='',
        alternative='greater', y_text=0.45):
    """
    Mean acc = acc overall seeds
    Overall Mean acc = Mean acc over all resolutions

    Parameters
    ----------
    dict_res_hvg
    dict_res_hig
    dict_res_std
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
        ind_mean_acc_std = np.nanargmax(dict_res_std[obs_name])
    else:
        ind_mean_acc_hvg = np.nanargmin(dict_res_hvg[obs_name])
        ind_mean_acc_hig = np.nanargmin(dict_res_hig[obs_name])
        ind_mean_acc_std = np.nanargmin(dict_res_std[obs_name])

    # Apply stats test to test whether the distribution differ significantly;
    # Wilcoxon signed-rank test is used for paired or related samples where each observation in one group
    # corresponds to a specific observation in the other group.
    result_hig_vs_hvg = wilcoxon(dict_res_hig[obs_name],
                                 dict_res_hvg[obs_name],
                                 alternative=alternative)
    result_hig_vs_std = wilcoxon(dict_res_hig[obs_name],
                                 dict_res_std[obs_name],
                                 alternative=alternative)

    fig, ax = plt.subplots()
    ax.grid(linewidth=0.5)
    ax.errorbar(x_values, dict_res_hvg[obs_name], yerr=dict_res_hvg[std_obs_name],
                label='HVG', fmt='-o', color='darkblue', alpha=0.7)
    ax.errorbar(x_values, dict_res_hig[obs_name], yerr=dict_res_hig[std_obs_name],
                label='HIG', fmt='-o', color='darkorange', alpha=0.7)
    ax.errorbar(x_values, dict_res_std[obs_name], yerr=dict_res_hig[std_obs_name],
                label='std', fmt='-o', color='darkred', alpha=0.7)
    # Mark highest overall acc for each method with vertical line
    ax.axvline(x=x_values[ind_mean_acc_hvg], ymin=0, ymax=5, color='darkblue', ls='--', lw=0.5)
    ax.axvline(x=x_values[ind_mean_acc_std], ymin=0, ymax=5, color='darkred', ls='--', lw=0.5)
    ax.axvline(x=x_values[ind_mean_acc_hig], ymin=0, ymax=5, color='darkorange', ls='--', lw=0.5)
    plt.xticks(list(plt.xticks()[0]) + [x_values[ind_mean_acc_hvg], x_values[ind_mean_acc_std],
                                        x_values[ind_mean_acc_hig]])
    # Add Statistic result
    ax.text(x=0.8, y=y_text, s="p{} = {:.2e}\np{} = {:.2e}".format(
        r'$_{HIG\; vs\; HVG}$', result_hig_vs_hvg[1], r'$_{HIG\; vs\; std}$', result_hig_vs_std[1]),
            fontsize=10, horizontalalignment='center', verticalalignment='center')
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
        'OverallMean{}_OptimalGenePercKeep_of_{}_{}.{}'.format(metric, key, default_title, fileformat)))
    plt.close(fig=fig)


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


class FSObject:
    def __init__(self, file, gene_names, method_name, save_folder, num_ticks):
        with open(file, 'rb') as ff:
            self.dict_method = pickle.load(ff)

        self.opt_geneperckeep = dict()
        self.gene_names = gene_names
        self.method = method_name

        self.method_color = {'HIG': 'darkorange', 'HVG':'darkblue', 'std': 'darkred'}

        self.highest_mean_resolution_acc = None
        self.normed_freq = {}
        self.normed_sorted_freq = {}

        self.output_dir = os.path.join(save_folder, method_name)
        os.makedirs(self.output_dir, exist_ok=True)

        self.num_ticks = num_ticks

    def get_optimal_gene_percentage(self, metric='acc'):
        x_values = np.unique(self.dict_method['gene percentage to keep'])
        # Get max overall_acc_overall_resolutions value
        ind_max_mean_metric = np.nanargmax(self.dict_method['mean_overall_{}_overall_resolutions'.format(metric)])

        self.opt_geneperckeep[metric] = x_values[ind_max_mean_metric]

    def get_overall_metric_overall_resolutions(self, metric='acc'):
        """

        Parameters
        ----------
        metric

        Returns
        -------

        """
        mean_overall_acc_overall_resolutions = []
        st_overall_acc_overall_resolutions = []

        df_res = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in self.dict_method.items()]))
        df_res = df_res[['seed', 'resolution', 'gene percentage to keep', metric]]

        for gene_to_kepp_per in df_res['gene percentage to keep'].unique():
            # Calculate acc of overall_resolutions per gene percentage to keep
            mean_overall_acc_overall_resolutions.append(
                df_res.loc[(df_res['gene percentage to keep'] == gene_to_kepp_per), metric].mean())
            st_overall_acc_overall_resolutions.append(
                df_res.loc[(df_res['gene percentage to keep'] == gene_to_kepp_per), metric].std())

        self.dict_method['mean_overall_{}_overall_resolutions'.format(metric)] = mean_overall_acc_overall_resolutions
        self.dict_method['std_overall_{}_overall_resolutions'.format(metric)] = st_overall_acc_overall_resolutions

    def get_mean_metric_overall_seeds_for_optimal_genepercentagetokeep(self, metric='acc', key=''):
        if metric not in self.opt_geneperckeep.keys():
            self.get_optimal_gene_percentage(metric=metric)

        mean_acc_overall_seeds = []
        std_acc_overall_seeds = []

        df_res = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in self.dict_method.items()]))

        df_res = df_res[['seed', 'resolution', 'gene percentage to keep', metric]]
        df_res = df_res.loc[df_res['gene percentage to keep'] == self.opt_geneperckeep[metric], :]

        for res in df_res['resolution'].unique():
            mean_acc_overall_seeds.append(df_res.loc[df_res['resolution'] == res, metric].mean())
            std_acc_overall_seeds.append(df_res.loc[df_res['resolution'] == res, metric].std())

        self.dict_method[
            'mean_{}_overall_seeds_optimal_geneperkeep_{}'.format(metric, key)] = mean_acc_overall_seeds
        self.dict_method[
            'std_{}_overall_seeds_optimal_geneperkeep_{}'.format(metric, key)] = std_acc_overall_seeds

    def get_probability_gene_topnpercent(self, metric):
        # Calculate probability that a genes falls into top 10% of feature selection method (HVG & HIG)
        df_res = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in self.dict_method.items()]))
        df_res = df_res[['gene percentage to keep', 'seed', 'resolution', 'index genes']]
        # overall seeds and resolutions for optimal gene percentage to keep value
        df_res = df_res.loc[df_res['gene percentage to keep'] == 10, :]

        gene_names = np.asarray(self.gene_names)
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
        ind_sorted = np.argsort(sorted_freq)[::-1]  # sorted_freq.sort(reverse=True)
        sorted_freq = np.asarray(sorted_freq)[ind_sorted]

        self.normed_freq[metric] = list(np.asarray(unzipped_freq[0])[ind_sorted])  # unzipped_freq[0]
        self.normed_sorted_freq[metric] = list(sorted_freq)

    def get_mean_std_acc_ch_index(self):
        mean_ch_index = []
        st_ch_index = []

        num_unique_resolutuions = np.unique(self.dict_method['resolution'])
        for resolution_val in num_unique_resolutuions:
            ind_clusters = np.where(np.asarray(self.dict_method['resolution']) == resolution_val)[0]
            # Get Mean and std  of Davies Boulding score over all gene percentage to keep
            nclusters_db_score = list(map(self.dict_method['Calinski-Harabasz Index'].__getitem__, ind_clusters))
            mean_ch_index.append(np.nanmean(nclusters_db_score))
            st_ch_index.append(np.nanstd(nclusters_db_score))

        self.dict_method['mean_calinski_haravasz_index'] = mean_ch_index
        self.dict_method['std_calinski_haravasz_index'] = st_ch_index

    # ============ Plots
    def plot_probability_gene_to_be_selected(self, metric):
        genes = list(self.normed_freq[metric])
        color = self.method_color[self.method]
        xticks = [0, 50, 100, 150, 200] + list(np.arange(500, len(genes), self.num_ticks))
        ind_xticks = list(xticks)

        fig, ax = plt.subplots()
        ax.grid(linewidth=0.5)
        ax.plot(np.arange(0, len(genes)), self.normed_sorted_freq[metric], color=color, ls='-')
        ax.set_xticks(xticks)
        ax.set_xticklabels(np.asarray(genes)[ind_xticks])
        ax.set_ylim(0, 1.02)
        ax.set_ylabel('Probability', fontsize=fontsize_labels)
        ax.set_xlabel('Genes selected by {}'.format(self.method), fontsize=fontsize_labels)
        # Increase ticks label size
        ax.tick_params(axis='both', which='major', labelsize=10)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45)
        sns.despine(ax=ax)
        plt.tight_layout()

        plt.savefig(os.path.join(self.output_dir, 'Genes selected by {}.{}'.format(self.method, fileformat)))
        plt.close()

    def plot_genepercentagetokeep_vs_overallaccoverallresolutions_method(self, metric):
        """
        Mean acc = acc overall seeds
        Overall Mean acc = Mean acc over all resolutions

        Parameters
        ----------

        Returns
        -------

        """
        color = self.method_color[self.method]
        x_values = np.unique(self.dict_method['gene percentage to keep'])
        # Get max overall_metric_overall_resolutions value
        if metric != 'Davies Bouldin Score':
            ind_mean_acc_method = np.nanargmax(self.dict_method['mean_overall_{}_overall_resolutions'.format(metric)])
            best_value = np.max(self.dict_method['mean_overall_{}_overall_resolutions'.format(metric)])
            y_pos = np.min(self.dict_method['mean_overall_{}_overall_resolutions'.format(metric)])
        else:
            ind_mean_acc_method = np.nanargmin(self.dict_method['mean_overall_{}_overall_resolutions'.format(metric)])
            best_value = np.min(self.dict_method['mean_overall_{}_overall_resolutions'.format(metric)])
            y_pos = np.max(self.dict_method['mean_overall_{}_overall_resolutions'.format(metric)])

        fig, ax = plt.subplots()
        ax.grid(linewidth=0.5)
        ax.errorbar(x_values,
                    self.dict_method['mean_overall_{}_overall_resolutions'.format(metric)],
                    yerr=self.dict_method['std_overall_{}_overall_resolutions'.format(metric)],
                    label=self.method, fmt='-o', color=color, alpha=0.7)
        # Mark highest overall acc for each method with vertical line
        ax.axvline(x=x_values[ind_mean_acc_method], ymin=0, ymax=1, color='darkblue', ls='--', lw=0.5)
        ax.text(x=x_values[ind_mean_acc_method], y=y_pos, s="{:.2f}".format(best_value),
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.5),
                horizontalalignment='center',
                verticalalignment='center')
        plt.xticks(list(plt.xticks()[0]) + [x_values[ind_mean_acc_method], x_values[ind_mean_acc_method]])
        if metric == 'acc':
            ax.set_ylim([45, 90])
        ax.set_xlim([0, 100.5])
        ax.set_xlabel('gene percentage to keep', fontsize=fontsize_labels)
        ax.set_ylabel('overall mean {}'.format(metric), fontsize=fontsize_labels)
        # Increase ticks label size
        ax.tick_params(axis='both', which='major', labelsize=fontsize_ticks)
        ax.legend()
        sns.despine(fig=fig, ax=ax)
        plt.tight_layout()

        fig.savefig(os.path.join(
            self.output_dir,
            'Mean{}_OverallResolutions_Overallseeds_{}.{}'.format(metric, self.method, fileformat)))
        plt.close(fig=fig)

    def plot_davies_bouldin_score(self):
        num_unique_resolutuions = np.unique(self.dict_method['resolution'])
        fig, ax = plt.subplots()
        ax.errorbar(num_unique_resolutuions, self.dict_method['mean_Davies Bouldin Score'],
                    yerr=self.dict_method['std_Davies Bouldin Score'], fmt='-o')
        ax.set_xlabel('Leiden resolution', fontsize=fontsize_labels)
        ax.set_ylabel('Davies Bouldin index', fontsize=fontsize_labels)
        # Increase ticks label size
        ax.tick_params(axis='both', which='major', labelsize=fontsize_ticks)
        sns.despine()
        plt.tight_layout()

        fig.savefig(os.path.join(
            self.output_dir,
            'DaviesBouldinScore_{}.{}'.format(self.method, fileformat)))
        plt.close(fig=fig)

    def plot_calinski_harabasz_index(self):
        """
        Higher value of CH index means the clusters are dense and well separated, although there is no
        “acceptable” cut-off value. We need to choose that solution which gives a peak or at least an
        abrupt elbow on the line plot of CH indices. On the other hand, if the line is smooth
        (horizontal or ascending or descending) then there is no such reason to prefer one solution over others.

        Parameters
        ----------

        Returns
        -------

        """
        num_unique_resolutuions = np.unique(self.dict_method['resolution'])
        fig, ax = plt.subplots()
        ax.errorbar(num_unique_resolutuions, self.dict_method['mean_calinski_haravasz_index'],
                    yerr=self.dict_method['std_calinski_haravasz_index'], fmt='-o')
        ax.set_xlabel('Leiden resolution', fontsize=fontsize_labels)
        ax.set_ylabel('Calinski-Harabasz Index', fontsize=fontsize_labels)
        # Increase ticks label size
        ax.tick_params(axis='both', which='major', labelsize=fontsize_ticks)
        sns.despine()
        plt.tight_layout()

        fig.savefig(os.path.join(
            self.output_dir,
            'CalinskiHarabaszIndex_{}.{}'.format(self.method, fileformat)))
        plt.close(fig=fig)

    def plot_errorbars(self):
        unique_percentage = np.unique(self.dict_method['gene percentage to keep'])
        num_unique_percentage = len(unique_percentage)
        num_unique_resolutions = np.unique(self.dict_method['resolution'])
        # per cluster
        max_mean_acc_cluster = dict.fromkeys(list(num_unique_resolutions))
        counter = 0
        for ind in range(0, len(num_unique_resolutions)):
            fig, ax = plt.subplots()
            ax.grid(linewidth=0.5)
            ax.errorbar(np.unique(self.dict_method['gene percentage to keep']),
                        self.dict_method['mean_overall_acc'][counter:counter + num_unique_percentage],
                        yerr=self.dict_method['std_overall_acc'][counter:counter + num_unique_percentage], fmt='-o')
            ax.set_ylabel('mean accuracy', fontsize=fontsize_labels)
            ax.set_xlabel('gene percentage to keep', fontsize=fontsize_labels)
            ax.set_title('Leiden resolution {}'.format(num_unique_resolutions[ind]), fontsize=fontsize_title)
            # Increase ticks label size
            ax.tick_params(axis='both', which='major', labelsize=fontsize_ticks)
            sns.despine()
            fig.savefig(os.path.join(
                self.output_dir,
                'MeanAccOverallSeeds_{}_Resolution_{}.{}'.format(
                    self.method, num_unique_resolutions[ind], fileformat)))
            plt.close(fig=fig)

            ind_max_mean_acc = np.nanargmax(
                self.dict_method['mean_overall_acc'][counter:counter + num_unique_percentage])
            max_mean_acc_cluster[num_unique_resolutions[ind]] = (
                self.dict_method['mean_overall_acc'][counter:counter + num_unique_percentage][ind_max_mean_acc],
                unique_percentage[ind_max_mean_acc])

            counter += num_unique_percentage

        self.highest_mean_resolution_acc = max_mean_acc_cluster


def save_top_robust_selected_genes(obj_hvg, obj_hig, obj_std, save_folder, metric='acc'):
    # save dataframe = TODO
    dict_selected_genes = {'HVG': obj_hvg.normed_freq[metric], 'HIG': obj_hig.normed_freq[metric],
                           'std': obj_std.normed_freq[metric]}
    unique_hvg = list(set(dict_selected_genes['HVG']) - set(dict_selected_genes['HIG']) - set(dict_selected_genes['std']))
    unique_hig = list(set(dict_selected_genes['HIG']) - set(dict_selected_genes['HVG']) - set(dict_selected_genes['std']))
    unique_std = list(set(dict_selected_genes['std']) - set(dict_selected_genes['HVG']) - set(dict_selected_genes['HIG']))
    dict_selected_genes['unique_HIG'] = unique_hig
    dict_selected_genes['unique_HVG'] = unique_hvg
    dict_selected_genes['unique_std'] = unique_std
    shared = list(set(dict_selected_genes['HVG']) & set(dict_selected_genes['HIG']) & set(dict_selected_genes['std']))
    # np.intersect1d(dict_selected_genes['HVG'], dict_selected_genes['HIG'], dict_selected_genes['std'])
    dict_selected_genes['Shared_genes'] = shared

    df_selected_genes = pd.DataFrame.from_dict(dict([(k, pd.Series(v)) for k, v in dict_selected_genes.items()]))

    df_selected_genes.to_excel(os.path.join(save_folder, 'Top10per_Selected_Genes.xlsx'))

    # # Save frequency, normed frequency, and color method to dataframe
    # dict_normed_freq = {'HVG': obj_hvg.normed_sorted_freq[metric], 'HIG': obj_hig.normed_sorted_freq[metric],
    #                     'std': obj_std.normed_sorted_freq[metric]}
    #
    # # Create plot all genes vs probability, set probability to zero for nonoccuring genes
    # all_genes_selected = np.unique(obj_hvg.normed_freq[metric] + obj_hig.normed_freq[metric] + obj_std.normed_freq[metric])
    # df_all = pd.DataFrame(index=all_genes_selected)
    # df_scanpy = pd.DataFrame(obj_hvg.normed_sorted_freq[metric], index=obj_hvg.normed_freq[metric], columns=['SCANPY'])
    # df_hrg = pd.DataFrame(obj_hig.normed_sorted_freq[metric], index=obj_hig.normed_freq[metric], columns=['HRG'])
    # df_sd = pd.DataFrame(obj_std.normed_sorted_freq[metric], index=obj_std.normed_freq[metric], columns=['SD'])
    #
    # # Combine dataframes
    # df_all = pd.concat([df_all, df_scanpy], axis=1)
    # df_all = pd.concat([df_all, df_hrg], axis=1)
    # df_all = pd.concat([df_all, df_sd], axis=1)
    #
    # xticks = list(np.arange(0, len(all_genes_selected), 500))
    # ind_xticks = list(xticks)
    #
    # fig, ax = plt.subplots()
    # ax.plot(np.arange(0, len(all_genes_selected)), df_all['SCANPY'], color=obj_hvg.method_color['HVG'], ls='-')
    # ax.plot(np.arange(0, len(all_genes_selected)), df_all['HRG'], color=obj_hig.method_color['HIG'], ls='-')
    # ax.plot(np.arange(0, len(all_genes_selected)), df_all['SD'], color=obj_std.method_color['std'], ls='-')
    # ax.set_xticks(xticks)
    # ax.set_xticklabels(np.asarray(all_genes_selected)[ind_xticks])
    # ax.set_ylim(0, 1.02)
    # ax.set_ylabel('Probability', fontsize=fontsize_labels)
    # ax.set_xlabel('Selected genes', fontsize=fontsize_labels)
    # # Increase ticks label size
    # ax.tick_params(axis='both', which='major', labelsize=10)
    # ax.set_xticklabels(ax.get_xticklabels(), rotation=45)
    # sns.despine(ax=ax)
    # plt.tight_layout()


if __name__ == '__main__':
    path_input = "/Volumes/CH__data/Projects/Eyerich_AG_projects/BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer/analysis/Molecular_subtypes/output/Leiden"
    output_dir = os.path.join(path_input, 'comparing_edgeR_HIG_HVG_HVGold', str(date.today()))
    os.makedirs(output_dir, exist_ok=True)

    # Load gene names
    adata = sc.read(os.path.join(
        data_root, 'Molecular_subtypes', 'input', 'h5_files', 'LESION',
        'Lesion_RNAseq_20210720_patient_meta_data_v04__CH_20220209_TherapyResponse_v04_v220715.h5'))
    var_names = list(adata.var_names)
    del adata

    hig_obj = FSObject(file=os.path.join(
        path_input, '2023-02-12', 'resolution_[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1. ]__featureselection_HIG',
        'Scores_GridSearch_9_minPercentage0.5_maxPercentage100.0_minResolution0.1_maxResolution1.0_.pkl'),
        gene_names=var_names, method_name='HIG', save_folder=output_dir, num_ticks=500)

    date_scanpy = '2025-07-09'  # '2023-02-13'
    hvg_obj = FSObject(file=os.path.join(
        path_input, date_scanpy, 'resolution_[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1. ]__featureselection_HVG',
        'Scores_GridSearch_9_minPercentage0.5_maxPercentage100.0_minResolution0.1_maxResolution1.0_.pkl'),
        gene_names=var_names, method_name='HVG', save_folder=output_dir, num_ticks=1000)

    # Load old HVG result where Scanpy operated on normalised gene counts (gene profiles)
    std_obj = FSObject(file=os.path.join(
        path_input, '2023-02-13', 'resolution_[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1. ]__featureselection_HVG',
        'Scores_GridSearch_9_minPercentage0.5_maxPercentage100.0_minResolution0.1_maxResolution1.0_.pkl'),
        gene_names=var_names, method_name='HVGold', save_folder=output_dir, num_ticks=1000)

    # ===================== start analysis
    dict_obj = {'HIG': hig_obj, 'HVG': hvg_obj, 'HVGold': std_obj}
    counts_methods = {}
    for tmp_obj_key in dict_obj.keys():
        tmp_obj = dict_obj[tmp_obj_key]
        tmp_obj.dict_method['acc'] = list(map(lambda x: x * 100, tmp_obj.dict_method['acc']))
        for metric_tmp in ['acc', 'Davies Bouldin Score', 'Calinski-Harabasz Index']:
            # Get mean metric overall seeds and resolutions
            tmp_obj.get_overall_metric_overall_resolutions(metric=metric_tmp)
            tmp_obj.get_mean_std_acc_ch_index()
            tmp_obj.get_optimal_gene_percentage(metric=metric_tmp)

            # ================ Comparing Feature selection methods HVG vs HIG using likelihood a gene falls into 10%
            # Calculate probability that a genes falls into top 10% of feature selection method (HVG & HIG)
            tmp_obj.get_probability_gene_topnpercent(metric=metric_tmp)
            tmp_obj.plot_probability_gene_to_be_selected(metric=metric_tmp)

            for key in ['acc', 'Davies Bouldin Score', 'Calinski-Harabasz Index']:
                # Using optimal Gene percentage to keep derived from metric (metric_tmp)
                # evaluated against another metric key
                # for same metric:
                # Find the highest resolution for gene percentage to keep parameter:
                # x = resolution y = mean acc over all seeds
                tmp_obj.get_mean_metric_overall_seeds_for_optimal_genepercentagetokeep(
                    metric=metric_tmp, key='best_{}'.format(key))

            # Identify the highest mean acc (mean by seeds) across gene percentage to keep
            tmp_obj.plot_errorbars()
            tmp_obj.plot_genepercentagetokeep_vs_overallaccoverallresolutions_method(metric=metric_tmp)
            tmp_obj.plot_calinski_harabasz_index()
            tmp_obj.plot_davies_bouldin_score()

            # Read out gene percentage to keep
            gene_per_keep = []
            for dict_item in tmp_obj.highest_mean_resolution_acc.items():
                gene_per_keep.append(dict_item[1][1])
            counts_method = collections.Counter(gene_per_keep)
            counts_methods[tmp_obj_key] = counts_method

    # TODO
    save_top_robust_selected_genes(obj_hvg=hvg_obj, obj_hig=hig_obj, obj_std=std_obj, save_folder=output_dir)

    # Take gene percentage to keep which occurred the most in counter
    # Draw plot gene percentage to keep vs mean overall acc comparing HVG and HIG per resolution
    plot_errorbars_HIG_vs_HVG_std(
        dict_res_hvg=hvg_obj.dict_method, dict_res_hig=hig_obj.dict_method, dict_res_std=std_obj.dict_method,
        save_folder=output_dir, default_title='HVG_vs_HIG_HVGold')

    plot_genepercentagetokeep_vs_overallaccoverallresolutions(
        dict_res_hvg=hvg_obj.dict_method, dict_res_hig=hig_obj.dict_method, dict_res_std=std_obj.dict_method,
        save_folder=output_dir, default_title='HVG_vs_HIG_HVGold')

    """ 
    Note: 
    - Optimal gene percentage to keep parameter is 6 for HVG and 7 for HIG
    - Mann WhitneyU test rejects null hypothesis -> Methods HVG and HIG differ significantly
    - HIG is significantly better than HVG -> use HIG for feature selection
    """
    plot_resolution_vs_overallacc_optimalgenekeep(
        dict_res_hvg=hvg_obj.dict_method, dict_res_hig=hig_obj.dict_method, dict_res_std=std_obj.dict_method,
        save_folder=output_dir, default_title='HVG_vs_HIG_HVGold', metric='acc', key='best_acc')

    # ================ Calinski Harabasz Index
    plot_genepercentagetokeep_vs_overallmetricoverallresolutions(
        dict_res_hvg=hvg_obj.dict_method, dict_res_hig=hig_obj.dict_method, dict_res_std=std_obj.dict_method,
        save_folder=output_dir, default_title='HVG_vs_HIG_HVGold',
        metric='Calinski-Harabasz Index', alternative='greater', y_text=2)

    plot_resolution_vs_overallmetric_optimalgenekeep(
        dict_res_hvg=hvg_obj.dict_method, dict_res_hig=hig_obj.dict_method, dict_res_std=std_obj.dict_method,
        save_folder=output_dir, default_title='HVG_vs_HIG_HVGold',
        metric='Calinski-Harabasz Index', key='best_Calinski-Harabasz Index')

    plot_resolution_vs_overallmetric_overallseeds_optimalgenekeep(
        dict_res_hvg=hvg_obj.dict_method, dict_res_hig=hig_obj.dict_method, dict_res_std=std_obj.dict_method,
        save_folder=output_dir, default_title='HVG_vs_HIG_HVGold',
        metric='Calinski-Harabasz Index', alternative='greater', y_text=2, key='best_Calinski-Harabasz Index')

    plot_resolution_vs_overallmetric_optimalgenekeep(
        dict_res_hvg=hvg_obj.dict_method, dict_res_hig=hig_obj.dict_method, dict_res_std=std_obj.dict_method,
        save_folder=output_dir, default_title='HVG_vs_HIG_HVGold',
        metric='Calinski-Harabasz Index', key='best_acc')

    plot_resolution_vs_overallmetric_overallseeds_optimalgenekeep(
        dict_res_hvg=hvg_obj.dict_method, dict_res_hig=hig_obj.dict_method, dict_res_std=std_obj.dict_method,
        save_folder=output_dir, default_title='HVG_vs_HIG_HVGold',
        metric='Calinski-Harabasz Index', alternative='greater', y_text=2, key='best_acc')

    # ================ Davies Bouldin Score
    plot_genepercentagetokeep_vs_overallmetricoverallresolutions(
        dict_res_hvg=hvg_obj.dict_method, dict_res_hig=hig_obj.dict_method, dict_res_std=std_obj.dict_method,
        save_folder=output_dir, default_title='HVG_vs_HIG_HVGold',
        metric='Davies Bouldin Score', alternative='less', y_text=3.4)

    plot_resolution_vs_overallmetric_optimalgenekeep(
        dict_res_hvg=hvg_obj.dict_method, dict_res_hig=hig_obj.dict_method, dict_res_std=std_obj.dict_method,
        save_folder=output_dir, default_title='HVG_vs_HIG_HVGold',
        metric='Davies Bouldin Score', key='best_Davies Bouldin Score')

    plot_resolution_vs_overallmetric_overallseeds_optimalgenekeep(
        dict_res_hvg=hvg_obj.dict_method, dict_res_hig=hig_obj.dict_method, dict_res_std=std_obj.dict_method,
        save_folder=output_dir, default_title='HVG_vs_HIG_HVGold',
        metric='Davies Bouldin Score', alternative='less', y_text=2.5, key='best_Davies Bouldin Score')

    plot_resolution_vs_overallmetric_optimalgenekeep(
        dict_res_hvg=hvg_obj.dict_method, dict_res_hig=hig_obj.dict_method, dict_res_std=std_obj.dict_method,
        save_folder=output_dir, default_title='HVG_vs_HIG_HVGold',
        metric='Davies Bouldin Score', key='best_acc')

    plot_resolution_vs_overallmetric_overallseeds_optimalgenekeep(
        dict_res_hvg=hvg_obj.dict_method, dict_res_hig=hig_obj.dict_method, dict_res_std=std_obj.dict_method,
        save_folder=output_dir, default_title='HVG_vs_HIG_HVGold',
        metric='Davies Bouldin Score', alternative='less', y_text=2.5, key='best_acc')
