import os
import pickle
import pandas as pd
import numpy as np


def combine_pickle_files(directory_files, output_filename):
    results = []
    for file_name in directory_files:
        with open(file_name, 'rb') as f:
            content = pickle.load(f)
            # convert to list
            results.append(content)

    # Convert output list of dictionary into dictionary of list
    dict_result = {key: [i[key] for i in results] for key in results[0]}

    # store in .pkl file
    with open(output_filename, 'wb') as out:
        pickle.dump(dict_result, out)


def get_mean_std_acc_dbscore(dict_result):
    # Get overall acc
    mean_seed_acc = []
    st_seed_acc = []
    mean_db_score = []
    st_db_score = []

    df = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in dict_result.items()]))
    for resolution_val in np.unique((dict_result['resolution'])):
        ind_clusters = np.where(np.asarray(dict_result['resolution']) == resolution_val)[0]
        # Get Mean and std  of Davies Boulding score over all gene percentage to keep
        nclusters_db_score = list(map(dict_result['Davies Bouldin Score'].__getitem__, ind_clusters))
        mean_db_score.append(np.mean(nclusters_db_score))
        st_db_score.append(np.std(nclusters_db_score))

        df_tmp = df[(df['resolution']).astype(str).str.contains(str(resolution_val))]
        # Get per seed the overall acc
        for val in np.unique((df_tmp['gene percentage to keep'])):
            ind_val = np.where(np.asarray(df_tmp['gene percentage to keep']) == val)[0]
            # Get Mean and std of Accuracy over all seeds
            seed_accs = list(df_tmp.iloc[ind_val]['acc'])
            mean_seed_acc.append(np.mean(seed_accs))
            st_seed_acc.append(np.std(seed_accs))

    dict_result['mean_overall_acc'] = mean_seed_acc
    dict_result['std_overall_acc'] = st_seed_acc
    dict_result['mean_Davies Bouldin Score'] = mean_db_score
    dict_result['std_Davies Bouldin Score'] = st_db_score

    return dict_result


if __name__ == '__main__':

    # directory_path = os.path.join(
    #     '/Volumes/CH__data/Projects/Eyerich_AG_projects',
    #     'BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer/analysis/Molecular_subtypes/output',
    #     'Leiden/2023-12-02/resolution_[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1. ]__featureselection_std')
    # output_file = os.path.join(directory_path,
    #                            'Scores_GridSearch_9_minPercentage0.5_maxPercentage100.0_minResolution0.1_maxResolution1.0_.pkl')
    analysis_date = '2025-07-11'  # '2023-12-02'  # '2025-07-09'
    selection_method = 'HVG'  # 'std'  # 'HVG'
    directory_path = os.path.join(
        '/Volumes/CH__data/Projects/Eyerich_AG_projects',
        'BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer/analysis/Molecular_subtypes/output',
        'Leiden/{}/resolution_[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1. ]__featureselection_{}'.format(
            analysis_date, selection_method))
    output_file = os.path.join(
        directory_path, 'Scores_GridSearch_9_minPercentage0.5_maxPercentage100.0_minResolution0.1_maxResolution1.0_.pkl')

    # len = 9090
    files = []
    for dirpath, dirnames, filenames in os.walk(directory_path):
        for filename in [f for f in filenames if f.endswith(".pkl")]:
            files.append(os.path.join(dirpath, filename))
    # usage
    combine_pickle_files(directory_files=files, output_filename=output_file)

    # from datetime import date
    # output_dir = os.path.join(directory_path, 'comparing_edgeR_std', str(date.today()))
    # os.makedirs(output_dir, exist_ok=True)
    with open(output_file, 'rb') as f:
        gs_std = pickle.load(f)

    gs_std = get_mean_std_acc_dbscore(gs_std)

    # store in .pkl file
    with open(output_file, 'wb') as out:
        pickle.dump(gs_std, out)
    #
    # highest_mean_resolution_acc_std = plot_errorbars(dict_res=gs_std, save_folder=output_dir, default_title='std')
    # plot_davies_bouldin_score(dict_res=gs_std, save_folder=output_dir, default_title='std')
    # gs_std = get_mean_std_acc_ch_index(dict_res=gs_std)
    # plot_calinski_harabasz_index(dict_res=gs_std, save_folder=output_dir, default_title='std')
    #
    #
    # # Read out gene percentage to keep
    # gene_per_keep_std = []
    # for dict_item in highest_mean_resolution_acc_std.items():
    #     gene_per_keep_std.append(dict_item[1][1])
    # counts_std = collections.Counter(gene_per_keep_std)
    #
    #
    # # Get mean acc overall seeds and resolutions
    # gs_std = get_overall_acc_overall_resolutions(dict_res=gs_std)
    # plot_genepercentagetokeep_vs_overallaccoverallresolutions_method(
    #     dict_res_method=gs_std, method='std', save_folder=output_dir, default_title='std')
    #
    # # Get optimal gene percentage to keep
    # x_values = np.unique(gs_std['gene percentage to keep'])
    # # Get max overall_acc_overall_resolutions value
    # ind_max_mean_acc_std = np.nanargmax(gs_std['mean_overall_acc_overall_resolutions'])
    # optimal_genekeep_std = x_values[ind_max_mean_acc_std]
    #
    # # Find highest resolution for gene percentage to keep parameter 74%:
    # # x = resolution y = mean acc over all seeds
    # gs_std = get_mean_acc_overall_seeds_for_optimal_genepercentagetokeep_value(
    #     dict_res=gs_std, opt_geneperckeep=optimal_genekeep_std)
    #
    # # ================ Calinski Harabasz Index
    # gs_std = get_overall_metric_overall_resolutions(dict_res=gs_std, metric='Calinski-Harabasz Index')
    #
    # # Using optimal Gene percentage to keep derived from Calinski-Harabasz Index
    # gs_hvg = get_mean_metric_overall_seeds_for_optimal_genepercentagetokeep(
    #     dict_res=gs_hvg, opt_geneperckeep=optimal_genekeep_hvg_chi, metric='Calinski-Harabasz Index', key='best_CHI')
