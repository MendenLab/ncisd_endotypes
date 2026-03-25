from scripts.utils import add_colors, add_endotypes

import scanpy as sc
import numpy as np
import pandas as pd
import os
from datetime import date
from operator import itemgetter

import matplotlib.pyplot as plt

plt.rcParams["font.family"] = "Arial"
fileformat = '.pdf'


class Evaluation:
    def __init__(self, save_folder):
        self._save_folder = save_folder

        self.probability_histogram = None

    def set_probability_histogram(self, true_labels, pred_labels):

        # Create crosstab showing the overlap between the two observations
        df_frequencies = pd.crosstab(true_labels, pred_labels)

        p = dict()
        # Get probabilities that a class resembles a true label (from clinician)
        for cluster_id in np.unique(pred_labels):
            # total number of samples of cluster_id
            l_label = np.count_nonzero(pred_labels == cluster_id)
            p[cluster_id] = []
            for true_class_label in list(df_frequencies.index):
                # total number of samples of disease in cluster_id
                k_label = df_frequencies.loc[true_class_label, cluster_id]
                # P(label|z) = P(z|label)P(label) / P(z) = K(label)/K
                # where K is the total number of samples in a cluster
                p[cluster_id].append((true_class_label, k_label / l_label))

        self.probability_histogram = pd.DataFrame.from_dict({k: dict(v) for k, v in p.items()}, orient='index')

    def get_match_probalities(self, pred_labels):
        assert self.probability_histogram is not None
        # Identify which predicted label matches to observation label
        df_prob_predicted_labels = self.probability_histogram.loc[pred_labels]
        # Read out most likely true label
        matched_obs_label_predicted = list(
            df_prob_predicted_labels.columns[np.argmax(np.asarray(df_prob_predicted_labels), axis=1)])

        return matched_obs_label_predicted, df_prob_predicted_labels

    def plot_probability_histogram_supplements(
            self, true_labels: np.ndarray, module_name: str, dataset: str, obs: str,
            colors: list = plt.get_cmap('jet'), str_int_truelabel_mapper: dict = None, ncol_legend: int = 1,
            figsize: tuple = (12, 4)):

        if str_int_truelabel_mapper is not None:
            mapper = {}
            labels = list(str_int_truelabel_mapper.keys())
            for ind_classlabels, val_classlabels in enumerate(np.unique(true_labels)):
                mapper[val_classlabels] = labels[ind_classlabels]
            str_labels = [mapper[label] for label in true_labels]
            mapper = {value: key for key, value in mapper.items()}
            dict_labels_colors = dict(zip(str_labels, colors))
            # reorder keys
            reordered_dict = {k: dict_labels_colors[k] for k in list(mapper.keys())}
            # read out new colors
            colors = list(itemgetter(*mapper)(reordered_dict))

        ax = self.probability_histogram.plot.bar(rot=0, stacked=True, color=colors, figsize=figsize, width=0.95)
        ax.set_ylabel('probability', fontsize=18)
        # ax.set_xlabel("{} clusters".format(module_name), fontsize=18)
        ax.set_xlabel("")
        ax.tick_params(axis='both', labelsize=15)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['left'].set_visible(True)
        ax.spines['bottom'].set_visible(True)
        # ax.set_xticklabels(list(itemgetter(*list(self.probability_histogram.index))(translate_number_label)))
        # Add legend
        legend_objs = ax.legend(title=obs.capitalize(), loc='upper center', bbox_to_anchor=(0.5, -0.05),
                                ncol=ncol_legend, prop={'size': 15}, frameon=False, title_fontsize=17)

        if (str_int_truelabel_mapper is not None) & (obs == 'Pattern'):
            # Save legend labels
            str_int_truelabel_mapper = {"UD" if k == 'undefined' else k: v for k, v in str_int_truelabel_mapper.items()}
            legends_labels = list(str_int_truelabel_mapper.keys())
            for ind_classlabels, _ in enumerate(np.unique(true_labels)):
                legend_objs.get_texts()[ind_classlabels].set_text(legends_labels[ind_classlabels])

        # Add percentages to bars
        # for p in ax.patches:
        #     width, height = p.get_width(), p.get_height()
        #     x, y = p.get_xy()
        #     if height >= 0.05:
        #         ax.text(x + width / 2,
        #                 y + height / 2,
        #                 '{:.0f} %'.format(height * 100),
        #                 horizontalalignment='center',
        #                 verticalalignment='center', fontsize=12)

        plt.savefig(os.path.join(
            self._save_folder, 'Supplements_Barplot_{}_vs_{}__{}_probabilities_percentage.pdf'.format(
                obs, module_name, dataset)), bbox_inches='tight')
        plt.close()

    def plot_probability_histogram(
            self, true_labels: np.ndarray, module_name: str, dataset: str, obs: str, translate_number_label: dict,
            colors: list = plt.get_cmap('jet'), str_int_truelabel_mapper: dict = None, ncol_legend: int = 1,
            figsize: tuple = (12, 4)):

        # Plot majority and second majority label in color and rest in grey
        # Get index of columns with highest value
        max_ind = np.argmax(self.probability_histogram, axis=1)
        # Get index of columns with second highest value
        second_highest_values = self.probability_histogram.apply(lambda row: row.nlargest(2).values[-1], axis=1)
        col_names = self.probability_histogram.columns
        # Get second highst colnames
        second_highest_colnames = []
        for ind, subtype_ind in enumerate(self.probability_histogram.index):
            second_highest_colnames.append(self.probability_histogram.loc[
                subtype_ind, self.probability_histogram.loc[
                             subtype_ind, :] == second_highest_values[subtype_ind]].index[0])
        # Old col names of highest and second highest values
        old_col_names = list(np.unique(list(np.unique(self.probability_histogram.iloc[:, max_ind].columns)) + list(
            np.unique(second_highest_colnames))))
        new_col_names = old_col_names + ['Others']

        df_tmp = pd.DataFrame(index=self.probability_histogram.index, columns=new_col_names)
        for ind, subtype_ind in enumerate(self.probability_histogram.index):
            second_highest_colname = self.probability_histogram.loc[
                subtype_ind, self.probability_histogram.loc[
                             subtype_ind, :] == second_highest_values[subtype_ind]].index[0]
            # Highest value
            df_tmp.loc[subtype_ind, col_names[max_ind[ind]]] = self.probability_histogram.loc[
                subtype_ind, col_names[max_ind[ind]]]
            # Second highest
            df_tmp.loc[subtype_ind, second_highest_colname] = self.probability_histogram.loc[
                subtype_ind, second_highest_colname]

            df_tmp.loc[subtype_ind, 'Others'] = 1 - df_tmp.loc[subtype_ind, col_names[max_ind[ind]]] - df_tmp.loc[
                subtype_ind, second_highest_colname]

        # Get colors
        if str_int_truelabel_mapper is not None:
            colors = list(itemgetter(*old_col_names)(str_int_truelabel_mapper))
            colors.append('snow')

        # Order columns for barchart
        df_tmp.columns = pd.CategoricalIndex(df_tmp.columns.values, ordered=True, categories=new_col_names)
        # Sort the columns (axis=1) by the new categorical ordering
        df_tmp = df_tmp.sort_index(axis=1)

        df1 = df_tmp.loc[['E5', 'E6', 'E7', 'E8', 'E9', 'E10'], :]
        df2 = df_tmp.loc[['E11', 'E12', 'E13'], :]
        df3 = df_tmp.loc[['E1', 'E2', 'E3', 'E4'], :]

        def plot_barchart(df_, axis_name, show_yaxis=True, show_legend=True):
            df_.plot.bar(
                rot=0, stacked=True, color=colors, figsize=figsize, width=0.95, edgecolor="black",
                linewidth=0.1, ax=axis_name, legend=show_legend)
            axis_name.set_xlabel("")
            axis_name.tick_params(axis='both', labelsize=18)
            axis_name.spines['right'].set_visible(False)
            axis_name.spines['top'].set_visible(False)
            axis_name.spines['bottom'].set_visible(True)
            if show_yaxis:
                axis_name.set_ylabel('probability', fontsize=18)
                axis_name.spines['left'].set_visible(True)
            else:
                axis_name.set_ylabel('', fontsize=18)
                axis_name.set_yticks([])
                axis_name.spines['left'].set_visible(False)

            # add percentages
            # for p in axis_name.patches:
            #     width, height = p.get_width(), p.get_height()
            #     x, y = p.get_xy()
            #     if height >= 0.05:
            #         axis_name.text(x + width / 2, y + height / 2, '{:.0f} %'.format(round(height * 100)),
            #                        horizontalalignment='center', verticalalignment='center', fontsize=15)

            if show_legend:
                # Add legend
                axis_name.legend(title=obs.capitalize(), loc='upper center', bbox_to_anchor=(0.5, -0.1),
                                 ncol=ncol_legend, prop={'size': 18}, frameon=False, title_fontsize=20)

            return axis_name

        # Split barplot into dendrogram clusters
        widths = [3, 1.5, 2]
        heights = [1]
        fig = plt.figure(figsize=(8, 6), constrained_layout=True)
        spec5 = fig.add_gridspec(ncols=3, nrows=1, width_ratios=widths, height_ratios=heights)
        # Pattern 2a-like clusters
        # fig, ax = plt.subplots(nrows=1, ncols=3)
        # ax = fig.add_axes([0.2, 0.1, 0.28, 0.72])
        ax_1 = fig.add_subplot(spec5[0, 0])
        ax_1 = plot_barchart(df_=df1, axis_name=ax_1, show_yaxis=True, show_legend=False)
        # ax_2 = fig.add_axes([0.50, 0.1, 0.14, 0.72])
        ax_2 = fig.add_subplot(spec5[0, 1])
        ax_2 = plot_barchart(df_=df2, axis_name=ax_2, show_yaxis=False, show_legend=False)
        # ax_3 = fig.add_axes([0.66, 0.1, 0.19, 0.72])
        ax_3 = fig.add_subplot(spec5[0, 2])
        ax_3 = plot_barchart(df_=df3, axis_name=ax_3, show_yaxis=False, show_legend=False)

        lines_labels = [ax.get_legend_handles_labels() for ax in fig.axes]
        lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]
        res = [idx for idx, val in enumerate(labels) if val not in labels[:idx]]

        fig.legend(np.asarray(lines)[res], np.asarray(labels)[res], title='',  # obs.capitalize()
                   loc='upper center', bbox_to_anchor=(0.5, -0.05),
                   ncol=ncol_legend, prop={'size': 24}, frameon=False)

        plt.savefig(os.path.join(
            self._save_folder, 'Barplot_{}_vs_{}__{}_probabilities_percentage.pdf'.format(
                obs, module_name, dataset)), bbox_inches='tight')
        plt.close()


def main(save_folder, general_directory):
    data_root = os.path.join(
        general_directory, "analysis", "Molecular_subtypes", "input", "h5_files", "LESION"
    )
    adata = sc.read(
        os.path.join(data_root, 'Lesion_RNAseq_20210720_patient_meta_data_v04__CH_20220209_TherapyResponse_v04_v220715_TherapyResponseCorrected__Endotypes_230620.h5'))
    endotype_colname = 'Endotypes'

    ncols = {'diag': 2, 'Pattern': 4}
    figsize = {'diag': (12, 6), 'Pattern': (12, 6)}
    color_mapper_diag = dict(zip([
        'lichen planus', 'lupus erythematosus', 'lichenoid drug reaction', 'eczema', 'prurigo simplex subacuta',
        'bullous pemphigoid', 'psoriasis', 'pityriasis rubra pilaris', 'morphea', 'venous ulcer',
        'systemic sclerosis', 'granuloma annulare', 'sarcoidosis', 'psoriasis pustulosa', 'pyoderma gangrenosum',
        'cutaneous lymphoma', 'cutaneous side effects of biologics', 'darier disease', 'keratosis lichenoides chronica',
        'erythrodermia', 'parapsoriasis', 'undefined'],
        ['orange', 'goldenrod', 'gold', 'maroon', 'firebrick', 'salmon', 'blue', 'dodgerblue', 'palevioletred',
         'mediumvioletred', 'deeppink', 'fuchsia', 'violet', 'darkviolet', 'mediumorchid', 'dimgrey',
         'darkgray', 'silver', 'grey', 'slategrey', 'lightgrey', 'whitesmoke']))
    color_mapper_pattern = dict(zip(['1', '2a', '2b', '3', '4a', '4b', '5', 'UD'],
                                    ['darkorange', 'darkred', 'tomato', 'royalblue', 'mediumvioletred',
                                     'magenta', 'darkviolet', 'grey']))
    str_int_truelabel_mapper = {'diag': color_mapper_diag, 'Pattern': color_mapper_pattern}

    for labels_name in ['diag', 'Pattern']:
        true_labels = adata.obs[labels_name].values
        if endotype_colname == 'Molecular Subtype res0.9':
            predicted_labels = np.asarray(adata.obs[endotype_colname].values) + 1
            order = list(np.arange(1, len(np.unique(predicted_labels))))
        else:
            predicted_labels = np.asarray(adata.obs[endotype_colname].map(lambda x: x.lstrip('E')).astype(int))
            order = ['E5', 'E6', 'E7', 'E8', 'E9', 'E10', 'E11', 'E12', 'E13', 'E1', 'E2', 'E3', 'E4']
        # Predicted to number labels dict
        translate_number_label = dict(zip(predicted_labels, adata.obs[endotype_colname]))
        eval_obs = Evaluation(save_folder=save_folder)
        eval_obs.set_probability_histogram(true_labels=true_labels, pred_labels=predicted_labels)
        # set index
        eval_obs.probability_histogram.index = list(adata.obs['Endotypes'].cat.categories)
        # reorder index like in Dendrogram
        eval_obs.probability_histogram = eval_obs.probability_histogram.reindex(order)

        # TODO remove text percentages
        # TODO norm by total number of samples of disease or in a pattern
        eval_obs.plot_probability_histogram(
            true_labels=true_labels, module_name='Leiden', dataset='application', obs=labels_name,
            colors=adata.uns['{}_colors'.format(labels_name)],
            str_int_truelabel_mapper=str_int_truelabel_mapper[labels_name],
            ncol_legend=ncols[labels_name], translate_number_label=translate_number_label, figsize=figsize[labels_name])
        eval_obs.plot_probability_histogram_supplements(
            true_labels=true_labels, module_name='Leiden', dataset='application', obs=labels_name,
            colors=adata.uns['{}_colors'.format(labels_name)], str_int_truelabel_mapper=None,
            ncol_legend=ncols[labels_name], figsize=figsize[labels_name])

        df_tmp = eval_obs.probability_histogram
        # df_tmp['Endotypes'] = list(itemgetter(*list(eval_obs.probability_histogram.index))(translate_number_label))
        df_tmp = (df_tmp * 100).round(2)
        df_tmp.to_excel(
            os.path.join(save_folder, 'Percentage_Endotype_per_{}.xlsx'.format(labels_name)))


if __name__ == '__main__':
    general_dir = os.path.join(
        "/Volumes",
        "CH__data",
        "Projects",
        "Eyerich_AG_projects",
        "BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer",
    )
    output_dir = os.path.join(general_dir, 'analysis', 'Molecular_subtypes', 'output', 'Figure_3BC_Histogram',
                              str(date.today()))
    os.makedirs(output_dir, exist_ok=True)
    main(save_folder=output_dir, general_directory=general_dir)
