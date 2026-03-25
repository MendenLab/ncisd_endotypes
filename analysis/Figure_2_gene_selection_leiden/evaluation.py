import numpy as np
import pandas as pd
import os
import itertools
from operator import itemgetter

import matplotlib.pyplot as plt

from sklearn import metrics
# from scripts.utils.DBCV import DBCV
import hdbscan


class Evaluation:
    def __init__(self, cfg, cluster):
        self._cfg = cfg

        self._cluster = cluster

        self.probability_histogram = None

    def fit(self, pred_labels, true_labels):
        # pred_labels = self._cluster.predict(data)

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
                # Bayes' Theorem Formula: P(label|z) = P(z|label)P(label) / P(z) = K(label)/K
                # where K is the total number of samples in a cluster
                p[cluster_id].append((true_class_label, k_label / l_label))

        self.probability_histogram = pd.DataFrame.from_dict({k: dict(v) for k, v in p.items()}, orient='index')

    def predict(self, pred_labels):
        assert self.probability_histogram is not None
        # pred_labels = self._cluster.predict(data)
        # Identify which predicted label matches to observation label
        df_prob_predicted_labels = self.probability_histogram.loc[pred_labels]
        # Read out most likely true label
        matched_obs_label_predicted = list(
            df_prob_predicted_labels.columns[np.argmax(np.asarray(df_prob_predicted_labels), axis=1)])

        return matched_obs_label_predicted, df_prob_predicted_labels

    def plot_probability_histogram(
            self, true_labels: np.ndarray, module_name: str, dataset: str,
            colors: list = plt.get_cmap('jet'), str_int_truelabel_mapper: dict = None):

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

        ax = self.probability_histogram.plot.bar(rot=0, stacked=True, color=colors)
        ax.set_title('Feature selection:{}, resolution: {}'.format(
            self._cfg.feature_selection.do_feature_selection, self._cfg.cluster.resolution), fontsize=18)
        ax.set_ylabel('probability', fontsize=18)
        ax.set_xlabel("{} clusters".format(module_name), fontsize=18)
        ax.tick_params(axis='both', labelsize=15)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['left'].set_visible(True)
        ax.spines['bottom'].set_visible(True)
        legend_objs = ax.legend(title='pattern', loc='center left', bbox_to_anchor=(1, 0.5), ncol=1,
                                prop={'size': 15}, frameon=False, title_fontsize=17)
        if str_int_truelabel_mapper is not None:
            # Save legend labels
            str_int_truelabel_mapper = {"UD" if k == 'undefined' else k: v for k, v in str_int_truelabel_mapper.items()}
            legends_labels = list(str_int_truelabel_mapper.keys())
            for ind_classlabels, _ in enumerate(np.unique(true_labels)):
                legend_objs.get_texts()[ind_classlabels].set_text(legends_labels[ind_classlabels])
        plt.tight_layout()
        plt.savefig(os.path.join(self._cfg.results.save_folder, 'Barplot_{}_vs_{}__{}_probabilities.pdf'.format(
            'Pattern', module_name, dataset)))
        plt.close()

    def get_confusion_matrix(
            self, pred_labels, true_labels, feature_selection_approach, observation_name, module_name, dataset,
            str_int_truelabel_mapper, normalize=True):
        # Predict labels
        matched_obs_label_predicted, _ = self.predict(pred_labels=pred_labels)

        # Calculate confusion matrix
        cm = metrics.confusion_matrix(true_labels, matched_obs_label_predicted)
        # Report score
        # Divide each row by its sum. Then calculate the mean of the diagonal
        # "By normalizing the row we will get some probability distribution of mistakes for the pattern"
        normed_cm = cm / np.sum(cm, axis=1)[:, None]
        # replace nan with 0
        normed_cm[np.isnan(normed_cm)] = 0
        # diag_score = np.mean(np.diagonal(normed_cm))
        diag_score = np.sum(np.diagonal(cm)) / np.sum(np.sum(cm))

        print('Diagonal score: ', diag_score)

        title = 'Feature selection:{}, resolution: {}\nConfusion matrix {} vs {}'.format(
            self._cfg.feature_selection.do_feature_selection, self._cfg.cluster.resolution,
            observation_name, module_name)
        axis_title = 'score: {:.2f}%'.format(diag_score * 100)
        unique_xticks = np.arange(
            0, np.max([len(np.unique(true_labels)), len(np.unique(matched_obs_label_predicted))]))

        if normalize:
            cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]

        fig, ax = plt.subplots(figsize=(7, 6))
        m = ax.imshow(cm, interpolation='nearest', cmap=plt.cm.Oranges)
        fig.suptitle(title, fontsize=18)
        ax.set_title(axis_title, fontsize=18)
        fig.colorbar(m, ax=ax)
        ax.set_xticks(unique_xticks)
        # Set y-ticks
        ax.set_yticks(unique_xticks)
        if str_int_truelabel_mapper is not None:
            str_int_truelabel_mapper = {"UD" if k == 'undefined' else k: v for k, v in str_int_truelabel_mapper.items()}
            # overwrite ticks
            ax.set_yticklabels(list(str_int_truelabel_mapper.keys()), fontsize=15)
            ax.set_xticklabels(list(str_int_truelabel_mapper.keys()), fontsize=15)

        thresh = cm.max() / 2.
        for i, j in itertools.product(range(cm.shape[0]), range(cm.shape[1])):
            ax.text(j, i, "{:.2f}".format(cm[i, j]), fontsize=12,
                    horizontalalignment="center", color="white" if cm[i, j] > thresh else "black")

        ax.set_ylabel('true label', fontsize=18)
        ax.set_xlabel('predicted label', fontsize=18)
        plt.tight_layout()

        fig.savefig(os.path.join(self._cfg.results.save_folder, '{}_vs_{}_{}_Confusion_Matrix_{}_GEx.pdf'.format(
            observation_name, module_name, dataset, feature_selection_approach)))
        plt.close(fig=fig)

        return matched_obs_label_predicted, diag_score

    def get_top_n_accuracy(self, pred_labels, true_labels):
        # Predict labels
        matched_obs_label_predicted, df_prob_predicted_labels = self.predict(pred_labels=pred_labels)
        # Top n Accuracy # TODO try to understand this ..
        top_n_acc = []
        for n in range(len(self.probability_histogram.columns)):
            ncorrect = 0
            for p, pred in zip(true_labels, np.asarray(df_prob_predicted_labels)):
                indices = np.argsort(pred)[::-1]

                if p in indices[:n+1]:
                    ncorrect += 1

            top_n_acc.append(ncorrect / len(true_labels))
            print("Top {} Accuracy: ".format(n+1), ncorrect / len(true_labels))

        return top_n_acc

    def get_ch_db_scores(self, data, predicted_labels):
        try:
            ch_score = metrics.calinski_harabasz_score(X=data, labels=predicted_labels)
        except ValueError:
            ch_score = np.nan
        try:
            db_score = metrics.davies_bouldin_score(X=data, labels=predicted_labels)
        except ValueError:
            db_score = np.nan
        try:
            silhouette_score = metrics.silhouette_score(X=data, labels=predicted_labels)
        except ValueError:
            silhouette_score = np.nan

        try:
            # High density within a cluster, and low density between clusters indicates good clustering assignments.
            # score between -1 to 1, with the larger the value the better clustering solution
            # dbcv_scores = DBCV(X=data, labels=predicted_labels)
            dbcv_scores = hdbscan.validity_index(X=data, labels=predicted_labels)
        except ValueError:
            dbcv_scores = np.nan

        return ch_score, db_score, silhouette_score, dbcv_scores

