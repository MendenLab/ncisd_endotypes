from scripts.utils.init import AnyRandom
import numpy as np

from scripts.gene_selection_HC.dataset import dataset
from scripts.gene_selection_HC.partitions import partitions
from scripts.gene_selection_HC.config import config, GridSearchParameters
from scripts.gene_selection_HC import tools
import pickle

from scripts.gene_selection_HC.h_kmeans import h_kmeans
from scripts.gene_selection_HC.evaluation import Evaluation
from scripts.gene_selection_HC.normalisation import Normalisation

import os
from matplotlib import pyplot as plt
plt.ion()


class Module:
    def __init__( self, random_state: AnyRandom = 0):
        self._cfg = config()
        self.random_state = random_state
        self._dataset = dataset( self._cfg.dataset )
        self._parts = partitions( self._cfg.part, self._dataset, self.random_state )

        self._normalisation = Normalisation(cfg=self._cfg, parts=self._parts)

        self._kmeans = None
        self._eval = None

    def train( self ):
        self._parts.set_train()
        X = self._parts.data

        self._kmeans = h_kmeans( self._cfg , normalisation=self._normalisation, random_state=self.random_state)
        self._kmeans.fit(X)

    def eval_train(self):
        self._parts.set_train()

        # Read out data
        X = self._parts.data
        true_labels = self._parts.pattern

        # Evaluate
        self._eval = Evaluation(cfg=self._cfg, kmeans=self._kmeans)
        self._eval.fit(data=X, true_labels=true_labels)
        self._eval.plot_probability_histogram(
            true_labels=true_labels, module_name='HKMeans',
            dataset='Train', str_int_truelabel_mapper=self._dataset.pattern_mapping)
        most_likely_true_label, _ = self._eval.get_confusion_matrix(
            data=X, true_labels=self._parts.pattern,
            feature_selection_approach=self._cfg.feature_selection.method, observation_name='Pattern',
            module_name='H_KMeans', normalize=True, dataset='Train',
            str_int_truelabel_mapper=self._dataset.pattern_mapping)

        # Plot evaluation
        print(most_likely_true_label)
        data_reduced = tools.get_embedding(data=X, random_state=self.random_state)
        self.plot_umap(
            data=data_reduced, module_name='H_KMeans', status='Train', observation_name='KMeans clusters',
            label=self._dataset.map_label_int_str(
                int_labels=most_likely_true_label, mapping=self._dataset.pattern_mapping))
        self.plot_umap(
            data=data_reduced, module_name='H_KMeans', status='Train', observation_name='Pattern',
            label=self._dataset.map_label_int_str(int_labels=true_labels, mapping=self._dataset.pattern_mapping),
            color=self._parts.color_pattern)

        predicted_test_labels = self._kmeans.predict(X)
        self.plot_umap(
            data=data_reduced, module_name='H_KMeans', status='Train', observation_name='KMeans labels',
            label=predicted_test_labels)

    def save( self ):
        with open('./scripts/hierarchical_kmeans/save_state.pkl', 'wb') as ff :
            save_data = {}
            save_data['kmeans'] = self._kmeans
            pickle.dump( save_data, ff )

    def load( self ):
        with open('./scripts/hierarchical_kmeans/save_state.pkl', 'rb') as ff :
            save_data = pickle.load( ff )
            self._kmeans = save_data['kmeans']

    def plot_umap(self, data, label, module_name, status, observation_name, color=None):
        fig, ax = plt.subplots()
        for ind_label, cluster_id in enumerate(np.unique(label)):
            inds = np.where( np.asarray(label) == cluster_id )[0]

            if len(inds) > 0:
                if color is not None:
                    ax.scatter(data[inds, 0], data[inds, 1], marker='.', label=cluster_id, c=color[inds])
                else:
                    ax.scatter(data[inds, 0], data[inds, 1], marker='.', label=cluster_id)
        ax.legend(title=observation_name, bbox_to_anchor=(1, 1.02), loc='upper left', ncol=1, prop={'size': 10},
                  frameon=False, title_fontsize=15)
        ax.set_ylabel('UMAP2', fontsize=18)
        ax.set_xlabel('UMAP1', fontsize=18)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        plt.tight_layout()
        fig.savefig(os.path.join(self._cfg.results.save_folder, 'UMAP_{}_{}_{}.png'.format(
            observation_name, module_name, status)))
        plt.close(fig=fig)

    def do_stuff( self ):
        self._parts.set_test()
        X = self._parts.data
        predicted_test_labels = self._kmeans.predict(X)

        print(predicted_test_labels)

        most_likely_true_label, diag_score = self._eval.get_confusion_matrix(
            data=X, true_labels=self._parts.pattern,
            feature_selection_approach=self._cfg.feature_selection.method, observation_name='Pattern',
            module_name='H_KMeans', normalize=True, dataset='Test',
            str_int_truelabel_mapper=None)

        data_reduced = tools.get_embedding(data=X, random_state=self.random_state)
        self.plot_umap(
            data=data_reduced, label=self._dataset.map_label_int_str(
                int_labels=most_likely_true_label, mapping=self._dataset.pattern_mapping),
            module_name='H_KMeans', status='Test', observation_name='Predicted clusters')

        return most_likely_true_label, diag_score


# def ensembl_labels():
#     most_likely_true_labels = []
#     diag_scores = []
#     for seed in self._cfg.seeds:
#         m = Module(random_state=seed)
#         m.train()
#         m.eval_train()
#         most_likely_true_label, diag_score = m.do_stuff()
#         most_likely_true_labels.append(most_likely_true_label)
#         diag_scores.append(diag_score)


def main():
    m = Module(random_state=0)
    m.train()
    print("evaluate training")
    m.eval_train()
    m.do_stuff()

    # TODO add normalisation

    # ensembl_labels()


if __name__ == '__main__':
    main()
