from scripts.utils.init import AnyRandom
import numpy as np
import pandas as pd
import pickle
import configparser

from scripts.gene_selection_leiden.dataset import dataset
from scripts.gene_selection_leiden.partitions import partitions
from scripts.gene_selection_leiden.config import config, GridSearchParameters
from scripts.gene_selection_leiden import tools

from scripts.gene_selection_leiden.base_leiden import Leiden
from scripts.gene_selection_leiden.evaluation import Evaluation
from scripts.gene_selection_leiden.normalisation import Normalisation

import os
from matplotlib import pyplot as plt

import itertools
from joblib import Parallel, delayed

plt.ion()


class Module:
    def __init__(self, random_state: AnyRandom = 0):
        self._cfg = config()
        self.random_state = random_state
        self._dataset = dataset(self._cfg.dataset)
        self._parts = partitions(self._cfg.part, self._dataset, self.random_state)

        self._normalisation = Normalisation()

        self._cluster = None
        self._eval = None

    def train(self):
        print("Train")
        self._parts.set_train()
        X = self._parts.data

        self._cluster = Leiden(self._cfg, parts=self._parts, normalisation=self._normalisation,
                               random_state=self.random_state)
        self._cluster.fit(X, gene_names=self._dataset.genes, coldata=self._parts.observations(),
                          batch_keys=self._cfg.dataset.batch_keys)

    def eval_train(self):
        print("Evaluate Train")
        self._parts.set_train()

        # Read out data
        X = self._parts.data
        true_labels = self._parts.pattern
        module_name = '{}'.format(self._cfg.cluster.method)

        # 2. Predict subtypes using train set
        predicted_labels = self._cluster.predict(X, coldata=self._parts.observations(),
                                                 batch_keys=self._cfg.dataset.batch_keys)

        # Evaluate
        self._eval = Evaluation(cfg=self._cfg, cluster=self._cluster)
        self._eval.fit(pred_labels=predicted_labels, true_labels=true_labels)
        self._eval.plot_probability_histogram(
            true_labels=true_labels, module_name=module_name,
            dataset=self._parts._model_status, str_int_truelabel_mapper=self._dataset.pattern_mapping,
            colors=self._parts.color_pattern)
        most_likely_true_label, _ = self._eval.get_confusion_matrix(
            pred_labels=predicted_labels, true_labels=self._parts.pattern,
            feature_selection_approach=self._cfg.feature_selection.method, observation_name='Pattern',
            module_name=module_name, normalize=True, dataset=self._parts._model_status,
            str_int_truelabel_mapper=self._dataset.pattern_mapping)

        # Apply normalisation after feature selection
        X_normed = self._normalisation.normalise(
            X=X, normalisation_method=self._cfg.dataset.normalise, parts=self._parts,
            coldata=self._parts.observations(), batch_keys=self._cfg.dataset.batch_keys)
        # Plot evaluation
        data_reduced = tools.get_embedding(data=X_normed, random_state=self.random_state)

        predicted_test_labels = self._cluster.predict(
            X, coldata=self._parts.observations(), batch_keys=self._cfg.dataset.batch_keys)
        self._plot_control_plots(
            data_embedded=data_reduced, module_name=module_name, true_labels=true_labels,
            predicted_labels=predicted_test_labels, most_likely_true_label=most_likely_true_label)

        # Create umap using only selected genes / weighted genes
        X_transformed = self._cluster.leiden[self._cfg.cluster.method]._feature_selection.transform(X)

        # Apply normalisation after feature selection
        X_t_normed = self._normalisation.normalise(
            X=X_transformed, normalisation_method=self._cfg.dataset.normalise, parts=self._parts,
            coldata=self._parts.observations(), batch_keys=self._cfg.dataset.batch_keys)

        data_geneselection_reduced = tools.get_embedding(data=X_t_normed, random_state=self.random_state)
        module_name = '{} genes selected'.format(self._cfg.cluster.method)
        self._plot_control_plots(
            data_embedded=data_geneselection_reduced, module_name=module_name, true_labels=true_labels,
            predicted_labels=predicted_test_labels, most_likely_true_label=most_likely_true_label)

    def score(self):
        print("Test")
        self._parts.set_test()
        X = self._parts.data
        # External criterion: Compare predicted labels to true labels
        true_labels = self._parts.pattern
        predicted_test_labels = self._cluster.predict(
            X, coldata=self._parts.observations(), batch_keys=self._cfg.dataset.batch_keys)

        module_name = '{}'.format(self._cfg.cluster.method)

        print(predicted_test_labels)

        # Internal criterion: Metrics if true labels are not known
        ch_score, db_score, silhouette_score, dbcv_scores = self._eval.get_ch_db_scores(
            data=X, predicted_labels=predicted_test_labels)

        self._eval.plot_probability_histogram(
            true_labels=self._parts.pattern, module_name=module_name, colors=self._parts.color_pattern,
            dataset=self._parts._model_status, str_int_truelabel_mapper=self._dataset.pattern_mapping)

        # Metric if true labels are known
        most_likely_true_label, pattern_score = self._eval.get_confusion_matrix(
            pred_labels=predicted_test_labels, true_labels=self._parts.pattern,
            feature_selection_approach=self._cfg.feature_selection.method, observation_name='Pattern',
            module_name=self._cfg.cluster.method, normalize=True, dataset=self._parts._model_status,
            str_int_truelabel_mapper=self._dataset.pattern_mapping)

        # Apply normalisation on raw counts
        X_normed = self._normalisation.normalise(
            X=X, normalisation_method=self._cfg.dataset.normalise, parts=self._parts,
            coldata=self._parts.observations(), batch_keys=self._cfg.dataset.batch_keys)

        data_reduced = tools.get_embedding(data=X_normed, random_state=self.random_state)
        module_name = '{}'.format(self._cfg.cluster.method)
        self._plot_control_plots(
            data_embedded=data_reduced, module_name=module_name, true_labels=true_labels,
            predicted_labels=predicted_test_labels, most_likely_true_label=most_likely_true_label)

        # Create umap using only selected genes / weights
        X_transformed = self._cluster.leiden[self._cfg.cluster.method]._feature_selection.transform(X)
        # Apply normalisation after feature selection
        X_t_normed = self._normalisation.normalise(
            X=X_transformed, normalisation_method=self._cfg.dataset.normalise, parts=self._parts,
            coldata=self._parts.observations(), batch_keys=self._cfg.dataset.batch_keys)

        data_geneselection_reduced = tools.get_embedding(data=X_t_normed, random_state=self.random_state)
        module_name = '{} genes selected'.format(self._cfg.cluster.method)
        self._plot_control_plots(
            data_embedded=data_geneselection_reduced, module_name=module_name, true_labels=true_labels,
            predicted_labels=predicted_test_labels, most_likely_true_label=most_likely_true_label)

        return most_likely_true_label, pattern_score, \
               self._cluster.leiden[self._cfg.cluster.method]._feature_selection._weights, db_score, ch_score

    def _plot_control_plots(self, data_embedded, module_name, true_labels, predicted_labels, most_likely_true_label):
        self.plot_umap(
            data=data_embedded, module_name=module_name, status=self._parts._model_status,
            observation_name="{} labels".format(self._cfg.cluster.method), label=predicted_labels)
        self.plot_umap(
            data=data_embedded, module_name=module_name, status=self._parts._model_status,
            observation_name='{} clusters'.format(self._cfg.cluster.method), label=self._dataset.map_label_int_str(
                int_labels=most_likely_true_label, mapping=self._dataset.pattern_mapping))
        self.plot_umap(
            data=data_embedded, module_name=module_name, status=self._parts._model_status,
            observation_name='Pattern',
            label=self._dataset.map_label_int_str(int_labels=true_labels, mapping=self._dataset.pattern_mapping),
            color=self._parts.color_pattern)
        self.plot_umap_deltapga_scores(
            data=data_embedded, module_name=module_name, status=self._parts._model_status,
            observation_name='Pattern', label=self._dataset.map_label_int_str(
                int_labels=true_labels, mapping=self._dataset.pattern_mapping),
            color=self._parts.color_pattern, figsize=(6, 6))
        if 'Diagnosis' in self._parts.observations().columns:
            self.plot_umap(
                data=data_embedded, module_name=module_name, status=self._parts._model_status,
                observation_name='Diagnosis',
                label=self._dataset.map_label_int_str(
                    int_labels=self._parts.diag, mapping=self._dataset.diag_mapping), color=self._parts.color_diag)
        if 'Sex_x' in self._parts.observations().columns:
            self.plot_umap(
                data=data_embedded, module_name=module_name, status=self._parts._model_status,
                observation_name='Gender',
                label=self._parts.observations()['Sex_x'].to_list())
        self.plot_umap(
            data=data_embedded, module_name=module_name, status=self._parts._model_status,
            observation_name='batchID', label=self._parts.observations()['batchID'].to_list())

    def plot_umap(self, data, label, module_name, status, observation_name, color=None):
        fig, ax = plt.subplots(figsize=(8, 6))
        for ind_label, cluster_id in enumerate(np.unique(label)):
            inds = np.where(np.asarray(label) == cluster_id)[0]

            if len(inds) > 0:
                if color is not None:
                    ax.scatter(data[inds, 0], data[inds, 1], marker='.', label=cluster_id, c=color[inds], s=160)
                else:
                    ax.scatter(data[inds, 0], data[inds, 1], marker='.', label=cluster_id, s=160)
        ax.legend(title=observation_name, bbox_to_anchor=(1, 1.02), loc='upper left', ncol=1, prop={'size': 15},
                  frameon=False, title_fontsize=17)
        ax.set_ylabel('UMAP2', fontsize=18)
        ax.set_xlabel('UMAP1', fontsize=18)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        plt.tight_layout()
        fig.savefig(os.path.join(self._cfg.results.save_folder, 'UMAP_{}_{}_{}.pdf'.format(
            observation_name, module_name, status)))
        plt.close(fig=fig)

    def plot_umap_deltapga_scores(
            self, data, observation_name,  module_name, status, label, color=None, figsize=(6, 6)):
        if ('Response week 12' in self._parts.observations().columns) & (
                'Therapeutic target' in self._parts.observations().columns):
            mask_nonresponder = (self._parts.observations()['Response week 12'] == 'Non-responder') & \
                                (self._parts.observations()['Therapeutic target'] == 'IL-23')
            mask_responder = (self._parts.observations()['Response week 12'] == 'Responder') & \
                             (self._parts.observations()['Therapeutic target'] == 'IL-23')
            mask_superresponder = (self._parts.observations()['Response week 12'] == 'Super-responder') & \
                                  (self._parts.observations()['Therapeutic target'] == 'IL-23')

            fig, ax = plt.subplots(figsize=figsize)
            for ind_label, cluster_id in enumerate(np.unique(label)):
                inds = np.where(np.asarray(label) == cluster_id)[0]

                if len(inds) > 0:
                    if color is not None:
                        ax.scatter(data[inds, 0], data[inds, 1], marker='.', label=cluster_id, c=color[inds], s=160)
                    else:
                        ax.scatter(data[inds, 0], data[inds, 1], marker='.', label=cluster_id, s=160)
            ax.legend(title=observation_name, bbox_to_anchor=(1, 1.02), loc='upper left', ncol=1, prop={'size': 15},
                      frameon=False, title_fontsize=17)

            # non-responder
            ax.plot(data[:, 0][mask_nonresponder],
                    data[:, 1][mask_nonresponder],
                    marker="^", color='red', label='Non Responder', markersize=50, markeredgecolor="black",
                    linestyle='none')
            # responder
            ax.plot(data[:, 0][mask_responder],
                    data[:, 1][mask_responder],
                    marker="^", color='gold', label='Responder', markersize=50, markeredgecolor="black", linestyle='none')
            # super responder
            ax.plot(data[:, 0][mask_superresponder],
                    data[:, 1][mask_superresponder],
                    marker="^", color='limegreen', label='Super Responder', markersize=50, markeredgecolor="black",
                    linestyle='none')

            ax.set_ylabel('UMAP2', fontsize=18)
            ax.set_xlabel('UMAP1', fontsize=18)
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            ax.spines['left'].set_visible(True)
            ax.spines['bottom'].set_visible(True)
            plt.tight_layout()
            fig.savefig(os.path.join(
                self._cfg.results.save_folder,
                'UMAP_Response week 12_{}_{}_{}.pdf'.format(observation_name, module_name, status)))
            plt.close(fig=fig)


class ModuleGridsearch:
    def __init__(self):
        self._gridsearch_config = GridSearchParameters()
        os.makedirs(self._gridsearch_config.results.save_folder, exist_ok=True)
        self.result = None

        self.min_genepercentage = np.min(self._gridsearch_config.feature_selection.genepercentage_keep)
        self.max_genepercentage = np.max(self._gridsearch_config.feature_selection.genepercentage_keep)

        self.min_resolution = np.min(self._gridsearch_config.cluster.resolution)
        self.max_resolution = np.max(self._gridsearch_config.cluster.resolution)

        # Names for file titles
        resolution_ = "minResolution{}_maxResolution{}".format(self.min_resolution, self.max_resolution)
        feature_percentage = "minPercentage{}_maxPercentage{}".format(self.min_genepercentage, self.max_genepercentage)
        nseeds_ = len(self._gridsearch_config.seeds)

        self.default_title = "GridSearch_{}_{}_{}_".format(nseeds_, feature_percentage, resolution_)

    def get_score(self, params):
        m = Module(random_state=params[0])
        m._cfg.results.save_folder = os.path.join(
            self._gridsearch_config.results.save_folder,
            "seed{}__percentage{}__resolution{}".format(params[0], params[1], params[2]))

        # pkl file name
        filename = os.path.join(
            m._cfg.results.save_folder, 'Scores_{}.pkl'.format(self.default_title))

        if os.path.isfile(filename):
            # load result
            with open(filename, 'rb') as ff:
                scores = pickle.load(ff)
        else:
            os.makedirs(m._cfg.results.save_folder, exist_ok=True)

            m._cfg.feature_selection.genepercentage_keep = params[1]
            m._cfg.cluster.resolution = params[2]
            print('Savefolder: ', m._cfg.results.save_folder)

            m.train()
            m.eval_train()
            most_likely_true_label, pattern_score, selected_genes, db_score, ch_score = m.score()

            print('=======================================\n')

            m._parts.set_train()
            train_ids = m._parts.sample_names()
            m._parts.set_test()
            test_ids = m._parts.sample_names()

            # Save params and results to dict
            scores = {'gene percentage to keep': params[1], 'seed': params[0], 'resolution': params[2],
                      'predicted_label': most_likely_true_label, 'acc': pattern_score,
                      'index genes': selected_genes, 'Davies Bouldin Score': db_score,
                      'Calinski-Harabasz Index': ch_score, 'train samples': train_ids, 'test samples': test_ids}

            # Save result to pickle file using pickle NOT pandas!
            with open(filename, 'wb') as ff:
                pickle.dump(scores, ff)

        return scores

    def run_gridsearch(self):
        param_one = list(np.asarray(self._gridsearch_config.seeds))
        param_two = list(np.asarray(self._gridsearch_config.feature_selection.genepercentage_keep))
        param_three = list(np.asarray(self._gridsearch_config.cluster.resolution))
        params_tuple = list(itertools.product(param_one, param_two, param_three))

        # Run job in parallel
        self.result = Parallel(n_jobs=self._gridsearch_config.njobs)(delayed(self.get_score)(p) for p in params_tuple)

        # Convert output list of dictionary into dictionary of list
        # self.result = {key: [i[key] for i in self.result] for key in self.result[0]}

    def get_mean_std_acc_dbscore(self):
        result = {key: [i[key] for i in self.result] for key in self.result[0]}
        # Get overall acc
        mean_seed_acc = []
        st_seed_acc = []
        mean_db_score = []
        st_db_score = []

        df = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in result.items()]))
        for resolution_val in np.unique((result['resolution'])):
            ind_clusters = np.where(np.asarray(result['resolution']) == resolution_val)[0]
            # Get Mean and std  of Davies Boulding score over all gene percentage to keep
            nclusters_db_score = list(map(result['Davies Bouldin Score'].__getitem__, ind_clusters))
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

        self.result['mean_overall_acc'] = mean_seed_acc
        self.result['std_overall_acc'] = st_seed_acc
        self.result['mean_Davies Bouldin Score'] = mean_db_score
        self.result['std_Davies Bouldin Score'] = st_db_score

    def save_parameters(self):
        # Convert output list of dictionary into dictionary of list
        self.result = {key: [i[key] for i in self.result] for key in self.result[0]}

        filename = os.path.join(
            self._gridsearch_config.results.save_folder,
            'Scores_{}.pkl'.format(self.default_title))

        # Save result to pickle file using pickle NOT pandas!
        with open(filename, 'wb') as ff:
            pickle.dump(self.result, ff)

    def save_config(self):
        cfg = configparser.ConfigParser()

        # cfg.read('config.ini')
        dict_gs = self._gridsearch_config.__dict__
        cfg.add_section('dataset')
        cfg['dataset'] = dict_gs['dataset']

        cfg.add_section('part')
        cfg['part'] = dict_gs['part']

        cfg.add_section('feature_selection')
        cfg.set('feature_selection', 'do_feature_selection', str(dict_gs['feature_selection']['do_feature_selection']))
        cfg.set('feature_selection', 'normalization', dict_gs['feature_selection']['normalization'])
        cfg.set('feature_selection', 'nclusters', str(dict_gs['feature_selection']['nclusters']))
        cfg.set('feature_selection', 'genepercentage_keep', str(dict_gs['feature_selection']['genepercentage_keep']))
        cfg.set('feature_selection', 'sort_rows', str(dict_gs['feature_selection']['sort_rows']))
        cfg.set('feature_selection', 'method', dict_gs['feature_selection']['method'])

        cfg.add_section('seeds')
        cfg.set('seeds', 'values', str(dict_gs['seeds']))

        cfg.add_section('cluster')
        cfg['cluster'] = dict_gs['cluster']

        cfg.add_section('results')
        cfg['results'] = dict_gs['results']

        filename = os.path.join(
            self._gridsearch_config.results.save_folder, 'config_{}.ini'.format(self.default_title))

        with open(filename, 'w') as ff:
            cfg.write(ff)

    def plot_errorbars(self):
        num_unique_percentage = len(np.unique(self.result['gene percentage to keep']))
        # per cluster
        counter = 0
        for ind in range(0, len(np.unique(np.asarray(self._gridsearch_config.cluster.resolution)))):
            fig, ax = plt.subplots()
            ax.errorbar(np.unique(self.result['gene percentage to keep']),
                        self.result['mean_overall_acc'][counter:counter+num_unique_percentage],
                        yerr=self.result['std_overall_acc'][counter:counter+num_unique_percentage], fmt='-o')
            ax.set_ylabel('Overall Accuracy')
            ax.set_xlabel('gene percentage to keep')
            fig.savefig(os.path.join(
                self._gridsearch_config.results.save_folder,
                'OverallAcc_{}_Resolution_{}.png'.format(self.default_title, np.asarray(
                    self._gridsearch_config.cluster.resolution)[ind])))
            plt.close(fig=fig)

            counter += num_unique_percentage

    def plot_percentage_vs_profiles(self):
        # Convert data to array
        scores_array = np.reshape(np.asarray(self.result['mean_overall_acc']),
                                  newshape=(
                                  np.asarray(self._gridsearch_config.feature_selection.genepercentage_keep).shape[0],
                                  np.asarray(self._gridsearch_config.cluster.resolution).shape[0]))
        scores_array = scores_array.transpose()

        # Plot gene percentage to keep vs number of gene profiles
        fig, ax = plt.subplots(figsize=(6, 6))
        if 1 in scores_array.shape:
            # CARE: x and y tick labels dont match parameters
            im = ax.imshow(scores_array, cmap=plt.cm.Reds, interpolation='none')
        else:
            im = ax.imshow(scores_array, cmap=plt.cm.Reds, interpolation='none',
                           extent=[self.min_genepercentage, self.max_genepercentage,
                                   self.max_resolution, self.min_resolution])
        ax.set_aspect(2)
        fig.colorbar(im, ax=ax)
        fig.savefig(os.path.join(
            self._gridsearch_config.results.save_folder,
            'FeaturesPerc_vs_Resolution_{}.png'.format(self.default_title)))
        plt.close(fig=fig)

    def plot_davies_bouldin_score(self):
        fig, ax = plt.subplots()
        ax.errorbar(np.asarray(self._gridsearch_config.cluster.resolution), self.result['mean_Davies Bouldin Score'],
                    yerr=self.result['std_Davies Bouldin Score'], fmt='-o')
        ax.set_xlabel('{} resolution'.format(str(self._gridsearch_config.cluster.method).capitalize()))
        ax.set_ylabel('Davies Bouldin score')

        fig.savefig(os.path.join(
            self._gridsearch_config.results.save_folder,
            'DaviesBouldinScore_{}.png'.format(self.default_title)))
        plt.close(fig=fig)


def main():
    gs = ModuleGridsearch()
    gs.run_gridsearch()
    gs.get_mean_std_acc_dbscore()
    gs.save_parameters()
    gs.save_config()
    gs.plot_errorbars()
    gs.plot_percentage_vs_profiles()
    gs.plot_davies_bouldin_score()


if __name__ == '__main__':
    main()
