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
import scanpy as sc

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

    def fit_predict(self):
        print("fit")
        # Predict subtypes using whole data set
        self._parts.set_application()
        # self._parts.set_train()
        X = self._parts.data

        self._cluster = Leiden(self._cfg, parts=self._parts, normalisation=self._normalisation,
                               random_state=self.random_state)
        self._cluster.fit(X, gene_names=self._dataset.genes, coldata=self._parts.observations(),
                          batch_keys=self._cfg.dataset.batch_keys)

        # 2. Predict subtypes using entire dataset
        predicted_labels = self._cluster.predict(X, coldata=self._parts.observations(),
                                                 batch_keys=self._cfg.dataset.batch_keys)

        return predicted_labels

    def eval_fit(self):
        print("Evaluate fit")
        # Predict subtypes using whole data set
        self._parts.set_application()
        # self._parts.set_train()
        X = self._parts.data
        true_labels = self._parts.pattern

        # 2. Predict subtypes using train set
        predicted_labels = self._cluster.predict(X, coldata=self._parts.observations(),
                                                 batch_keys=self._cfg.dataset.batch_keys)

        # Set up eval obj
        # Evaluate
        self._eval = Evaluation(cfg=self._cfg, cluster=self._cluster)

        # Metrics if true labels are not known
        ch_score, db_score, silhouette_score, dbcv_scores = self._eval.get_ch_db_scores(
            data=X, predicted_labels=predicted_labels)

        # Read out igraph object with membership
        igraph_obj = self._cluster.leiden[self._cfg.cluster.method]

        # 3. Plot evaluation
        # 3.1 Without feature selection
        # 3.1.1 Apply normalisation
        X_normed = self._normalisation.normalise(
            X=X, normalisation_method=self._cfg.dataset.normalise, parts=self._parts,
            coldata=self._parts.observations(), batch_keys=self._cfg.dataset.batch_keys)
        # 3.1.2 Get embedding
        data_reduced = tools.get_embedding(data=X_normed, random_state=self.random_state)
        module_name = '{}'.format(self._cfg.cluster.method)
        self._plot_control_plots(
            data_embedded=data_reduced, module_name=module_name, predicted_labels=predicted_labels,
            true_labels=true_labels)
        # 3.2 With feature selection
        # 3.2.1 Selected / weight genes
        X_transformed = self._cluster.leiden[self._cfg.cluster.method]._feature_selection.transform(X)
        # 3.2.2 Apply normalisation after feature selection
        X_t_normed = self._normalisation.normalise(
            X=X_transformed, normalisation_method=self._cfg.dataset.normalise, parts=self._parts,
            coldata=self._parts.observations(), batch_keys=self._cfg.dataset.batch_keys)
        # 3.2.3 Get embedding
        data_geneselection_reduced = tools.get_embedding(data=X_t_normed, random_state=self.random_state)
        module_name = '{} genes selected'.format(self._cfg.cluster.method)
        self._plot_control_plots(
            data_embedded=data_geneselection_reduced, module_name=module_name, predicted_labels=predicted_labels,
            true_labels=true_labels)

        global_indices = self._cluster.leiden[self._cfg.cluster.method]._feature_selection.global_indices

        return self._cluster.leiden[self._cfg.cluster.method]._feature_selection._weights, \
               db_score, ch_score, silhouette_score, igraph_obj, dbcv_scores, list(self._parts._dataset.sample_names), \
               data_geneselection_reduced, global_indices

    def _plot_control_plots(self, data_embedded, module_name, predicted_labels,
                            most_likely_true_label=None, true_labels=None):
        self.plot_umap(
            data=data_embedded, module_name=module_name, status=self._parts._model_status,
            observation_name="{} labels".format(self._cfg.cluster.method), label=predicted_labels)
        if most_likely_true_label is not None:
            self.plot_umap(
                data=data_embedded, module_name=module_name, status=self._parts._model_status,
                observation_name='{} clusters'.format(self._cfg.cluster.method), label=self._dataset.map_label_int_str(
                    int_labels=most_likely_true_label, mapping=self._dataset.pattern_mapping))
        if true_labels is not None:
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

        self.plot_umap_deltapga_scores(
            data=data_embedded, module_name=module_name, status=self._parts._model_status,
            observation_name="{} labels".format(self._cfg.cluster.method), label=predicted_labels, figsize=(6, 6))
        self.plot_umap(
            data=data_embedded, module_name=module_name, status=self._parts._model_status,
            observation_name='Diagnosis',
            label=self._dataset.map_label_int_str(
                int_labels=self._parts.diag, mapping=self._dataset.diag_mapping), color=self._parts.color_diag)
        self.plot_umap(
            data=data_embedded, module_name=module_name, status=self._parts._model_status,
            observation_name='Gender',
            label=list(self._parts.observations()['Sex_x'].values))
        self.plot_umap(
            data=data_embedded, module_name=module_name, status=self._parts._model_status,
            observation_name='batchID', label=list(self._parts.observations()['batchID'].values))

    def plot_umap(self, data, label, module_name, status, observation_name, color=None):
        fig, ax = plt.subplots(figsize=(8, 6))
        for ind_label, cluster_id in enumerate(np.unique(label)):
            inds = np.where(np.asarray(label) == cluster_id)[0]

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
        fig.savefig(os.path.join(self._cfg.results.save_folder, 'UMAP_{}_{}_{}.pdf'.format(
            observation_name, module_name, status)))
        plt.close(fig=fig)

    def plot_umap_deltapga_scores(
            self, data, observation_name,  module_name, status, label, color=None, figsize=(6, 6)):
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
                    ax.scatter(data[inds, 0], data[inds, 1], marker='.', label=cluster_id, c=color[inds])
                else:
                    ax.scatter(data[inds, 0], data[inds, 1], marker='.', label=cluster_id)
        ax.legend(title=observation_name, bbox_to_anchor=(1, 1.02), loc='upper left', ncol=1, prop={'size': 10},
                  frameon=False, title_fontsize=15)

        # non-responder
        ax.plot(data[:, 0][mask_nonresponder],
                data[:, 1][mask_nonresponder],
                marker="^", color='red', label='Non Responder', markersize=10, markeredgecolor="black",
                linestyle='none')
        # responder
        ax.plot(data[:, 0][mask_responder],
                data[:, 1][mask_responder],
                marker="^", color='gold', label='Responder', markersize=10, markeredgecolor="black", linestyle='none')
        # super responder
        ax.plot(data[:, 0][mask_superresponder],
                data[:, 1][mask_superresponder],
                marker="^", color='limegreen', label='Super Responder', markersize=10, markeredgecolor="black",
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


class ModuleOptimiseResolution:
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
        opt_param = m._cfg.feature_selection.genepercentage_keep
        m._cfg.results.save_folder = os.path.join(
            self._gridsearch_config.results.save_folder,
            "seed{}__optpercentage{}__resolution{}".format(params[0], opt_param, params[1]))
        os.makedirs(m._cfg.results.save_folder, exist_ok=True)

        m._cfg.cluster.resolution = params[1]
        print('Savefolder: ', m._cfg.results.save_folder)

        predicted_labels = m.fit_predict()
        selected_genes, db_score, ch_score, silhouette_score, igraph_obj, dbcv_scores, \
        muc_ids, data_geneselection_reduced, global_indices = m.eval_fit()

        print('=======================================\n')

        m._parts.set_train()
        train_ids = m._parts.sample_names()
        m._parts.set_test()
        test_ids = m._parts.sample_names()
        # Save params and results to dict
        scores = {'gene percentage to keep': opt_param, 'seed': params[0], 'resolution': params[1],
                  'Davies Bouldin index': db_score, 'Calinski-Harabasz index': ch_score,
                  'Silhouette coefficient': silhouette_score, 'Density Based clustering validation': dbcv_scores,
                  'index selected genes': selected_genes, 'MUC IDs': muc_ids, 'UMAP': data_geneselection_reduced,
                  'Gene names': list(m._dataset.genes[selected_genes]), 'Molecular subtypes': predicted_labels,
                  'train samples': train_ids, 'test samples': test_ids, 'all_genes': m._dataset.genes.tolist(),
                  'global_indices': global_indices}
                  # 'igraph object': igraph_obj}

        return scores

    def run_gridsearch(self):
        param_one = list(np.asarray(self._gridsearch_config.seeds))
        param_two = list(np.asarray(self._gridsearch_config.cluster.resolution))
        params_tuple = list(itertools.product(param_one, param_two))

        # Run job in parallel
        self.result = Parallel(n_jobs=self._gridsearch_config.njobs)(delayed(self.get_score)(p) for p in params_tuple)

        # Convert output list of dictionary into dictionary of list
        self.result = {key: [i[key] for i in self.result] for key in self.result[0]}

    # def save_parameters(self):
    #     igraph_obj = self.result.popitem()
    #     filename = os.path.join(
    #         self._gridsearch_config.results.save_folder, 'Graph_{}.graphml'.format(self.default_title))
    #     # TODO check how you can reconstruct graph
    #     igraph_obj[1][0].graph.write_graphml(filename)
    #
    #     # test
    #     adata = sc.read(os.path.join(igraph_obj[1][0]._cfg.dataset.path,
    #                                  'LESION_TNF_cohort_20231225_new_responder_definition.h5'))
    #     adata.obs['cohort'] = 'Michigan'
    #     adata.obs['Sex_x'] = adata.obs['Sex_x'].str.lower()
    #     # Subset to all measured genes in BRAIN cohort
    #     df_genes = pd.read_csv(os.path.join(igraph_obj[1][0]._cfg.dataset.path, 'LNL_Gene_information.csv'), index_col=0)
    #     mask_genes = adata.var.index.isin(df_genes.hgnc_symbol)
    #     mask_samples = adata.obs['sampleType'] == 'd'
    #     adata = adata[mask_samples, mask_genes].copy()
    #     # TODO Get raw counts - in my case stored in 'counts' but in your case it could be stored in adata.X instead
    #     X = adata.to_df(layer='counts')
    #     # add missing genes - need same dimension in test and train set for MinMaxScaler()
    #     # Find missing features in test set compared to the train set
    #     missing_features = set(df_genes.hgnc_symbol) - set(adata.var_names)  # 889
    #     # Add missing features to the test set with default values (e.g., zeros, mean of HKG)
    #     # replace with expression mean value of HKGs
    #     data_mean_hkg = X[['TBP', 'SDHAF2', 'TSR3', 'GAPDH']].mean(axis=1).astype(int)
    #     tmp = np.tile(np.array([data_mean_hkg]).transpose(), (1, len(missing_features)))
    #     df_tmp = pd.DataFrame(data=tmp, columns=list(missing_features), index=X.index)
    #     X = pd.concat([X, df_tmp], axis=1)  # Add missing feature with zeros (or any other default value)
    #     X = X[df_genes.hgnc_symbol]
    #     # Check if test set genes are ordered equally as train set genes -- important for normalisation and scaling
    #     assert np.all(X.columns == df_genes.hgnc_symbol), "Features don't align with Eyerich dataset"
    #
    #     # predict cluster
    #     clusters = igraph_obj[1][0].predict(np.asarray(X), coldata=adata.obs, batch_keys='Sex_x')
    #     self.result['Michigan_labels'] = clusters
    #
    #     # Integrate with Eyerich data set first
    #     X_eyerich = igraph_obj[1][0]._parts.data
    #     obs_eyerich = igraph_obj[1][0]._parts.observations()
    #     obs_eyerich['cohort'] = 'Eyerich'
    #     df_eyerich = pd.DataFrame(X_eyerich, columns=igraph_obj[1][0]._parts._dataset.genes, index=obs_eyerich.index)
    #     # combine eyerich and michigan
    #     df_combined = pd.concat([df_eyerich, X])
    #     df_obs_combined = pd.concat([obs_eyerich[['Sex_x', 'cohort']], adata.obs[['Sex_x', 'cohort']]])
    #     # predict cluster labels
    #     clusters_combined = igraph_obj[1][0].predict(
    #         np.asarray(df_combined), coldata=df_obs_combined, batch_keys=['Sex_x', 'cohort'])
    #     self.result['Eyerich_Michigan_labels'] = clusters_combined
    #
    #     filename = os.path.join(
    #         self._gridsearch_config.results.save_folder, 'Scores_{}.pkl'.format(self.default_title))
    #     # Save result to pickle file using pickle NOT pandas!
    #     with open(filename, 'wb') as ff:
    #         pickle.dump(self.result, ff)
    #
    #     # visualise together with Eyerich cohort
    #     X_eyerich = igraph_obj[1][0]._parts.data
    #     # Apply feature selection
    #     if igraph_obj[1][0]._feature_selection is not None:
    #         X_t_eyerich = igraph_obj[1][0]._feature_selection.transform(X_eyerich)
    #     else:
    #         X_t_eyerich = X_eyerich
    #
    #     # Apply feature selection
    #     if igraph_obj[1][0]._feature_selection is not None:
    #         X_t = igraph_obj[1][0]._feature_selection.transform(np.asarray(X))
    #     else:
    #         X_t = X
    #
    #     # Apply normalisation after feature selection
    #     X_t_normed = igraph_obj[1][0].normalisation.normalise(
    #         X=X_t, normalisation_method=igraph_obj[1][0]._cfg.dataset.normalise, parts=igraph_obj[1][0]._parts,
    #         coldata=adata.obs, batch_keys='Sex_x')
    #     igraph_obj[1][0].create_embedding(X_t_normed, label=clusters, key='Michigan')
    #
    #     # combine both cohorts
    #     X_t_combined = np.concatenate([X_t_eyerich, X_t])
    #     X_t_normed_combined = igraph_obj[1][0].normalisation.normalise(
    #         X=X_t_combined, normalisation_method=igraph_obj[1][0]._cfg.dataset.normalise, parts=igraph_obj[1][0]._parts,
    #         coldata=df_obs_combined, batch_keys=['Sex_x', 'cohort'])
    #     labels = [val + 1 for val in igraph_obj[1][0]._leiden.membership] + list(clusters)
    #     igraph_obj[1][0].create_embedding(X_t_normed_combined, label=labels, key='Eyerich_Michigan')

    def save_parameters(self):
        # igraph_obj = self.result.popitem()
        # igraph_obj = self.result.pop('igraph object', None)
        # filename = os.path.join(
        #     self._gridsearch_config.results.save_folder, 'Graph_{}.graphml'.format(self.default_title))
        # # TODO check how you can reconstruct graph
        # igraph_obj[1][0].graph.write_graphml(filename)
        #
        # # test
        # adata = sc.read(os.path.join(igraph_obj[1][0]._cfg.dataset.path,
        #                              'Pso_cohort_20240311_new_responder_definition.h5'))
        # adata.obs['batchID'] = '20'
        # adata.obs['cohort'] = 'Kiel'
        # adata.obs['Sex_x'] = adata.obs['Sex_x'].astype(str).replace({"1.0": "m", "2.0": 'f'}).str.lower()
        # # Subset to all measured genes in BRAIN cohort
        # df_genes = pd.read_csv(os.path.join(igraph_obj[1][0]._cfg.dataset.path, 'LNL_Gene_information.csv'), index_col=0)
        # mask_genes = adata.var.index.isin(df_genes.hgnc_symbol)
        # mask_samples = adata.obs['sampleType'] == 'd'
        # adata = adata[mask_samples, mask_genes].copy()
        # # TODO Get raw counts - in my case stored in 'counts' but in your case it could be stored in adata.X instead
        # X = adata.to_df(layer='counts')
        # # add missing genes - need same dimension in test and train set for MinMaxScaler()
        # # Find missing features in test set compared to the train set
        # missing_features = set(df_genes.hgnc_symbol) - set(adata.var_names)  # 889
        # # Add missing features to the test set with default values (e.g., zeros, mean of HKG)
        # # replace with expression mean value of HKGs
        # data_mean_hkg = X[['TBP', 'SDHAF2', 'TSR3', 'GAPDH']].mean(axis=1).astype(int)
        # tmp = np.tile(np.array([data_mean_hkg]).transpose(), (1, len(missing_features)))
        # df_tmp = pd.DataFrame(data=tmp, columns=list(missing_features), index=X.index)
        # X = pd.concat([X, df_tmp], axis=1)  # Add missing feature with zeros (or any other default value)
        # X = X[df_genes.hgnc_symbol]
        # # Check if test set genes are ordered equally as train set genes -- important for normalisation and scaling
        # assert np.all(X.columns == df_genes.hgnc_symbol), "Features don't align with Eyerich dataset"
        #
        # # predict cluster
        # clusters = igraph_obj[1][0].predict(np.asarray(X), coldata=adata.obs, batch_keys='Sex_x')
        # self.result['Kiel_labels'] = clusters
        #
        # # Integrate with Eyerich data set first
        # X_eyerich = igraph_obj[1][0]._parts.data
        # obs_eyerich = igraph_obj[1][0]._parts.observations()
        # obs_eyerich['cohort'] = 'Eyerich'
        # df_eyerich = pd.DataFrame(X_eyerich, columns=igraph_obj[1][0]._parts._dataset.genes, index=obs_eyerich.index)
        # # combine eyerich and michigan
        # df_combined = pd.concat([df_eyerich, X])
        # df_obs_combined = pd.concat([obs_eyerich[['Sex_x', 'cohort', 'batchID']],
        #                              adata.obs[['Sex_x', 'cohort', 'batchID']]])
        # # predict cluster labels
        # clusters_combined = igraph_obj[1][0].predict(
        #     np.asarray(df_combined), coldata=df_obs_combined, batch_keys=['Sex_x', 'batchID'])
        # self.result['Eyerich_Kiel_labels'] = clusters_combined

        filename = os.path.join(
            self._gridsearch_config.results.save_folder, 'Scores_{}.pkl'.format(self.default_title))
        # Save result to pickle file using pickle NOT pandas!
        with open(filename, 'wb') as ff:
            pickle.dump(self.result, ff)

        # # visualise together with Eyerich cohort
        # X_eyerich = igraph_obj[1][0]._parts.data
        # # Apply feature selection
        # if igraph_obj[1][0]._feature_selection is not None:
        #     X_t_eyerich = igraph_obj[1][0]._feature_selection.transform(X_eyerich)
        # else:
        #     X_t_eyerich = X_eyerich
        #
        # # Apply feature selection
        # if igraph_obj[1][0]._feature_selection is not None:
        #     X_t = igraph_obj[1][0]._feature_selection.transform(np.asarray(X))
        # else:
        #     X_t = X
        #
        # # Apply normalisation after feature selection
        # X_t_normed = igraph_obj[1][0].normalisation.normalise(
        #     X=X_t, normalisation_method=igraph_obj[1][0]._cfg.dataset.normalise, parts=igraph_obj[1][0]._parts,
        #     coldata=adata.obs, batch_keys=['Sex_x'])
        # igraph_obj[1][0].create_embedding(X_t_normed, label=clusters, key='Kiel')
        #
        # # combine both cohorts
        # X_t_combined = np.concatenate([X_t_eyerich, X_t])
        # X_t_normed_combined = igraph_obj[1][0].normalisation.normalise(
        #     X=X_t_combined, normalisation_method=igraph_obj[1][0]._cfg.dataset.normalise, parts=igraph_obj[1][0]._parts,
        #     coldata=df_obs_combined, batch_keys=['Sex_x', 'batchID'])
        # labels = [val + 1 for val in igraph_obj[1][0]._leiden.membership] + list(clusters)
        # igraph_obj[1][0].create_embedding(X_t_normed_combined, label=labels, key='Eyerich_Kiel')

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


def main():
    gs = ModuleOptimiseResolution()
    gs.run_gridsearch()
    gs.save_parameters()
    gs.save_config()


if __name__ == '__main__':
    main()
