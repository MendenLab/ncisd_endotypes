import numpy as np
import pandas as pd
import scanpy as sc

from scripts.gene_selection_leiden.feature_selection import (HVGFeatureSelection, GeneProfileFeatureSelection,
                                                             StdFeatureSelection)
from sklearn.neighbors import NearestNeighbors
import igraph as ig
import leidenalg as la
import umap

import os

import matplotlib.pyplot as plt


class Leiden:
    def __init__(self, cfg, parts, normalisation, random_state):
        self._cfg = cfg
        self._parts = parts
        self.normalisation = normalisation
        self.resolution = self._cfg .cluster.resolution
        self.random_state = random_state
        self.partition_type = la.RBConfigurationVertexPartition

        self.leiden = None

        self.labels_ = None

    def fit(self, X_, coldata: pd.DataFrame, batch_keys: [str, list], gene_names=None):
        leiden_cluster = base_leiden(cfg=self._cfg, parts=self._parts, normalisation=self.normalisation,
                                     random_state=self.random_state,
                                     partition_type=self.partition_type, gene_names=gene_names)
        leiden_cluster.set_fit()
        labels = leiden_cluster.fit(X_, coldata=coldata, batch_keys=batch_keys)

        out = dict()
        out['labels'] = labels
        out['leiden'] = leiden_cluster
        self.leiden = out

    def predict(self, X: np.ndarray, coldata: pd.DataFrame, batch_keys: [str, list]):
        """ Merges cluster labels from hierarchical Leiden clustering into one array

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)

        Returns
        -------
        labels : ndarray of shape (n_samples,)

        """
        assert self.leiden is not None

        node = self.leiden
        node['leiden'].set_predict()
        # Transformations of data stored inside of Leiden model
        self.labels_ = node['leiden'].predict(X, coldata=coldata, batch_keys=batch_keys)

        return self.labels_


class base_leiden:
    def __init__(self, cfg, parts, normalisation, partition_type, random_state, gene_names=None):
        self._cfg = cfg
        self._parts = parts
        self.normalisation = normalisation
        self._random_state = random_state
        self._save_folder = self._cfg.results.save_folder
        self._resolution = self._cfg.cluster.resolution
        self._nn = self._cfg.cluster.nearest_neighbors
        self.partition_type = partition_type

        if self._cfg.feature_selection.do_feature_selection:
            if 'HVG' in self._cfg.feature_selection.method:
                self._feature_selection = HVGFeatureSelection(
                    self._cfg.feature_selection, random_state=self._random_state, save_folder=self._save_folder,
                    gene_names=gene_names)
            elif 'std' in self._cfg.feature_selection.method:
                self._feature_selection = StdFeatureSelection(
                    self._cfg.feature_selection, random_state=self._random_state, save_folder=self._save_folder,
                    gene_names=gene_names)
            else:
                self._feature_selection = GeneProfileFeatureSelection(
                    self._cfg.feature_selection, random_state=self._random_state, save_folder=self._save_folder,
                    gene_names=gene_names)

        else:
            self._feature_selection = None

        self._leiden = None
        self.graph = None

        self._is_fit = True

    def set_fit(self):
        self._is_fit = True

    def set_predict(self):
        self._is_fit = False

    def _get_igraph_from_adjacency(self, adjacency, directed=None):
        """ Get igraph graph from adjacency matrix. """

        sources, targets = adjacency.nonzero()
        weights = adjacency[sources, targets]
        if isinstance(weights, np.matrix):
            weights = weights.A1
        g = ig.Graph(directed=directed)
        g.add_vertices(adjacency.shape[0])  # this adds adjacency.shape[0] vertices
        g.add_edges(list(zip(sources, targets)))
        try:
            g.es['weight'] = weights
        except KeyError:
            pass
        if g.vcount() != adjacency.shape[0]:
            print(
                f'The constructed graph has only {g.vcount()} nodes. '
                'Your adjacency matrix contained redundant nodes.')
        return g

    def transform_data_to_knngraph(self, X_t: np.ndarray, coldata: pd.DataFrame, batch_keys: [str, list]):
        # Apply normalisation after feature selection
        X_t_normed = self.normalisation.normalise(
            X=X_t, normalisation_method=self._cfg.dataset.normalise, parts=self._parts,
            coldata=coldata, batch_keys=batch_keys)

        # Apply Leiden clustering
        # 1. First build KNN graph TODO replace with kernel
        neigh = NearestNeighbors(
            n_neighbors=1, radius=1.0, algorithm='auto', leaf_size=30, metric='minkowski',
            p=2, metric_params=None, n_jobs=None)
        neigh.fit(X_t_normed)
        knn_graph = neigh.kneighbors_graph(X_t_normed, n_neighbors=self._nn, mode='connectivity')
        # 2. Convert to igraph
        graph = self._get_igraph_from_adjacency(adjacency=knn_graph, directed=True)

        setattr(graph, "knn_graph", knn_graph)
        setattr(graph, "neighbors", neigh)

        return graph

    def fit(self, X: np.array, coldata: pd.DataFrame, batch_keys: [str, list]):
        # Apply feature selection
        if self._feature_selection is not None:
            # Apply normalisation before feature selection
            X_normed = self.normalisation.normalise(
                X=X, normalisation_method=self._cfg.dataset.normalise, parts=self._parts,
                coldata=coldata, batch_keys=batch_keys)
            # Identify gene profiles using normed counts as input
            self._feature_selection.fit(X_normed)
            # Select features on raw count matrix
            X_t = self._feature_selection.transform(X, fit=self._is_fit)
        else:
            X_t = X

        # Apply normalisation after feature selection, get KNN graph
        knn_graph = self.transform_data_to_knngraph(X_t=X_t, coldata=coldata, batch_keys=batch_keys)

        # Apply Leiden clustering until optimal clustering result is achieved
        self._leiden = la.find_partition(
            knn_graph, self.partition_type, resolution_parameter=self._resolution, seed=self._random_state,
            weights=np.array(knn_graph.es['weight']).astype(np.float64), n_iterations=-1)

        # Add clustering to graph
        setattr(knn_graph, "clustering", self._leiden)
        # save graph
        self.graph = knn_graph

        # Add 1 to start cluster names at 1 instead of 0
        return np.array(self._leiden.membership) + 1

    def predict(self, X, coldata: pd.DataFrame, batch_keys: [str, list]):
        assert self._leiden is not None

        # Apply feature selection
        if self._feature_selection is not None:
            X_t = self._feature_selection.transform(X)
        else:
            X_t = X

        # Apply normalisation after feature selection
        X_t_normed = self.normalisation.normalise(
            X=X_t, normalisation_method=self._cfg.dataset.normalise, parts=self._parts,
            coldata=coldata, batch_keys=batch_keys)

        # Apply matching of new sample to closest existing graph node
        dist, ind_similar_vertices = self.graph.neighbors.kneighbors(X_t_normed, n_neighbors=1, return_distance=True)
        clusters = np.asarray([self.graph.clustering.membership[i] for i in ind_similar_vertices.T[0]])

        # Add 1 to start cluster names at 1 instead of 0
        clusters = clusters + 1

        if self._is_fit:
            self.create_embedding(data=X_t_normed, label=clusters)

        return clusters

    def create_embedding(self, data, label, key=''):
        reducer = umap.UMAP(random_state=self._random_state, n_neighbors=self._nn)
        data_transformed = reducer.fit_transform(data)
        colors = sc.pl.palettes.default_20

        fig, ax = plt.subplots()
        for cluster_id in np.unique(label):
            inds = np.where(np.asarray(label) == cluster_id)[0]

            if len(inds) > 0:
                ax.plot(data_transformed[inds, 0], data_transformed[inds, 1], '.',
                        label=cluster_id, c=colors[cluster_id])
        ax.legend(title='Clusters', bbox_to_anchor=(1, 1.02), loc='upper left', ncol=1, prop={'size': 10},
                  frameon=False, title_fontsize=15)
        ax.set_ylabel('UMAP2', fontsize=18)
        ax.set_xlabel('UMAP1', fontsize=18)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        plt.tight_layout()
        fig.savefig(os.path.join(self._save_folder, 'UMAP_Clusters_{}.pdf'.format(key)), bbox_inches='tight')
        plt.close(fig=fig)
