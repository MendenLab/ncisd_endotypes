import numpy as np
import scanpy as sc
import anndata
from sklearn.cluster import KMeans
from itertools import accumulate
import operator

from matplotlib import pyplot as plt
import pandas as pd
import seaborn as sns

import os

plt.ion()


class base_feature_selection:
    def __init__(self, cfg, random_state, save_folder, gene_names=None):
        self._cfg = cfg

        self._means = None
        self._weights = None
        self._powers = None

        self.norm_type = cfg.normalization
        self.n_clusters = cfg.nclusters
        self.keep = cfg.genepercentage_keep  # cfg.clusters_keep
        self.sort_rows = cfg.sort_rows

        self._random_state = random_state

        self._save_folder = save_folder
        self.gene_names = gene_names

        self.global_indices = None

    def fit(self, X):
        raise NotImplemented

    def gene_normalise(self, X, sort_rows=True):
        # We normalize gene over the population
        out = X.copy()
        out_unsorted = X.copy()

        # TODO Not sorted gene profiles
        for idx, row in enumerate(X):
            gene_over_population = row.copy()

            if self.norm_type == 'mean':
                gene_over_population -= np.mean(gene_over_population)
            elif self.norm_type == 'scale':
                gene_over_population = gene_over_population - (
                        np.min(gene_over_population) / (np.max(gene_over_population) - np.min(gene_over_population)))
            elif self.norm_type == 'std':
                gene_over_population = (gene_over_population - np.mean(gene_over_population)) / np.std(
                    gene_over_population)

            if sort_rows:
                # By sorting the scaled expression values of a gene over all samples we remove the patient and
                # disease bias and ensure that we cluster genes independent of the samples but by their behavior
                # gene_over_population = np.sort(gene_over_population)
                unsorted_gene_over_population = np.copy(gene_over_population)
                sort_index = np.argsort(gene_over_population)
                gene_over_population = gene_over_population[sort_index]
            else:
                unsorted_gene_over_population = np.copy(gene_over_population)

            out[idx] = gene_over_population
            out_unsorted[idx] = unsorted_gene_over_population

        df = pd.DataFrame(out)
        df = df.iloc[[0, 1000, 17815], :].T
        df = np.array(df)
        df_genes = pd.DataFrame(data=df, columns=['intermediate', 'highly variable', 'constant'])
        df_genes = df_genes[['highly variable', 'intermediate', 'constant']]
        # Plot three examples of
        fig, ax = plt.subplots(figsize=(4, 4))
        sns.lineplot(data=df_genes, ax=ax, palette=['indianred', 'burlywood', 'dimgrey'], legend=True)
        sns.move_legend(ax, loc="best")
        # ax.get_legend()
        # ax.legend(labels=['intermediate', 'highly variable', 'constant'], loc='best', prop={'size': 14})
        # ax.get_legend().remove()
        ax.set_ylabel('scaled gene counts', fontsize=18)
        if sort_rows:
            ax.set_xlabel('sorted samples', fontsize=18)
        else:
            ax.set_xlabel('samples', fontsize=18)
        plt.setp(ax.get_xticklabels(), fontsize=14)
        plt.setp(ax.get_yticklabels(), fontsize=14)
        sns.despine(fig=fig, ax=ax)
        plt.savefig(os.path.join(
            self._save_folder, 'Sorted_{}_GeneProfiles_lineplot.pdf'.format(sort_rows)), bbox_inches='tight')
        plt.close(fig=fig)

        sns.set(font_scale=2, style='white')
        fig, ax = plt.subplots()
        sns.heatmap(df, yticklabels=False, xticklabels=False,
                    center=0, vmin=-5, vmax=5, cmap="RdBu_r", ax=ax,
                    cbar_kws={'label': 'normed gene expression'})
        ax.set_ylabel('genes', fontsize=18)
        ax.set_xlabel('samples')
        plt.savefig(os.path.join(
            self._save_folder, 'Sorted_{}_GeneProfiles_heatmap.pdf'.format(sort_rows)), bbox_inches='tight')
        plt.close(fig=fig)

        return out, out_unsorted

    def transform(self, X, fit=False):
        assert self._weights is not None, "Not initiated"

        if (self._cfg.method == 'HIG') | ('HVG' in self._cfg.method) | ('std' in self._cfg.method):
            # If we use percentage, we already have the global indices
            X_selected = X[:, self._weights]

            print('Fit: {}; Data shape after gene selection: {}'.format(fit, X_selected.shape))
            return X_selected
        else:
            X_weighted = np.zeros_like(X)

            for idx, (m, w, p) in enumerate(zip(self._means, self._weights, self._powers)):
                # TODO np.exp(-(X[:, idx] - m) ** 2 / (2 * stds ** 2))
                X_weighted[:, idx] = ((X[:, idx] - m) * (w ** (p + 1)))

            return X_weighted


class HVGFeatureSelection(base_feature_selection):
    def fit(self, X):
        X = X.copy().T

        n_genes = X.shape[0]
        n_top_genes = int(self.keep * n_genes / 100)

        #
        # For each row in X (aka genes), we will normalize them and sort their values from
        # smallest to the largest. This would give us the gene profiles!
        #

        # X_norm = self.gene_normalise(X=X, sort_rows=self.sort_rows)

        # std = np.std(X_norm, axis=0)
        # mean = np.mean(X_norm, axis=0)
        std = np.std(X, axis=0)
        mean = np.mean(X, axis=0)
        min_disp = np.min(std)
        max_disp = np.max(std)
        min_mean = np.min(mean)
        max_mean = np.max(mean)

        # mask = X_norm != 0
        # X_norm = X_norm[mask].reshape(X_norm.shape)

        # adata = anndata.AnnData(X=X_norm.T)
        adata = anndata.AnnData(X=X.T)
        sc.pp.highly_variable_genes(
            adata=adata, n_top_genes=n_top_genes, min_disp=min_disp, max_disp=max_disp,
            min_mean=min_mean, max_mean=max_mean, flavor='cell_ranger')

        # Returns a mask with true (HVG) and false
        self._weights = adata.var['highly_variable']


class StdFeatureSelection(base_feature_selection):
    def fit(self, X: np.ndarray):
        # transpose count matrix
        X = X.copy().T

        n_genes = X.shape[0]
        n_top_genes = int(self.keep * n_genes / 100)

        #
        # For each row in X (aka genes), we will normalize them and sort their values from
        # smallest to the largest. This would give us the gene profiles!
        #
        # normed_data = self.gene_normalise(X=X, sort_rows=self.sort_rows)
        # normed_data = stats.zscore(X, axis=1)
        # normed_data_df = pd.DataFrame(normed_data)
        normed_data_df = pd.DataFrame(X)

        # Compute the standard deviation along axis=1
        std_deviation = normed_data_df.std(axis=1)

        # Get the top n columns with the highest standard deviation
        top_columns = std_deviation.nlargest(n_top_genes).index.tolist()

        # self.global_indices = std_deviation.nlargest(X.shape[1]).index.tolist()

        # Returns a mask with true (biological relevant genes) and false
        self._weights = top_columns


class GeneProfileFeatureSelection(base_feature_selection):
    def cluster(self, X: np.ndarray, n_clusters: int):
        """

        Parameters
        ----------
        X : numpy.ndarray
        n_clusters : int

        Returns
        -------
        cluster_indices : list of gene profile clusters containing the indices of each cluster starting at 0

        """
        # Apply Kmeans on gene profiles
        kmeans = KMeans(n_clusters=n_clusters, random_state=self._random_state).fit(X)

        labels = kmeans.labels_

        cluster_indices = []

        # append indices of clusters starting at cluster 0
        for i in range(self.n_clusters):
            inds = np.where(labels == i)[0]
            cluster_indices.append(inds)

        return cluster_indices

    def reorder_clusters(self, X, cluster_indices):
        # We are using average std of the cluster members to rack the clusters

        stds = []
        means = []

        for indices in cluster_indices:
            std_gene_cluster = np.std(X[indices], axis=1)
            stds.append(np.mean(std_gene_cluster))
            means.append(np.mean(np.mean(X[indices], axis=1)))

        stds = np.array(stds)
        sorted_indices = np.argsort(stds)[::-1]  # Sorting std from the largest to smallest

        # sort cluster indices by sorted largest to smallest mean of std
        cluster_indices_sorted = [cluster_indices[i] for i in sorted_indices]

        # Re-ordering the elements of each cluster based on their std values
        cluster_indices_reordered = cluster_indices_sorted.copy()
        for idx, indices in enumerate(cluster_indices_sorted):
            cluster_stds = np.std(X[indices], axis=1)
            argsort = np.argsort(cluster_stds)[::-1]
            cluster_indices_reordered[idx] = indices[argsort]

        # Show sorted cluster profiles, highlight number of genes in cluster via color
        num_genes_cluster = [len(genes) for genes in cluster_indices_reordered]
        df_scatter = pd.DataFrame.from_dict(
            {'Cluster': np.arange(0, stds.shape[0]), 'Std': stds[sorted_indices], 'n_genes': num_genes_cluster})

        sns.set(font_scale=2, style='white')
        n_genes = X.shape[0]
        num_ele_lists = [len(v) for v in cluster_indices_reordered]
        cumsum_gene_list = self.cumulative_sum(num_ele_lists)
        perc_num_genes = int(self.keep * n_genes / 100)
        ind_closest_val = np.argmin([abs(val - perc_num_genes) for val in cumsum_gene_list])

        fig, ax = plt.subplots(figsize=(6, 4))
        sns.scatterplot(data=df_scatter, x='Cluster', y='Std', hue='n_genes', palette='RdBu', ax=ax)
        # 7% of genes == taking almost all genes from 18 Clusters into account for further analysis
        ax.axvline(x=ind_closest_val, c='k', ls='dashed', zorder=1)
        norm = plt.Normalize(df_scatter['n_genes'].min(), df_scatter['n_genes'].max())
        sm = plt.cm.ScalarMappable(cmap="RdBu", norm=norm)
        sm.set_array([])
        ax.set_xlabel('gene profile clusters', fontsize=18)
        ax.set_ylabel('standard deviation', fontsize=18)
        ax.set_xticks([0, ind_closest_val, self.n_clusters], [1, ind_closest_val + 1, self.n_clusters + 1])
        sns.despine(ax=ax)

        # Remove the legend and add a colorbar
        ax.get_legend().remove()
        ax.figure.colorbar(sm, label='number of genes')
        plt.tight_layout()
        sns.set(font_scale=1)
        plt.savefig(os.path.join(self._save_folder, 'cluster_vs_std.pdf'))
        plt.close('all')

        # plt.scatter(np.arange(0, stds.shape[0]), stds[sorted_indices], c=num_genes_cluster)

        # Plot expression profiles of first and last cluster in heatmap?
        df_first = pd.DataFrame(X[cluster_indices_reordered[0]])
        df_last = pd.DataFrame(X[cluster_indices_reordered[6]])
        df_first_last = pd.concat([df_first, df_last])

        # sns.set(font_scale=2, style='white')
        # fig, ax = plt.subplots()
        # sns.heatmap(df_first_last, yticklabels=False, xticklabels=False,
        #             center=0, vmin=-5, vmax=5, cmap="RdBu_r", ax=ax,
        #             cbar_kws={'label': 'normed gene expression'})
        # ax.set_xlabel('samples', fontsize=18)
        # ax.set_ylabel('genes', fontsize=18)
        # plt.savefig(os.path.join(self._save_folder, 'Clustered_GeneProfiles_first_sixth.pdf'), bbox_inches='tight')
        # plt.close(fig=fig)
        # sns.set(font_scale=1)

        # bool_similar = []
        # for indices, val in enumerate(cluster_indices_reordered):
        #     bool_similar.append(np.all(cluster_indices_sorted[indices] == cluster_indices_reordered[indices]))
        # np.all(bool_similar)  # True

        # for indices in cluster_indices_reordered:
        #     print(len(indices))

        global_indices = np.concatenate(cluster_indices_reordered)
        self.global_indices = cluster_indices_reordered

        # if we use number of gene profile clusters use this code
        # global_indices = np.arange(X.shape[0])

        # d_all = pd.DataFrame(X[global_indices], index=global_indices)
        # sns.set(font_scale=2, style='white')
        # fig, ax = plt.subplots()
        # sns.heatmap(d_all, yticklabels=False, xticklabels=False,
        #             center=0, vmin=-5, vmax=5, cmap="RdBu_r", ax=ax, cbar_kws={'label': 'normed gene expression'})
        # ax.set_xlabel('samples', fontsize=18)
        # ax.set_ylabel('genes', fontsize=18)
        # plt.savefig(os.path.join(self._save_folder, 'Sorted_Clustered_GeneProfiles_all.pdf'), bbox_inches='tight')
        # plt.close(fig=fig)

        return cluster_indices_reordered, global_indices

    def fit(self, X):
        # Transpose data in order to cluster genes
        X = X.copy().T

        #
        #
        # There are three steps
        #
        #

        #
        # For each row in X (aka genes), we will normalize them and sort their values from
        # smallest to the largest. This would give us the gene profiles!
        #

        X_norm, X_norm_unsorted = self.gene_normalise(X=X, sort_rows=self.sort_rows)

        #
        # We cluster the gene profiles using standard k-means
        #

        print('Number of gene profile clusters: {}'.format(self.n_clusters))
        cluster_indices = self.cluster(X_norm, n_clusters=self.n_clusters)

        #
        # Re-order the clusters based on how much variation their members have over the
        # population on average.
        #

        cluster_indices, global_indices = self.reorder_clusters(X, cluster_indices)

        # df = pd.DataFrame(X_norm[self._weights])
        # sns.heatmap(df, yticklabels=False, xticklabels=False, center=0, vmin=-5, vmax=5, cmap="RdBu_r")

        # Cumsum plot
        self.plot_cumsum(cluster_indices=cluster_indices, n_genes=X_norm.shape[0], save_folder=self._save_folder)

        # First cluster sorted gene profiles
        for cluster_number in np.arange(0, self.n_clusters, 1):
            save_folder_tmp = os.path.join(self._save_folder, 'clustered_gene_profiles')
            os.makedirs(save_folder_tmp, exist_ok=True)
            df_genes = pd.DataFrame(data=X_norm[cluster_indices[cluster_number]].T,
                                    columns=np.asarray(self.gene_names)[cluster_indices[cluster_number]])
            # Plot
            self.lineplot_gene_profile_clusters(df_genes, cluster_number, save_folder_tmp, sort_rows=self.sort_rows)
            # boxplot of clusters to show they don't have same std
            if cluster_number == 0:
                df_genes = df_genes.melt()
                self.boxplot(df_long=df_genes, save_folder=self._save_folder,
                             title='Cluster_{}_Sorted_{}_GeneProfiles_boxplot.pdf'.format(
                                 cluster_number + 1, self.sort_rows))

            df_genes = pd.DataFrame(data=X_norm_unsorted[cluster_indices[cluster_number]].T,
                                    columns=np.asarray(self.gene_names)[cluster_indices[cluster_number]])
            self.lineplot_gene_profile_clusters(df_genes, cluster_number, save_folder_tmp, sort_rows=False)

        # plot boxplot of expected highly variable, intermediate and constant genes
        constant_gene = ['TBP', 'GAPDH', 'SDHAF2', 'TSR3']  # definitely constant genes
        hvg = ['NOS2', 'IL36G', 'CCL27']
        # intermediate = ['S100A7', 'DEFB4A', 'LCE3A']
        d_all = pd.DataFrame(X_norm[global_indices], index=global_indices)

        df_boxplot = d_all.reset_index()
        df_boxplot.index = np.asarray(self.gene_names)[global_indices]
        # df_boxplot_preselected = df_boxplot.loc[:, df_boxplot.columns != 'index'].T[hvg + intermediate + constant_gene]
        df_boxplot_preselected = df_boxplot.loc[:, df_boxplot.columns != 'index'].T[hvg + constant_gene]
        df_boxplot_preselected = df_boxplot_preselected.melt()
        df_boxplot_preselected['genes'] = 'Unknown'
        df_boxplot_preselected.loc[df_boxplot_preselected['variable'].isin(hvg), 'genes'] = 'highly variable'
        # df_boxplot_preselected.loc[df_boxplot_preselected['variable'].isin(intermediate), 'genes'] = 'intermediate'
        df_boxplot_preselected.loc[df_boxplot_preselected['variable'].isin(constant_gene), 'genes'] = 'constant'
        self.boxplot(df_long=df_boxplot_preselected, save_folder=self._save_folder,
                     title='Preselected_GeneProfiles.pdf', hue='genes',
                     palette=['indianred', 'dimgrey'])  # ['indianred', 'burlywood', 'dimgrey']

        ix = np.concatenate([cluster_indices[0][:2],
                             cluster_indices[int(self.n_clusters/2)][:2],
                             cluster_indices[-1][:2]])
        df_boxplot_genes = df_boxplot.loc[df_boxplot['index'].isin(ix), df_boxplot.columns != 'index'].T
        df_boxplot_genes = df_boxplot_genes.melt()
        df_boxplot_genes['genes'] = 'Unknown'
        df_boxplot_genes.loc[df_boxplot_genes['variable'].isin(
            list(df_boxplot[df_boxplot['index'].isin(cluster_indices[0][:2])].index)), 'genes'] = 'highly variable'
        df_boxplot_genes.loc[df_boxplot_genes['variable'].isin(
            list(df_boxplot[df_boxplot['index'].isin(
                cluster_indices[int(self.n_clusters/2)][:2])].index)), 'genes'] = 'intermediate'
        df_boxplot_genes.loc[df_boxplot_genes['variable'].isin(
            list(df_boxplot[df_boxplot['index'].isin(cluster_indices[-1][:2])].index)), 'genes'] = 'constant'
        self.boxplot(df_long=df_boxplot_genes,
                     save_folder=self._save_folder,
                     title='Highly_variable_intermediate_constant_GeneProfiles.pdf', hue='genes',
                     palette=['indianred', 'burlywood', 'dimgrey'])

        #
        # Pick several of the top clusters as our feature selection
        #

        # Using number of clusters of gene profiles to read out informative genes
        # gene_indices = np.concatenate([cluster_indices[i] for i in self.keep])
        # self._weights = global_indices[gene_indices]
        # Using percentage of genes instead of merging clusters
        n_genes = X_norm.shape[0]
        print('Gene percentage to keep: {}'.format(self.keep))
        self._weights = global_indices[np.arange(0, int(self.keep * n_genes / 100))]

    @staticmethod
    def cumulative_sum(input_list):
        # Use the accumulate() function to perform a cumulative sum of the elements in the list
        cumulative_sum_iter = accumulate(input_list, operator.add)
        # Convert the iterator to a list and return it
        return list(cumulative_sum_iter)

    def plot_cumsum(self, cluster_indices, n_genes, save_folder):
        num_ele_lists = [len(v) for v in cluster_indices]
        cumsum_gene_list = self.cumulative_sum(num_ele_lists)

        perc_num_genes = int(self.keep * n_genes / 100)
        ind_closest_val = np.argmin([abs(val - perc_num_genes) for val in cumsum_gene_list])

        sns.set(font_scale=2, style='white')
        fig, ax = plt.subplots(figsize=(4, 4))
        ax.plot(cumsum_gene_list, color='black', lw=1)
        ax.axhline(y=int(self.keep * n_genes / 100), ls='dashed', lw=1, color='k')
        ax.axvline(x=ind_closest_val, ls='dashed', lw=1, color='k')
        plt.setp(ax.get_xticklabels(), fontsize=12)
        plt.setp(ax.get_yticklabels(), fontsize=12)
        ax.set_ylabel('Cumulative number of genes', fontsize=18)
        ax.set_xlabel('gene profile cluster', fontsize=18)
        ax.set_xticks([0, ind_closest_val, self.n_clusters], [1, ind_closest_val + 1, self.n_clusters + 1])
        sns.despine(fig=fig, ax=ax)
        plt.savefig(os.path.join(save_folder, 'CumSum_KMeans_plot.pdf'), bbox_inches='tight')
        plt.close(fig=fig)

    @staticmethod
    def lineplot_gene_profile_clusters(df_genes, cluster_number, save_folder_tmp, sort_rows):
        sns.set(font_scale=2, style='white')
        # Plot
        fig, ax = plt.subplots(figsize=(4, 4))
        for gene in df_genes.columns:
            ax.plot(df_genes[gene], label=gene, ls='-')
        if len(df_genes.columns) < 20:
            ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), prop={'size': 12}, ncol=3, frameon=False)
        ax.set_ylabel('scaled gene counts', fontsize=14)
        if sort_rows:
            ax.set_xlabel('sorted samples of cluster {}'.format(cluster_number + 1), fontsize=14)
        else:
            ax.set_xlabel('not sorted samples of cluster {}'.format(cluster_number + 1), fontsize=14)
        ax.set_title("number of genes: {}".format(df_genes.shape[1]), fontsize=12)
        plt.setp(ax.get_xticklabels(), fontsize=12)
        plt.setp(ax.get_yticklabels(), fontsize=12)
        sns.despine(fig=fig, ax=ax)
        plt.savefig(os.path.join(save_folder_tmp, 'Cluster_{}_Sorted_{}_GeneProfiles_lineplot.pdf'.format(
            cluster_number + 1, sort_rows)),
                    bbox_inches='tight')
        plt.close(fig=fig)

    @staticmethod
    def boxplot(df_long: pd.DataFrame, save_folder: str, title: str, hue: str = 'variable', palette: list = None):
        sns.set(font_scale=2, style='white')
        props = {'boxprops': {'edgecolor': 'k'}, 'medianprops': {'color': 'k'},
                 'whiskerprops': {'color': 'k'}, 'capprops': {'color': 'k'}}
        fig, ax = plt.subplots(figsize=(8, 4))
        sns.boxplot(df_long, x='variable', y='value', ax=ax, hue=hue, palette=palette, dodge=False, **props)
        ax.set_ylabel('scaled gene counts', fontsize=18)
        ax.set_xlabel('')
        plt.setp(ax.get_xticklabels(), fontsize=14)
        plt.setp(ax.get_yticklabels(), fontsize=14)
        if palette != None:
            ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False)
        else:
            ax.get_legend().remove()
        sns.despine(fig=fig)
        plt.savefig(os.path.join(save_folder, title), bbox_inches='tight')
        plt.close(fig=fig)
