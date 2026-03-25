from scripts.utils import add_colors, add_endotypes, normalisation
from scripts.feature_engineering import scaling

import numpy as np
import scanpy as sc
import pandas as pd
from datetime import date

from scipy.spatial.distance import pdist
import scipy.stats as stats
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import dendrogram, linkage

import os

from matplotlib import pyplot as plt
import seaborn as sns

import rpy2.robjects as ro
import rpy2.robjects.packages as rpackages
from rpy2.robjects import numpy2ri, pandas2ri

numpy2ri.activate()
pandas2ri.activate()

base = rpackages.importr("base")
ape = rpackages.importr("ape")
factoextra = rpackages.importr("factoextra")
dplyr = rpackages.importr('dplyr')
ggtree = rpackages.importr('ggtree')
ggplot2 = rpackages.importr('ggplot2')

fileformat = '.pdf'


def main(save_folder, input_dir):
    adata = sc.read(
        os.path.join(input_dir,
                     'Lesion_RNAseq_20210720_patient_meta_data_v04__CH_20220209_TherapyResponse_v04_v220715.h5'))
    adata = add_endotypes.load_molecular_subtypes(adata=adata, data_root=input_dir)

    # Sample: PseudoID 10350, MUC4259
    # new diag: lichenoid drug reaction
    # new pattern: pattern 1
    adata.obs.loc[adata.obs['diag'] == 'undefined', 'Pattern'] = '1'
    adata.obs['diag'] = adata.obs['diag'].astype(str)
    adata.obs.loc[adata.obs['diag'] == 'undefined', 'diag'] = 'lichenoid drug reaction'
    adata.obs['diag'] = adata.obs['diag'].astype('category')

    adata, _ = add_colors.diag_order_lesion(adata=adata)
    adata, _ = add_colors.sdiag_order_lesion(adata=adata)
    adata, _ = add_colors.pattern_order_lesion(adata=adata)
    adata, _ = add_colors.endotype_order_lesion(adata=adata, obs_name='Molecular Subtype res0.9')

    merge_methods = ['mean', 'median']
    # 'hamming': for binary data
    metric_methods = ['braycurtis', 'canberra', 'chebyshev', 'cityblock',
        'correlation', 'cosine', 'dice', 'euclidean',
        'jaccard', 'jensenshannon', 'kulczynski1',
        'mahalanobis', 'matching', 'minkowski', 'rogerstanimoto',
        'russellrao', 'seuclidean', 'sokalmichener', 'sokalsneath',
        'sqeuclidean', 'yule']
    linkage_methods = ['complete', 'average', 'weighted', 'ward', 'single']

    # Read out dataframe
    # Should we use only our 7% of genes?
    df = adata.to_df(layer='counts')[list(adata.var_names[adata.var['geneselection']])]
    # Apply normalisation and batch correction
    array_gex = normalisation.edger_normalise(X=np.asarray(df), coldata=adata.obs, batch_keys=['batchID', 'Sex_x'])
    df_normed = pd.DataFrame(array_gex, columns=adata.var_names[adata.var['geneselection']],
                             index=adata.obs.index)
    # Scale data
    df_normed_zscores = stats.zscore(df_normed, axis=0)
    adata.obs['Molecular Subtype res0.9'] = adata.obs['Molecular Subtype res0.9'].astype(int)
    df_normed_zscores['cluster'] = adata.obs['Molecular Subtype res0.9'] + 1
    df_normed_zscores['cluster'] = df_normed_zscores['cluster'].astype(str)

    df_normed_zscores['cluster'] = df_normed_zscores['cluster'].replace(
        {'1': 'E12', '2': 'E7', '3': 'E13', "4": 'E3', "5": 'E1', "6": 'E9', "7": 'E8',
         "8": 'E10', "9": 'E11', "10": 'E5', "11": 'E6', "12": 'E2', "13": 'E4'})

    for merge in merge_methods:
        # Merge samples of a cluster into one column by using mean, median, ect
        if merge == 'mean':
            df_merged = df_normed_zscores.groupby(['cluster']).mean()  # Used in manuscript
        elif merge == 'median':
            df_merged = df_normed_zscores.groupby(['cluster']).median()

        df_merged.to_csv(os.path.join(save_folder, 'Normed_Dataframe_{}_gene_selection.csv'.format(merge)))
        # for linkage_method in linkage_methods:
        #     for metric_method in metric_methods:
                # # Merge samples of a cluster into one column by using mean, median, ect
                # if merge == 'mean':
                #     df_merged = df_normed_zscores.groupby(['cluster']).mean()
                # elif merge == 'median':
                #     df_merged = df_normed_zscores.groupby(['cluster']).median()
                #
                # df_merged.to_csv(os.path.join(save_folder, 'Normed_Dataframe_{}_gene_selection.csv'.format(merge)))

                # # Build distance matrix uisng e.g. euclidean distance
                # # dist_matrix = distance_matrix(X1,X1)
                # try:
                #     X = pdist(df_merged, metric=metric_method)  # squareform(pdist(df_merged, metric=metric_method))
                # except ValueError:
                #     continue
                #
                # # apply linkage using ward
                # try:
                #     Z = linkage(X, method=linkage_method, metric=metric_method)
                # except ValueError:
                #     continue

                # ro.r.assign('save_folder', save_folder)
                # ro.r.assign('df', df_merged)
                # ro.r.assign("metric_method",  metric_method)
                # ro.r.assign("linkage_method", linkage_method)
                # ro.r("res.dist <- dist(df, method=metric_method)")
                # ro.r("Z <- hclust(d=res.dist, method=linkage_method)")
                # ro.r("p1 <- base::plot(ape::as.phylo(Z), type = 'fan', label.offset = 1, cex = 0.7)")
                # # ro.r("pdf(file=file.path(save_folder, paste0('Dendrogram', '.pdf')), width=6, height=6)")
                # # ro.r("print(p1)")
                # ro.r("dev.copy(pdf, file.path(save_folder, paste0('Dendrogram', '.pdf')))")
                # ro.r("dev.off()")
                #
                # ro.r("circ <- factoextra::fviz_dend(Z, cex = 0.8, lwd = 0.8, k = 4, rect = TRUE, k_colors = c('darkorange', 'orange', 'blue', 'red'), rect_border = 'gray', rect_fill = TRUE, type = 'circular')")
                # ro.r("print(circ)")
                # ro.r("dev.copy(pdf, file.path(save_folder, paste0('Circular_Dendrogram', '.pdf')))")
                # ro.r("dev.off()")
                #
                # ro.r(
                #     "rectangular <- factoextra::fviz_dend(Z, cex = 0.8, lwd = 0.8, k = 4, rect = TRUE, k_colors = c('darkorange', 'orange', 'blue', 'red'), rect_border = 'gray', rect_fill = TRUE)")
                # ro.r("print(rectangular)")
                # ro.r("dev.copy(pdf, file.path(save_folder, paste0('Rectangular_Dendrogram', '.pdf')))")
                # ro.r("dev.off()")
                #
                # # With ggtree
                # # Tutorial https://yulab-smu.top/treedata-book/chapter5.html
                # ro.r("p.tree <- ggtree::ggtree(ape::as.phylo(Z), layout='fan') + geom_hilight( node=15, fill='darkorange', alpha=.6) + geom_hilight(node=17, fill='red', alpha=.6) + geom_hilight(node=18, fill='blue', alpha=.6) + geom_tiplab() + geom_strip('13', '12', barsize=1, color='orange', label='Pattern 1 like',  offset.text=11, offset=10, angle = -60) + geom_strip('8', '7', barsize=1, color='blue', label = 'Pattern 3 like', offset.text=11, offset=10, angle = 0) + geom_strip('11', '10', barsize=1, color='red', label = 'Pattern 2a like', offset.text=11, offset=10, angle = 60) ")
                # ro.r("print(p.tree)")
                # ro.r("dev.copy(pdf, file.path(save_folder, paste0('Circular_Dendrogram_tree', '.pdf')))")
                # ro.r("dev.off()")

                # # Plot dendorgram (round)
                # fig, ax = plt.subplots(figsize=(4, 3))
                # dn = dendrogram(Z, ax=ax)
                # sns.despine(fig=fig, ax=ax)
                # plt.savefig(os.path.join(save_folder, 'Dendrogram_{}_{}_{}.pdf'.format(
                #     merge, metric_method, linkage_method)))
                # plt.close(fig=fig)


if __name__ == '__main__':
    output_dir = os.path.join('/Volumes/CH__data/Projects/Eyerich_AG_projects',
                              'BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer/analysis/output',
                              'Dendrogram', str(date.today()))
    os.makedirs(output_dir, exist_ok=True)

    data_root = os.path.join('/Volumes/CH__data/Projects/Eyerich_AG_projects',
                             'BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer/data/h5_files/LESION')

    main(save_folder=output_dir, input_dir=data_root)
