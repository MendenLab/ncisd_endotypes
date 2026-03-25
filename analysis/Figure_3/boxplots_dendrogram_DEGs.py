from scripts.utils import add_endotypes, normalisation

import numpy as np
import scanpy as sc
import pandas as pd
from itertools import cycle
from datetime import date

import os
import glob

from matplotlib import pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages, FigureCanvasPdf
from statannotations.Annotator import Annotator


fileformat = '.pdf'


def load_deg_lists(input_path, log2fc_cut, padj_cut):
    full_path_1vs1 = [os.path.join(input_path, f) for f in
                      glob.glob(os.path.join(input_path, "*.csv"))]

    degs_list = []
    # df_dict = {}
    degs_dict = {}
    degs_dict_df = {}
    for file in full_path_1vs1:
        splitted = file.split(os.sep)[-1].split('__')[0:3]
        splitted[0] = str.replace(splitted[0], "DEGs_", "")
        comparison = "_".join(splitted)
        df_tmp = pd.read_csv(file)

        # Mark all DEGs
        mask = ((df_tmp.log2FoldChange > log2fc_cut) | (df_tmp.log2FoldChange < -log2fc_cut)) & (
                df_tmp.padj < padj_cut)
        df_tmp = df_tmp.loc[mask, :]

        # Save to dictionary
        # df_dict = df_tmp[['log2FoldChange', 'padj', 'hgnc_symbol']]
        degs_list.extend(list(df_tmp['hgnc_symbol']))
        degs_dict[comparison] = list(df_tmp['hgnc_symbol'])
        degs_dict_df[comparison] = df_tmp[['hgnc_symbol', 'log2FoldChange', 'padj']]

    return degs_list, degs_dict, degs_dict_df


def main(save_folder, input_dir, input_dir_degs_dendrogram):
    adata = sc.read(
        os.path.join(
            input_dir,
            'Lesion_RNAseq_20210720_patient_meta_data_v04__CH_20220209_TherapyResponse_v04_v220715_TherapyResponseCorrected__Endotypes_230620.h5'))
    # adata = add_endotypes.load_molecular_subtypes(adata=adata, data_root=input_dir)

    obs_name = 'Endotypes'

    log2fc_cut = 0.5
    padj_cut = 0.05
    min_classes = 1

    # Create Heatmap with samples on x-axis and DEGs on y axis, apply also Hierarchical clustering

    # 1. Load all DEGs for 1 vs 1 paired comparison
    degs_list, degs_dict, degs_dict_df = load_deg_lists(
        input_path=input_dir_degs_dendrogram, log2fc_cut=log2fc_cut, padj_cut=padj_cut)

    # normalise data
    bulk_data_temp = adata.copy()
    df_gex = bulk_data_temp.to_df(layer='counts')
    # Apply normalisation and batch correction
    array_gex = normalisation.edger_normalise(
        X=np.asarray(df_gex), coldata=bulk_data_temp.obs, logcpm=True, batch_keys=['batchID', 'Sex_x'])

    palette = {'E9_ vs _E10': ['lightcoral', 'indianred'], 'E7_ vs _E8_E9_E10': ['firebrick', 'tomato'],
               'E3_ vs _E4': ['papayawhip', 'navajowhite'], 'E1_E2_ vs _E3_E4': ['sandybrown', 'peachpuff'],
               'E11_ vs _E12_E13': ['dodgerblue', 'steelblue'], 'E5_E6_ vs _E7_E8_E9_E10': ['coral', 'orangered'],
               'E5_E6_E7_E8_E9_E10_ vs _E11_E12_E13_E1_E2_E3_E4': ['orangered', 'chocolate'],
               'E11_E12_E13_ vs _E1_E2_E3_E4': ['royalblue', 'orange'], 'E1_ vs _E2': ['blanchedalmond', 'tan'],
               'E12_ vs _E13': ['lightblue', 'lightskyblue'], 'E5_ vs _E6': ['lightsalmon', 'mistyrose'],
               'E8_ vs _E9_E10': ['salmon', 'crimson'],
               'E1_ vs _E6': ["#006ddb", "#490092"]}

    # New: 23.01.2025
    genes_to_plot = {
        'E5_E6_E7_E8_E9_E10_ vs _E11_E12_E13_E1_E2_E3_E4': ['FLG2', 'KRT2', 'KRT6C', 'PI3'],
        'E11_E12_E13_ vs _E1_E2_E3_E4': ['S100A7A', 'SERPINB4', 'FN1', 'CXCL9'],
        'E1_E2_ vs _E3_E4': ['CXCL8', 'MMP1', 'KRT2', 'CXCL9']}

    # OLD:
    # genes_to_plot = {
    #     'E5_E6_E7_E8_E9_E10_ vs _E11_E12_E13_E1_E2_E3_E4': ['IL1B', 'CXCL8', 'CD300E',  'KRT26', 'CD274', 'CCL3',
    #                                                         'SOX5', 'ZNF273'],
    #     'E11_E12_E13_ vs _E1_E2_E3_E4': ['VCAM1', 'CXCL8', 'LCN2', 'NOS2', 'PLA2G4D'],
    #     'E1_E2_ vs _E3_E4': ['IL1B', 'IL6', 'IFNG', 'CARD11', 'IL37', 'CXCL8', 'HIF1A', 'MMP1', 'MMP3', 'CXCL1',
    #                          'CSF2', 'IL1RL1',  'CARD9', 'CASP14']}

    rename_comparisons = {'E1 E2 E3 E4': 'E1-4', 'E11 E12 E13': 'E11-13', 'E5 E6 E7 E8 E9 E10': 'E5-E10',
                          'E11 E12 E13 E1 E2 E3 E4': 'E1-4/E11-13'}

    # palette = {'9_ vs _10': ["#920000", "#E69F00"], '7_ vs _8_9_10': ["#b66dff", 'tomato'],
    #            '3_ vs _4': ["#004949", "#009292"], '1_2_ vs _3_4': ['sandybrown', 'peachpuff'],
    #            '11_ vs _12_13': ["#D55E00", 'steelblue'], '5_6_ vs _7_8_9_10': ['coral', 'orangered'],
    #            '5_6_7_8_9_10_ vs _11_12_13_1_2_3_4': ['orangered', 'chocolate'],
    #            '11_12_13_ vs _1_2_3_4': ['royalblue', 'orange'], '1_ vs _2': ["#006ddb", "#b6dbff"],
    #            '12_ vs _13': ["#8B4513", "#999999"], '5_ vs _6': ["#ff6db6", "#490092"],
    #            '8_ vs _9_10': ["#000000", 'crimson']}

    for comparison in degs_dict_df.keys():
        save_folder_comp = os.path.join(save_folder, comparison)
        os.makedirs(save_folder_comp, exist_ok=True)

        # group_1 = ['E' + s for s in comparison.split('vs')[0].split('_')[:-1]]
        # group_2 = ['E' + s for s in comparison.split('vs')[1].split('_')[1:]]
        group_1 = comparison.split('vs')[0].split('_')[:-1]
        group_2 = comparison.split('vs')[1].split('_')[1:]

        df_tmp = degs_dict_df[comparison]
        df_tmp.to_excel(os.path.join(save_folder_comp, 'DEGs_{}_vs_{}_Age_Sex_batchID_subtypes.xlsx'.format(
                "_".join(group_1), "_".join(group_2))), index=False)
        # Sort dataframe by log2fc and padj
        # df_tmp = df_tmp.sort_values(['log2FoldChange', 'padj'], ascending=[False, True])
        df_tmp = df_tmp.sort_values(['padj'], ascending=[True])

        # Boxplot of top genes of comparison
        df = pd.DataFrame(array_gex, columns=adata.var_names, index=adata.obs.index)
        if comparison in genes_to_plot.keys():
            genes_filtered = genes_to_plot[comparison]
        else:
            # Select top 10 genes and last 10
            df_tmp_top = df_tmp.iloc[:10, :]
            df_tmp_last = df_tmp.tail(10)
            df_tmp_top_last = pd.concat([df_tmp_top, df_tmp_last])
            df_tmp_top_last.loc[:, 'log10padj'] = -np.log10(df_tmp_top_last['padj'])

            genes = list(df_tmp_top_last.hgnc_symbol)
            # Read out genes which can be found also in lesion adata object
            genes_filtered = list(adata.var_names[adata.var_names.isin(genes)])
        df = df.loc[:, genes_filtered]
        # Group by pairs and remove rest
        df['cluster_id'] = bulk_data_temp.obs[obs_name]
        df = df.loc[df.cluster_id.isin(group_1) | df.cluster_id.isin(group_2), :]
        df['cluster_id'] = df['cluster_id'].astype(str)
        if " ".join(group_1) in rename_comparisons:
            name_group_1 = rename_comparisons[" ".join(group_1)]
        else:
            name_group_1 = " ".join(group_1)
        if " ".join(group_2) in rename_comparisons:
            name_group_2 = rename_comparisons[" ".join(group_2)]
        else:
            name_group_2 = " ".join(group_2)
        df.loc[df.cluster_id.isin(group_1), 'cluster_id'] = name_group_1
        df.loc[df.cluster_id.isin(group_2), 'cluster_id'] = name_group_2
        df['cluster_id'] = df['cluster_id'].astype('category')

        col_names = ['cluster_id'] + genes_filtered
        long_df = df[col_names].melt(id_vars=['cluster_id'], var_name='gene', value_name='value')
        # Reorder categories to match DEG comparison
        long_df['cluster_id'] = long_df['cluster_id'].cat.reorder_categories([name_group_1, name_group_2])

        # hatches = ["o", "o", "xx", "xx"]  # ["//", "xx", "o", "*", "\\"]

        # colors = ["#006ddb", "#b6dbff", "#004949", "#009292", "#ff6db6", "#490092",
        #           "#b66dff", "#000000", "#920000", "#E69F00", "#D55E00", "#8B4513", "#999999"]

        if comparison in genes_to_plot.keys():
            genes = long_df['gene'].unique().tolist()
            pairs = [((gene, long_df['cluster_id'].cat.categories.to_list()[0]),
                      (gene, long_df['cluster_id'].cat.categories.to_list()[1])) for gene in genes]

            annotations = []
            for gene in genes:
                if df_tmp['hgnc_symbol'].isin([gene]).sum() == 1:
                    annotations.append("padj={:.2e}\nlog2FC={:.2f}".format(
                    df_tmp.loc[df_tmp['hgnc_symbol'] == gene, 'padj'].values[0],
                    df_tmp.loc[df_tmp['hgnc_symbol'] == gene, 'log2FoldChange'].values[0]))
                else:
                    annotations.append("ns")

            fig, ax = plt.subplots(figsize=(6, 3))
            sns.boxplot(data=long_df, x='gene', y='value', ax=ax, hue='cluster_id', palette=palette[comparison],
                        showfliers=False)
            sns.stripplot(data=long_df, x='gene', y='value', hue='cluster_id',
                          palette='dark:black', ax=ax, dodge=True, ec='k', linewidth=1, size=2, legend=False)
            # add statsvalues (padj + log2FC)
            annot = Annotator(ax, pairs=pairs, data=long_df, x='gene', y='value', hue='cluster_id')
            annot.set_custom_annotations(annotations)
            annot.annotate()
            ax.set_ylabel('normed counts', fontsize=16)
            ax.set_xlabel('', fontsize=16)
            plt.setp(ax.get_xticklabels(), rotation=0, fontsize=14)
            plt.setp(ax.get_yticklabels(), fontsize=14)
            # # Put a legend to the right of the current axis
            ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1), title="", ncol=2, prop={'size': 13},
                      frameon=False, title_fontsize=15)
            sns.despine(fig=fig, ax=ax, trim=False)
            plt.tight_layout()
            fig.savefig(os.path.join(save_folder_comp, '{}_vs_{}_boxplot.pdf'.format(
                "_".join(group_1), "_".join(group_2))), bbox_inches='tight')
            plt.close(fig=fig)

            with PdfPages(os.path.join(save_folder_comp, '{}_vs_{}_genes_boxplot.pdf'.format(
                    "_".join(group_1), "_".join(group_2)))) as pages:
                for ind, gene in enumerate(genes_filtered):
                    pair = [((gene, long_df['cluster_id'].cat.categories.to_list()[0]),
                             (gene, long_df['cluster_id'].cat.categories.to_list()[1]))]
                    if df_tmp['hgnc_symbol'].isin([gene]).sum() == 1:
                        annotation = ["padj={:.2e}\nlog2FC={:.2f}".format(
                            df_tmp.loc[df_tmp['hgnc_symbol'] == gene, 'padj'].values[0],
                            df_tmp.loc[df_tmp['hgnc_symbol'] == gene, 'log2FoldChange'].values[0])]
                    else:
                        annotation = ["ns"]

                    fig, ax = plt.subplots(figsize=(3, 3))
                    sns.boxplot(data=long_df.loc[long_df['gene'] == gene], x='gene', y='value', hue='cluster_id',
                                palette=palette[comparison], ax=ax, showfliers=False)
                    sns.swarmplot(data=long_df.loc[long_df['gene'] == gene], x='gene', y='value', hue='cluster_id',
                                  palette='dark:black', ax=ax, dodge=True, legend=False, size=3)
                    # add statsvalues (padj + log2FC)
                    annot = Annotator(ax=ax, pairs=pair, data=long_df.loc[long_df['gene'] == gene],
                                      x='gene', y='value', hue='cluster_id')
                    annot.set_custom_annotations(annotation)
                    annot.annotate()
                    ax.set_title(gene, fontsize=18)

                    # # iterate through the patches for each subplot
                    # for i, patch in enumerate(ax.patches):
                    #     # Set a different hatch for each bar
                    #     patch.set_hatch(hatches[i])
                    #     patch.set_edgecolor('k')

                    # l = ax.legend()
                    # for lp, hatch in zip(l.get_patches(), hatches):
                    #     lp.set_hatch(hatch)
                    #     lp.set_edgecolor('k')

                    ax.set_ylabel('normed counts', fontsize=18)
                    ax.set_xlabel('', fontsize=18)
                    ax.set_xticks([])
                    plt.setp(ax.get_yticklabels(), fontsize=14)
                    # # Put a legend to the right of the current axis
                    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1), title="",
                              ncol=2, prop={'size': 13}, frameon=False, title_fontsize=15)
                    sns.despine(fig=fig, ax=ax, trim=False)
                    plt.tight_layout()

                    fig.savefig(os.path.join(save_folder_comp, '{}_vs_{}_{}_boxplot.pdf'.format(
                        "_".join(group_1), "_".join(group_2), gene)), bbox_inches='tight')

                    canvas = FigureCanvasPdf(fig)
                    canvas.print_figure(pages)

                    plt.close(fig=fig)


if __name__ == '__main__':
    general_dir = os.path.join(
        '/Volumes', 'CH__data', 'Projects', 'Eyerich_AG_projects', 'BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer')
    output_dir = os.path.join(general_dir, 'analysis', 'Molecular_subtypes', 'output',
                              'Boxplots_Dendrogram_DEGs', str(date.today()))
    os.makedirs(output_dir, exist_ok=True)

    # data_root = '/Users/christina.hillig/PycharmProjects/Eyerich_ERCGrant/Molecular_subtypes/input/data_freeze'
    data_root = os.path.join(general_dir, 'analysis', 'Molecular_subtypes', 'input', 'h5_files', 'LESION')
    # degs_dendorgram = os.path.join(general_dir, 'analysis', 'Molecular_subtypes',
    #                                'output/DGE_analysis/Endotypes_Dendrogram/2023-06-22')
    degs_dendorgram = os.path.join(general_dir, 'analysis', 'Molecular_subtypes',
                                   'output/DGE_analysis/Endotypes_Dendrogram/2024-01-17')

    main(save_folder=output_dir, input_dir=data_root, input_dir_degs_dendrogram=degs_dendorgram)
