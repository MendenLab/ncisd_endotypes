"""Plot diseases in UMAP

"""
from scripts.utils import add_colors

import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import os
import numpy as np
from datetime import date

fileformat = '.pdf'


def export_legend(legend, save_folder, ncol, filename="legend", expand=[-5, -5, 5, 5]):
    # https://stackoverflow.com/questions/4534480/get-legend-as-a-separate-picture-in-matplotlib
    fig = legend.figure
    sns.despine(left=True, bottom=True, fig=fig)
    fig.canvas.draw()
    # bbox = legend.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    bbox = legend.get_window_extent()
    bbox = bbox.from_extents(*(bbox.extents + np.array(expand)))
    bbox = bbox.transformed(fig.dpi_scale_trans.inverted())
    # fig.tight_layout()
    fig.savefig(os.path.join(save_folder, "{}_{}{}".format(filename, ncol, fileformat)),
                dpi='figure', bbox_inches=bbox)
    plt.close(fig=fig)


def main(save_folder):
    adata = sc.read(
        '/Users/christina.hillig/PycharmProjects/Eyerich_ERCGrant/Molecular_subtypes/input/data_freeze/Lesion_RNAseq_20210720_patient_meta_data_v04__CH_20220209_TherapyResponse_v04_v220715.h5')

    # Obervable to plot
    obs = 'diag'
    ncol = 3

    adata, _ = add_colors.diag_order_lesion(adata=adata)

    # Select only those diagnosis which are used in dendrogram
    diags = ['psoriasis', 'eczema', 'lichen planus', 'pityriasis ruba pilaris',
             'darier disease', 'lupus erythematosus', 'pyoderma gangrenosum']

    adata = adata[adata.obs['diag'].isin(diags)].copy()

    fig, ax = plt.subplots(figsize=(6, 6))
    sc.pl.umap(adata, color=obs, palette=adata.uns['{}_colors'.format(obs)], use_raw=False, size=160,
               show=False, ax=ax, title='', legend_loc='center right')

    # Put a legend below current axis
    # ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05), ncol=3, frameon=False, fontsize=16, markerscale=3.)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout()
    plt.savefig(os.path.join(save_folder, "UMAP_{}.pdf".format(obs)), bbox_inches='tight')
    plt.close(fig=fig)

    colordict = dict(zip(adata.obs[obs].cat.categories, adata.uns['{}_colors'.format(obs)]))
    f = lambda m, c: plt.plot([], [], marker='o', color=c, ls="none")[0]
    handles = [f("s", adata.uns['{}_colors'.format(obs)][i]) for i in range(len(adata.uns['{}_colors'.format(obs)]))]
    labels = list(colordict.keys())
    legend = plt.legend(handles, labels, loc=3, framealpha=1, frameon=False, ncol=ncol)
    export_legend(legend, save_folder=save_folder, ncol=ncol)


if __name__ == '__main__':
    output_dir = os.path.join('/Users/christina.hillig/PycharmProjects/Eyerich_ERCGrant/Molecular_subtypes/output',
                              'Figure_3B_dendrogram_legend', str(date.today()))
    os.makedirs(output_dir, exist_ok=True)

    main(save_folder=output_dir)
