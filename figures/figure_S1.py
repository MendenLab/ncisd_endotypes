"""Plot diseases in UMAP

"""
from scripts.utils import add_colors

import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import os
import numpy as np
from datetime import date
import pandas as pd

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
    fig.savefig(os.path.join(save_folder, "Figure_1_{}_ncol{}{}".format(
        filename, ncol, fileformat)), dpi='figure', bbox_inches=bbox)
    plt.close(fig=fig)


def main(save_folder):
    adata = sc.read(os.path.join(
        '/Volumes', 'CH__data', 'Projects', 'Eyerich_AG_projects', 'BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer',
        'analysis', 'Molecular_subtypes', 'input', 'h5_files', 'NON_LESION',
        'LNL_RNAseq_20210720_patient_meta_data_v04__CH__Endotypes_230620.h5'))

    # Obervables to plot
    for obs, ncol in (['diag', 4], ['Pattern', 1]):
        if obs == 'diag':
            adata, lnl_order = add_colors.diag_order_lesion_nonlesion(adata=adata)
        else:
            adata, lnl_order = add_colors.pattern_order_lesion_nonlesion(adata=adata)
        adata.obs[obs] = adata.obs[obs].cat.reorder_categories(lnl_order)
        mapp_diag_color = dict(zip(list(adata.obs[obs].cat.categories), adata.uns['{}_colors'.format(obs)]))
        colors = [mapp_diag_color[val] for val in adata.obs[obs]]

        df_coordinates = pd.DataFrame.from_dict(
            {'x': adata.obsm['X_umap'].T[0], 'y': adata.obsm['X_umap'].T[1],
             obs: adata.obs[obs], 'color': colors})
        df_coordinates.to_excel(os.path.join(save_folder, 'Figure_S1_UMAP_{}.xlsx'.format(obs)))

        fig, ax = plt.subplots(figsize=(6, 6))
        sc.pl.umap(adata, color=obs, palette=list(adata.uns['{}_colors'.format(obs)]), use_raw=False, size=160,
                   show=False, ax=ax, title='', legend_loc='center right')

        # Put a legend below current axis
        # ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05), ncol=3, frameon=False, fontsize=16, markerscale=3.)

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        fig.tight_layout()
        plt.savefig(os.path.join(save_folder, "Figure_S1_UMAP_{}.pdf".format(obs)))
        plt.close(fig=fig)

        colordict = dict(zip(adata.obs[obs].cat.categories, adata.uns['{}_colors'.format(obs)]))
        f = lambda m, c: plt.plot([], [], marker='o', color=c, ls="none")[0]
        handles = [f("s", adata.uns['{}_colors'.format(obs)][i]) for i in range(len(adata.uns['{}_colors'.format(obs)]))]
        labels = list(colordict.keys())
        legend = plt.legend(handles, labels, loc=3, framealpha=1, frameon=False, ncol=ncol)
        export_legend(legend, save_folder=save_folder, ncol=ncol, filename="legend_{}".format(obs))


if __name__ == '__main__':
    output_dir = os.path.join('/Volumes', 'CH__data', 'Projects', 'Eyerich_AG_projects',
                              'BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer', 'analysis', 'Molecular_subtypes',
                              'output', 'Figure_S1_UMAP_diag_Pattern', str(date.today()))
    os.makedirs(output_dir, exist_ok=True)

    main(save_folder=output_dir)
