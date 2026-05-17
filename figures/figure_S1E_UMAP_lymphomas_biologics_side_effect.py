"""Plot diseases in UMAP

"""
from analysis.utils import add_colors

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
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
    fig.savefig(os.path.join(save_folder, "Figure_S1E_{}_ncol{}{}".format(
        filename, ncol, fileformat)), dpi='figure', bbox_inches=bbox)
    plt.close(fig=fig)


def main(save_folder):
    adata = sc.read(os.path.join(
        '/Volumes', 'CH__data', 'Projects', 'Eyerich_AG_projects', 'BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer',
        'analysis', 'Molecular_subtypes', 'input', 'h5_files', 'LESION',
        'Lesion_RNAseq_20210720_patient_meta_data_v04__CH_20220209_TherapyResponse_v04_v220715_TherapyResponseCorrected__Endotypes_230620.h5'))

    # Two-layer UMAP:
    # - all cells in grey background
    # - lymphoma + biologics side-effect cells overlaid in color
    # - colored by sdiag (subtype/clinical phenotype), not diag
    # - diagnoses to highlight
    highlight_diag = [
        'cutaneous lymphoma',
        'cutaneous side effects of biologics'
    ]

    adata.obs['sdiag'] = adata.obs['sdiag'].astype('object')
    adata.obs['sdiag'] = adata.obs['sdiag'].replace('nan', np.nan).fillna('N/A')
    adata.obs['sdiag'] = adata.obs['sdiag'].astype('category')

    # old categories + colors
    cats = list(adata.obs['sdiag'].cat.categories)
    colors = list(adata.uns['sdiag_colors'])

    # try to find where NaN used to be (often last or missing)
    # if NaN existed as a category before replacement, it is usually dropped now,
    # so we assume it was the last color or explicitly stored
    na_color = colors[-1]  # safest fallback assumption
    # rebuild categories ensuring N/A is included
    new_cats = cats
    # ensure N/A is present
    if 'N/A' not in new_cats:
        new_cats.append('N/A')
    # assign color mapping
    color_map = dict(zip(new_cats, [na_color] + colors[:len(new_cats) - 1]))
    adata.uns['sdiag_colors'] = [
        color_map[c] for c in adata.obs['sdiag'].cat.categories
    ]

    # mask for highlighted cells
    mask = adata.obs['diag'].isin(highlight_diag)

    # UMAP coordinates
    umap = adata.obsm['X_umap']

    # Scanpy categorical color mapping
    sdiag_palette = dict(
        zip(
            adata.obs['sdiag'].cat.categories,
            adata.uns['sdiag_colors']
        )
    )

    fig, ax = plt.subplots(figsize=(9, 6))
    # --- background layer ---
    ax.scatter(
        umap[:, 0],
        umap[:, 1],
        c='gainsboro',
        s=40,
        alpha=1,
        linewidths=0
    )
    # --- overlay highlighted groups ---
    for sdiag in adata.obs.loc[mask, 'sdiag'].unique():
        submask = mask & (adata.obs['sdiag'] == sdiag)

        ax.scatter(
            umap[submask, 0],
            umap[submask, 1],
            s=40,
            color=sdiag_palette[sdiag],
            label=sdiag,
            linewidths=0
        )
    ax.set_xlabel('UMAP1', fontsize=12)
    ax.set_ylabel('UMAP2', fontsize=12)
    ax.set_xticks([])
    ax.set_yticks([])
    sns.despine(ax=ax)

    grey_point = Line2D(
        [0], [0],
        marker='o',
        color='w',
        label='Other diagnosis',
        markerfacecolor='gainsboro',
        markersize=8
    )
    handles, labels = ax.get_legend_handles_labels()

    ax.legend(
        handles=[grey_point] + handles,
        labels=['Other diagnosis'] + labels,
        bbox_to_anchor=(1.05, 1),
        loc='upper left',
        frameon=False,
        title='Subtype/clinical phenotype'
    )
    plt.tight_layout()
    plt.savefig(os.path.join(save_folder, "Figure_S1E_UMAP_sdiag.pdf"))
    plt.close(fig=fig)

    # Build dataframe for reproducibility
    df_coordinates = pd.DataFrame(
        adata.obsm['X_umap'],
        columns=['UMAP1', 'UMAP2'],
        index=adata.obs_names
    )
    df_coordinates['diag'] = adata.obs['diag'].values
    df_coordinates['sdiag'] = adata.obs['sdiag'].values

    highlight_diag = [
        'cutaneous lymphoma',
        'cutaneous side effects of biologics'
    ]
    df_coordinates['highlight'] = df_coordinates['diag'].isin(highlight_diag)

    df_coordinates['plot_group'] = df_coordinates['sdiag'].copy()
    df_coordinates['plot_group'] = df_coordinates['plot_group'].astype(object)
    df_coordinates.loc[~df_coordinates['highlight'], 'plot_group'] = 'Other diagnosis'

    df_coordinates.to_excel(os.path.join(save_folder, 'Figure_S1E_UMAP_sdiag.xlsx'))


if __name__ == '__main__':
    output_dir = os.path.join('/Volumes', 'CH__data', 'Projects', 'Eyerich_AG_projects',
                              'BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer', 'analysis', 'Molecular_subtypes',
                              'output', 'Figure_S1E_UMAP_subdiagnosis', str(date.today()))
    os.makedirs(output_dir, exist_ok=True)

    main(save_folder=output_dir)
