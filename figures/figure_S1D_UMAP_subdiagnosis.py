"""Plot diseases in UMAP

"""
from analysis.utils import add_colors

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
    fig.savefig(os.path.join(save_folder, "Figure_S1D_{}_ncol{}{}".format(
        filename, ncol, fileformat)), dpi='figure', bbox_inches=bbox)
    plt.close(fig=fig)


def sdiag_order_lesion(adata):
    # 'eczemaundefined', 'nummular eczemaundefined', 'lupus erythematosusundefined',
    # 'psoriasis pustulosa palmoplantarisundefined', 'psoriasis pustulosaundefined'
    # sdiag_order = [
    #     'chilblain lupus', 'chronic discoid lupus erythematosus', 'subacute cutaneous lupus erythematosus',
    #     'lupus erythematosus',
    #     'erythrodermia', 'asteatotic eczema', 'atopic dermatitis', 'hyperkeratotic rhagadiform eczema of the hands',
    #     'nummular eczema', 'rosacea', 'seborrheic eczema', 'eczema',
    #     'plaque psoriasis',
    #     'psoriasis guttata', 'psoriasis inversa', 'psoriasis palmoplantaris', 'psoriasis pustulosa',
    #     'psoriasis pustulosa palmoplantaris', 'N/A']
    #
    # adata.obs['sdiag'] = adata.obs['sdiag'].astype('category')
    # adata.obs['sdiag'] = adata.obs['sdiag'].cat.reorder_categories(sdiag_order)
    # adata.uns['sdiag_colors'] = [
    #     'sandybrown', 'peru', 'bisque', 'peachpuff',
    #     'lightcoral', 'indianred', 'orangered', 'brown', 'salmon', 'tomato', 'red', 'maroon',
    #     'midnightblue',
    #     'deepskyblue', 'blue', 'dodgerblue', 'steelblue',
    #     'darkviolet', 'grey']

    lupus_colors = [
        "#e76f51",  # strong orange-red (darkest)
        "#f4a261",  # soft orange
        "#ffd166",  # yellow-orange
        "#fff3b0"  # light yellow (lightest)
    ]

    eczema_colors = [
        "#660000",  # darkest
        "#7f0000",
        "#99000d",
        "#a50f15",
        "#cb181d",
        "#e41a1c",
        "#f08080",
        "#fbb4ae"  # lightest
    ]

    psoriasis_colors = [
        "#08306b",  # darkest navy
        "#08519c",
        "#3182bd",
        "#6baed6",
        "#9ecae1",
        "#deebf7"  # lightest
    ]

    na_color = ["#bdbdbd"]

    sdiag_order = [
        'chilblain lupus',
        'chronic discoid lupus erythematosus',
        'subacute cutaneous lupus erythematosus',
        'lupus erythematosus',

        'erythrodermia',
        'asteatotic eczema',
        'atopic dermatitis',
        'hyperkeratotic rhagadiform eczema of the hands',
        'nummular eczema',
        'rosacea',
        'seborrheic eczema',
        'eczema',

        'plaque psoriasis',
        'psoriasis guttata',
        'psoriasis inversa',
        'psoriasis palmoplantaris',
        'psoriasis pustulosa',
        'psoriasis pustulosa palmoplantaris',

        'N/A'
    ]

    adata.obs['sdiag'] = adata.obs['sdiag'].astype('category')
    # print(set(adata.obs["sdiag"].cat.categories) - set(sdiag_order))
    # print(set(sdiag_order) - set(adata.obs["sdiag"].cat.categories))
    adata.obs['sdiag'] = adata.obs['sdiag'].cat.reorder_categories(sdiag_order)

    adata.uns['sdiag_colors'] = (
            lupus_colors +
            eczema_colors +
            psoriasis_colors +
            na_color
    )

    return adata, sdiag_order


def main(save_folder):
    adata = sc.read(os.path.join(
        '/Volumes', 'CH__data', 'Projects', 'Eyerich_AG_projects', 'BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer',
        'analysis', 'Molecular_subtypes', 'input', 'h5_files', 'LESION',
        'Lesion_RNAseq_20210720_patient_meta_data_v04__CH_20220209_TherapyResponse_v04_v220715_TherapyResponseCorrected__Endotypes_230620.h5'))

    # Obervable to plot
    obs = 'sdiag'
    ncol = 2

    # Put into one group (only plaques were biopsied):
    # Rename to "plaque psoriasis and psoriasis arthritis" and "plaques psoriasis and psoriasis inversa" to "plaque psoriasis"
    adata.obs[obs] = adata.obs[obs].replace({
        "plaque psoriasis and psoriasis arthritis": "plaque psoriasis",
        "plaque psoriasis and psoriasis inversa": "plaque psoriasis",
        # Generalized pustular psoriasis and psoriasis pustulosa are the same
        "generalized pustular psoriasis": "psoriasis pustulosa",
        # Replace nan with N/A
        "nan": "N/A"
    })
    adata.obs[obs] = adata.obs[obs].cat.remove_unused_categories()

    # Read out colors for subtypes
    adata, _ = sdiag_order_lesion(adata=adata)
    mapp_diag_color = dict(zip(list(adata.obs[obs].cat.categories), adata.uns['{}_colors'.format(obs)]))
    colors = [mapp_diag_color[val] for val in adata.obs[obs]]

    df_coordinates = pd.DataFrame.from_dict(
        {'x': adata.obsm['X_umap'].T[0], 'y': adata.obsm['X_umap'].T[1],
         obs: adata.obs[obs], 'color': colors})
    df_coordinates.to_excel(os.path.join(save_folder, 'Figure_S1D_UMAP_{}.xlsx'.format(obs)))

    fig, ax = plt.subplots(figsize=(6, 6))
    sc.pl.umap(adata, color=obs, use_raw=False, size=160,
               show=False, ax=ax, title='', legend_loc='center right')

    # Put a legend below current axis
    # ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05), ncol=3, frameon=False, fontsize=16, markerscale=3.)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout()
    plt.savefig(os.path.join(save_folder, "Figure_S1D_UMAP_{}.pdf".format(obs)))
    plt.close(fig=fig)

    colordict = dict(zip(adata.obs[obs].cat.categories, adata.uns['{}_colors'.format(obs)]))
    f = lambda m, c: plt.plot([], [], marker='o', color=c, ls="none")[0]
    handles = [f("s", adata.uns['{}_colors'.format(obs)][i]) for i in range(len(adata.uns['{}_colors'.format(obs)]))]
    labels = list(colordict.keys())
    legend = plt.legend(handles, labels, loc=3, framealpha=1, frameon=False, ncol=ncol)
    export_legend(legend, save_folder=save_folder, ncol=ncol)

    ct = pd.crosstab(adata.obs["sdiag"], adata.obs["Pattern"])
    row_colors = pd.Series(
        adata.uns["sdiag_colors"],
        index=adata.obs["sdiag"].cat.categories
    ).loc[ct.index]

    fig, ax = plt.subplots(figsize=(10, 6))
    # draw colored background
    ax.imshow(np.zeros_like(ct), cmap="Greys", aspect="auto")
    # apply row colors manually
    for i, row in enumerate(ct.index):
        ax.add_patch(plt.Rectangle(
            (-0.5, i - 0.5),
            ct.shape[1],
            1,
            color=row_colors[row],
            alpha=0.35
        ))
    # add numbers
    for i in range(ct.shape[0]):
        for j in range(ct.shape[1]):
            ax.text(
                j, i,
                str(ct.iloc[i, j]),
                ha="center",
                va="center",
                fontsize=12
            )
    # ticks
    ax.set_xticks(range(ct.shape[1]))
    ax.set_xticklabels(ct.columns, rotation=0, fontsize=12)

    ax.set_yticks(range(ct.shape[0]))
    ax.set_yticklabels(ct.index, fontsize=12)

    ax.set_title("Subtypes vs. Pattern")
    plt.tight_layout()
    plt.show()
    plt.savefig(os.path.join(save_folder, "Figure_S1D_Table_{}_vs_Pattern.pdf".format(obs)))
    plt.close(fig=fig)


if __name__ == '__main__':
    output_dir = os.path.join('/Volumes', 'CH__data', 'Projects', 'Eyerich_AG_projects',
                              'BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer', 'analysis', 'Molecular_subtypes',
                              'output', 'Figure_S1D_UMAP_subdiagnosis', str(date.today()))
    os.makedirs(output_dir, exist_ok=True)

    main(save_folder=output_dir)
