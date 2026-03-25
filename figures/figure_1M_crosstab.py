from scripts.utils import add_colors
from scripts.GEx_classifier.utils import utils

import matplotlib.pyplot as plt
import matplotlib
import scanpy as sc
import os
from datetime import date
# %matplotlib notebook

import pandas as pd

import seaborn as sns
# from pySankey.sankey import sankey


def main(save_folder):
    adata = sc.read(
        os.path.join(
            '/Volumes', 'CH__data', 'Projects', 'Eyerich_AG_projects', 'BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer',
            'analysis', 'Molecular_subtypes', 'input', 'h5_files', 'LESION',
            'Lesion_RNAseq_20210720_patient_meta_data_v04__CH_20220209_TherapyResponse_v04_v220715_TherapyResponseCorrected__Endotypes_230620.h5'))

    # adata, _ = add_colors.diag_order_lesion(adata=adata)
    # adata, _ = add_colors.sdiag_order_lesion(adata=adata)
    # adata, _ = add_colors.pattern_order_lesion(adata=adata)

    df = pd.crosstab(adata.obs['Pattern'], adata.obs['diag'], normalize=False)
    my_colors = []
    for c in adata.uns['diag_colors']:
        my_colors.append(matplotlib.colors.hex2color(c))

    figsize = (9.5, 3)

    fig, ax = plt.subplots(figsize=figsize)
    p = df.plot.bar(
        stacked=True, color=my_colors, ax=ax)
    p.get_legend().remove()
    ax.set_xticklabels(['1', '2a', '2b', '3', '4a', '4b', '5', 'UD'], rotation=0)
    # ax.legend(loc='center left', bbox_to_anchor=(1, 0.5),
    #           fancybox=False, shadow=False, ncol=2, prop={'size': 12})
    ax.set_ylabel('samples', fontsize=24)
    ax.set_xlabel('pattern', fontsize=24)
    ax.xaxis.set_tick_params(labelsize=20)
    ax.yaxis.set_tick_params(labelsize=20)
    sns.despine(fig=fig, ax=ax)
    plt.tight_layout()
    plt.savefig(os.path.join(save_folder, 'Figure_1M_Crosstab_pattern_diag.pdf'))
    plt.close('all')

    # Save to dataframe
    df.loc[len(df)] = my_colors
    df.rename(index={8: "colors"})
    df.to_excel(os.path.join(save_folder, 'Figure_1M_crosstab_Pattern_diag.xlsx'))

    # Crosstab sub-diseases
    df = pd.crosstab(adata.obs['Pattern'], adata.obs['sdiag'], normalize=False)
    my_colors = []
    for c in adata.uns['sdiag_colors']:
        my_colors.append(matplotlib.colors.hex2color(c))

    fig, ax = plt.subplots(figsize=figsize)  # (w, h)
    df.plot.bar(stacked=True, color=my_colors, ax=ax)
    ax.get_legend().remove()
    # ax.legend(loc='center left', bbox_to_anchor=(1, 0.5),
    #           fancybox=False, shadow=False, ncol=2, prop={'size': 12})
    ax.set_xticklabels(['1', '2a', '2b', '3', '4a', '4b', '5', 'UD'], rotation=0)
    ax.set_ylabel('counts', fontsize=24)
    ax.set_xlabel('pattern', fontsize=24)
    ax.xaxis.set_tick_params(labelsize=20)
    ax.yaxis.set_tick_params(labelsize=20)
    sns.despine(fig=fig, ax=ax)
    plt.tight_layout()
    plt.savefig(os.path.join(save_folder, 'Crosstab_pattern_subdiseases.pdf'))
    plt.close('all')


if __name__ == '__main__':
    output_dir = os.path.join('/Volumes', 'CH__data', 'Projects', 'Eyerich_AG_projects',
                              'BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer', 'analysis', 'Molecular_subtypes',
                              'output', 'Figure_1M_crosstab', str(date.today()))
    os.makedirs(output_dir, exist_ok=True)

    main(save_folder=output_dir)
