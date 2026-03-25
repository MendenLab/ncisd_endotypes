from scripts.utils import add_colors
import matplotlib.pyplot as plt
import scanpy as sc
import numpy as np

import os
from datetime import date
# %matplotlib notebook


def diag_order_lesion(adata):
    diag_order = [
        'lichen planus', 'lupus erythematosus', 'eczema', 'prurigo simplex subacuta', 'bullous pemphigoid', 'psoriasis',
        'pityriasis rubra pilaris', 'morphea', 'venous ulcer', 'systemic sclerosis', 'granuloma annulare',
        'sarcoidosis', 'psoriasis pustulosa', 'pyoderma gangrenosum', 'cutaneous lymphoma',
        'cutaneous side effects of biologics', 'darier disease', 'keratosis lichenoides chronica', 'erythrodermia',
        'parapsoriasis', 'undefined']
    adata.uns['diag_colors'] = [
        'orange', 'goldenrod', 'maroon', 'firebrick', 'salmon', 'blue', 'dodgerblue', 'forestgreen', 'darkolivegreen',
        'seagreen', 'springgreen', 'palegreen', 'darkviolet', 'mediumorchid', 'dimgrey', 'darkgray', 'silver', 'grey',
        'slategrey', 'lightgrey', 'whitesmoke']

    return adata, diag_order


def diag_order_lesion_nonlesion(adata):
    diag_order = [
        'lichen planus', 'lupus erythematosus', 'eczema', 'prurigo simplex subacuta', 'bullous pemphigoid', 'psoriasis',
        'pityriasis rubra pilaris', 'morphea', 'venous ulcer', 'systemic sclerosis', 'granuloma annulare',
        'sarcoidosis', 'psoriasis pustulosa', 'pyoderma gangrenosum', 'cutaneous lymphoma',
        'cutaneous side effects of biologics', 'darier disease', 'keratosis lichenoides chronica', 'erythrodermia',
        'parapsoriasis', 'undefined', 'non lesional']
    adata.uns['diag_colors'] = [
        'orange', 'goldenrod', 'maroon', 'firebrick', 'salmon', 'blue', 'dodgerblue', 'forestgreen', 'darkolivegreen',
        'seagreen', 'springgreen', 'palegreen', 'darkviolet', 'mediumorchid', 'dimgrey', 'darkgray', 'silver', 'grey',
        'slategrey', 'lightgrey', 'whitesmoke', 'cornsilk']

    return adata, diag_order


def get_outer_ring(adata, diag_order):
    group_outer = adata.obs['diag'].value_counts(dropna=False)
    group_outer = group_outer.reindex(diag_order)
    group_outer_name = list(group_outer.index)

    return group_outer, group_outer_name


# def get_middle_ring(adata):
#     group_middle = adata.obs['Endotypes'].value_counts(dropna=False)
#
#     from operator import itemgetter
#     map_endotypes_pattern = list(zip(adata.obs['Endotypes'], adata.obs['Pattern']))
#     sorted(map_endotypes_pattern, key=itemgetter(1))
#
#     group_middle = group_middle.reindex(adata.obs['Endotypes'].cat.categories)
#
#     group_middle_name = list(group_middle.index)
#
#     # sort by diag or disease
#
#     return group_middle, group_middle_name


def get_inner_ring(adata):
    group_inner = adata.obs['Pattern'].value_counts(dropna=False)
    group_inner = group_inner.reindex(adata.obs['Pattern'].cat.categories)
    group_inner_name = list(group_inner.index)

    return group_inner, group_inner_name


def plot_nested_3_layer_piechart(adata, group_outer, group_middle, group_inner,
                                 group_outer_name, group_middle_name, group_inner_name, save_folder, key):

    radius = 2
    width = 0.3

    explode = (0, ) * len(group_outer.values)

    # First Ring (outside)
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.axis('equal')
    # First Ring (outside)
    mypie, _ = ax.pie(group_outer.values, radius=radius,
                      colors=adata.uns['diag_colors'], counterclock=False, startangle=90,
                      explode=explode, wedgeprops={'linewidth': 0.5, 'edgecolor': 'white'})
    plt.setp(mypie, width=width, edgecolor='white')

    # Second Ring (middle)
    mypie2, _ = ax.pie(group_middle.values, radius=radius - width,
                       labels=group_middle_name, labeldistance=0.85,
                       colors=adata.uns['Endotypes_colors'], counterclock=False, startangle=90,
                       wedgeprops={'linewidth': 0.5, 'edgecolor': 'white'})
    plt.setp(mypie2, width=width, edgecolor='white')

    # Third Ring (inside)
    mypie3, _ = ax.pie(group_inner.values, radius=radius - (2 * width),
                       labels=group_inner_name, labeldistance=0.5,
                       colors=adata.uns['Pattern_colors'], counterclock=False, startangle=90,
                       wedgeprops={'linewidth': 0.5, 'edgecolor': 'white'})
    plt.setp(mypie3, width=width, edgecolor='white')

    plt.margins(0, 0)
    ax.legend(mypie, group_outer_name, loc='upper center', bbox_to_anchor=(0.5, -0.4), frameon=False,
              fancybox=False, shadow=False, ncol=4, prop={'size': 6})
    plt.tight_layout()

    plt.savefig(os.path.join(save_folder, 'Nested_3_layer_Piechart_{}.pdf'.format(key)))
    plt.close(fig=fig)


def plot_nested_piechart(adata, group_outer, group_inner, group_outer_name, group_inner_name, save_folder, key):
    explode = (0, ) * len(group_outer.values)
    # First Ring (outside)
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.axis('equal')
    mypie, _ = ax.pie(group_outer.values, radius=1.3,
                      colors=adata.uns['diag_colors'], counterclock=False, startangle=90,
                      explode=explode, wedgeprops={'linewidth': 0.5, 'edgecolor': 'white'})
    plt.setp(mypie, width=0.3, edgecolor='white')
    # Second Ring (inside)
    mypie2, _ = ax.pie(group_inner.values, radius=1.3 - 0.3,
                       labels=group_inner_name, labeldistance=0.7,
                       colors=adata.uns['Pattern_colors'], counterclock=False, startangle=90,
                       wedgeprops={'linewidth': 0.5, 'edgecolor': 'white'})
    plt.setp(mypie2, width=0.4, edgecolor='white')
    plt.margins(0, 0)
    ax.legend(mypie, group_outer_name, loc='center right', bbox_to_anchor=(-0.05, 0.5),
              fancybox=False, shadow=False, ncol=1, prop={'size': 6})
    plt.tight_layout()

    plt.savefig(os.path.join(save_folder, 'Nested_Piechart_{}.pdf'.format(key)))
    plt.close(fig=fig)


def plot_nested_piechart_reversed(adata, group_outer, group_inner, group_outer_name, save_folder, key):
    explode = (0, ) * len(group_outer.values)
    # First Ring (outside)
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.axis('equal')
    mypie, _ = ax.pie(group_outer.values, radius=1.3,
                      colors=adata.uns['Pattern_colors'], counterclock=False, startangle=90,
                      explode=explode, wedgeprops={'linewidth': 0.5, 'edgecolor': 'white'})
    plt.setp(mypie, width=0.3, edgecolor='white')
    # Second Ring (inside)
    mypie2, _ = ax.pie(group_inner.values, radius=1.3 - 0.3,
                       labeldistance=0.7, colors=adata.uns['diag_colors'], counterclock=False, startangle=90,
                       wedgeprops={'linewidth': 0.5, 'edgecolor': 'white'})
    plt.setp(mypie2, width=0.4, edgecolor='white')
    plt.margins(0, 0)

    bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=0.72)
    kw = dict(arrowprops=dict(arrowstyle="-"), bbox=bbox_props, zorder=0, va="center",  fontsize=10)

    for i, p in enumerate(mypie):
        ang = (p.theta2 - p.theta1) / 2. + p.theta1
        y = np.sin(np.deg2rad(ang))
        x = np.cos(np.deg2rad(ang))
        horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
        connectionstyle = "angle,angleA=0,angleB={}".format(ang)
        kw["arrowprops"].update({"connectionstyle": connectionstyle})

        # To make every second arrow longer
        ax.annotate(group_outer_name[i], xy=(x, y), xytext=((1.4 + (i % 2) * 0.4) * np.sign(x), 1.4 + y),
                    horizontalalignment=horizontalalignment, **kw)
    # plt.tight_layout()

    plt.savefig(os.path.join(save_folder, 'Nested_Piechart_reversed_{}.pdf'.format(key)), bbox_inches='tight')
    plt.close(fig=fig)


def plot_nested_piechart_reversed_lesion(adata, group_outer, group_inner, group_outer_name, save_folder, key):
    explode = (0, ) * len(group_outer.values)
    # First Ring (outside)
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.axis('equal')
    mypie, _ = ax.pie(group_outer.values, radius=1.3,
                      colors=adata.uns['Pattern_colors'], counterclock=False, startangle=90,
                      explode=explode, wedgeprops={'linewidth': 0.5, 'edgecolor': 'white'})
    plt.setp(mypie, width=0.3, edgecolor='white')
    # Second Ring (inside)
    mypie2, _ = ax.pie(group_inner.values, radius=1.3 - 0.3,
                       labeldistance=0.7, colors=adata.uns['diag_colors'], counterclock=False, startangle=90,
                       wedgeprops={'linewidth': 0.5, 'edgecolor': 'white'})
    plt.setp(mypie2, width=0.4, edgecolor='white')
    plt.margins(0, 0)

    bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=0.72)
    kw = dict(arrowprops=dict(arrowstyle="-"), bbox=bbox_props, zorder=0, va="center",  fontsize=18)

    for i, p in enumerate(mypie):
        ang = (p.theta2 - p.theta1) / 2. + p.theta1
        y = np.sin(np.deg2rad(ang))
        x = np.cos(np.deg2rad(ang))
        horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
        connectionstyle = "angle,angleA=0,angleB={}".format(ang)
        kw["arrowprops"].update({"connectionstyle": connectionstyle})

        # To make every second arrow shorter
        ax.annotate(group_outer_name[i], xy=(x, y), xytext=((1.4 + (i % 2) * 0.8) * np.sign(x), 1.4 * y),
                    horizontalalignment=horizontalalignment, **kw)

    plt.savefig(os.path.join(save_folder, 'Nested_Piechart_reversed_{}.pdf'.format(key)), bbox_inches='tight')
    plt.close(fig=fig)


def plot_nested_3_layer_piechart_reversed_lesion(
        adata, group_outer, group_middle, group_inner, group_outer_name, group_middle_name, group_inner_name,
        save_folder, key):
    radius = 1.5
    width = 0.3

    explode = (0, ) * len(group_outer.values)

    # First Ring (outside)
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.axis('equal')
    # First Ring (outside)
    mypie, _ = ax.pie(group_outer.values, radius=radius,
                      colors=adata.uns['Pattern_colors'], counterclock=False, startangle=90,
                      explode=explode, wedgeprops={'linewidth': 0.5, 'edgecolor': 'white'})
    plt.setp(mypie, width=width, edgecolor='white')

    # Second Ring (middle)
    mypie2, _ = ax.pie(group_middle.values, radius=radius - width,
                       labels=group_middle_name, labeldistance=0.85,
                       colors=adata.uns['Endotypes_colors'], counterclock=False, startangle=90,
                       wedgeprops={'linewidth': 0.5, 'edgecolor': 'white'})
    plt.setp(mypie2, width=width, edgecolor='white')

    # Third Ring (inside)
    mypie3, _ = ax.pie(group_inner.values, radius=radius - (2 * width), labeldistance=0.7,
                       colors=adata.uns['diag_colors'], counterclock=False, startangle=90,
                       wedgeprops={'linewidth': 0.5, 'edgecolor': 'white'})
    plt.setp(mypie3, width=width, edgecolor='white')

    bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=0.72)
    kw = dict(arrowprops=dict(arrowstyle="-"), bbox=bbox_props, zorder=0, va="center",  fontsize=18)

    for i, p in enumerate(mypie):
        ang = (p.theta2 - p.theta1) / 2. + p.theta1
        y = np.sin(np.deg2rad(ang))
        x = np.cos(np.deg2rad(ang))
        horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
        connectionstyle = "angle,angleA=0,angleB={}".format(ang)
        kw["arrowprops"].update({"connectionstyle": connectionstyle})

        # To make every second arrow shorter
        ax.annotate(group_outer_name[i], xy=(x, y), xytext=((1.4 + (i % 2) * 0.8) * np.sign(x), 1.4 * y),
                    horizontalalignment=horizontalalignment, **kw)

    plt.savefig(os.path.join(save_folder,
                             'Nested_3_layer_Piechart_reversed_{}.pdf'.format(key)), bbox_inches='tight')
    plt.close(fig=fig)


def main(save_folder):
    adata = sc.read(os.path.join(
        '/Volumes', 'CH__data', 'Projects', 'Eyerich_AG_projects', 'BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer',
        'analysis', 'Molecular_subtypes', 'input', 'h5_files', 'LESION',
        'Lesion_RNAseq_20210720_patient_meta_data_v04__CH_20220209_TherapyResponse_v04_v220715_TherapyResponseCorrected__Endotypes_230620.h5'))

    file_name = '20210720_patient_meta_data_v04__CH__Endotypes_230620'
    h5file_name = os.path.join('/Volumes', 'CH__data', 'Projects', 'Eyerich_AG_projects',
                               'BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer', 'analysis', 'Molecular_subtypes',
                               'input', 'h5_files', 'NON_LESION', 'LNL_RNAseq_{}.h5'.format(file_name))
    adata_lnl = sc.read(h5file_name)

    # number of patients in dataset after filtering for bad quality samples
    print("Number of patients in L and NL dataset after filtering: ", len(adata_lnl.obs['Pseudo ID'].unique()))
    print("Number of samples in L and NL dataset after filtering: ", len(adata_lnl.obs.index.tolist()))

    # Values for piechart
    # reorder diagnosis and pattern
    adata, diag_order = add_colors.diag_order_lesion(adata=adata)
    adata_lnl, lnl_diag_order = add_colors.diag_order_lesion_nonlesion(adata=adata_lnl)
    adata, _ = add_colors.pattern_order_lesion(adata=adata)
    adata_lnl, _ = add_colors.pattern_order_lesion_nonlesion(adata=adata_lnl)

    # Get outer ring - diseases
    group_outer, group_outer_name = get_outer_ring(adata=adata, diag_order=diag_order)
    group_outer_lnl, group_outer_name_lnl = get_outer_ring(adata=adata_lnl, diag_order=lnl_diag_order)

    # Middle ring - Endotypes - TODO how to show that endotypes are randomly distributed?
    # group_middle, group_middle_name = get_middle_ring(adata=adata)  # only lesion samples

    # Inner ring -  Pattern
    group_inner, group_inner_name = get_inner_ring(adata=adata)
    # rename pattern 1-5 to Pattern 1 to Pattern 5
    group_inner_name = ['Pattern ' + s if s != 'UD' else s for s in group_inner_name]
    group_inner_lnl, group_inner_name_lnl = get_inner_ring(adata=adata_lnl)
    # rename pattern 1-5 to Pattern 1 to Pattern 5
    group_inner_name_lnl = ['Pattern ' + s if ((s != 'UD') and (s != 'NL')) else s for s in group_inner_name_lnl]

    # One Ring (Diseases + non-lesional)
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.axis('equal')
    mypie, _ = ax.pie(group_outer_lnl.values, radius=1.3,
                      colors=adata_lnl.uns['diag_colors'], counterclock=False, startangle=90,
                      wedgeprops={'linewidth': 1, 'edgecolor': 'white'})
    plt.setp(mypie, width=0.3, edgecolor='white')
    ax.legend(mypie, group_outer_name_lnl, loc='upper center', bbox_to_anchor=(0.5, -0.05),
              fancybox=False, shadow=False, ncol=3, prop={'size': 10}, frameon=False)
    plt.tight_layout()
    plt.savefig(os.path.join(save_folder, 'Disease_Piechart.pdf'))
    plt.close(fig=fig)

    # plot_nested_3_layer_piechart(
    #     adata=adata, group_outer=group_outer, group_middle=group_middle, group_inner=group_inner,
    #     group_outer_name=group_outer_name, group_middle_name=group_middle_name, group_inner_name=group_inner_name,
    #     save_folder=save_folder, key="L")

    plot_nested_piechart(
        adata=adata, group_outer=group_outer, group_inner=group_inner,
        group_outer_name=group_outer_name, group_inner_name=group_inner_name, save_folder=save_folder,
        key="L")

    plot_nested_piechart(
        adata=adata_lnl, group_outer=group_outer_lnl, group_inner=group_inner_lnl,
        group_outer_name=group_outer_name_lnl, group_inner_name=group_inner_name_lnl, save_folder=save_folder,
        key="L_NL")

    plot_nested_piechart_reversed_lesion(
        adata=adata, group_outer=group_inner, group_inner=group_outer,
        group_outer_name=group_inner_name, save_folder=save_folder,
        key="L")

    plot_nested_piechart_reversed(
        adata=adata_lnl, group_outer=group_inner_lnl, group_inner=group_outer_lnl,
        group_outer_name=group_inner_name_lnl, save_folder=save_folder,
        key="L_NL")

    # plot_nested_3_layer_piechart_reversed_lesion(
    #     adata=adata, group_outer=group_inner, group_middle=group_middle, group_inner=group_outer,
    #     group_outer_name=group_inner_name, group_middle_name=group_middle_name, group_inner_name=group_outer_name,
    #     save_folder=save_folder, key="L")


if __name__ == '__main__':
    output_dir = os.path.join('/Volumes', 'CH__data', 'Projects', 'Eyerich_AG_projects',
                              'BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer', 'analysis', 'Molecular_subtypes',
                              'output', 'Figure_1A', str(date.today()))
    os.makedirs(output_dir, exist_ok=True)

    main(save_folder=output_dir)
