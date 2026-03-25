import numpy as np
import scanpy as sc
import re


def diag_order_lesion_nonlesion(adata):
    diag_order = [
        'lichen planus', 'lupus erythematosus', 'lichenoid drug reaction', 'eczema', 'prurigo simplex subacuta',
        'bullous pemphigoid', 'psoriasis', 'pityriasis rubra pilaris', 'morphea', 'venous ulcer', 'systemic sclerosis',
        'granuloma annulare', 'sarcoidosis', 'psoriasis pustulosa', 'pyoderma gangrenosum', 'cutaneous lymphoma',
        'cutaneous side effects of biologics', 'darier disease', 'keratosis lichenoides chronica', 'erythrodermia',
        'parapsoriasis', 'non lesional']
    adata.uns['diag_colors'] = [
        'orange', 'goldenrod', 'gold', 'maroon', 'firebrick', 'salmon', 'blue',
        'dodgerblue', 'palevioletred', 'mediumvioletred', 'deeppink', 'fuchsia',
        'violet', 'darkviolet', 'mediumorchid', 'dimgrey',
        'darkgray', 'silver', 'grey', 'slategrey',
        'lightgrey', 'pink']  # 'cornsilk' 'aquamarine'

    return adata, diag_order


def sdiag_order_lesion_nonlesion(adata):
    # 'eczemaundefined', 'nummular eczemaundefined', 'lupus erythematosusundefined',
    # 'psoriasis pustulosa palmoplantarisundefined', 'psoriasis pustulosaundefined'
    sdiag_order = [
        'chilblain lupus', 'chronic discoid lupus erythematosus', 'subacute cutaneous lupus erythematosus',
        'lupus erythematosus',
        'erythrodermia', 'asteatotic eczema', 'atopic dermatitis', 'hyperkeratotic rhagadiform eczema of the hands',
        'nummular eczema', 'rosacea', 'seborrheic eczema', 'eczema',
        'plaque psoriasis and psoriasis arthritis', 'plaque psoriasis and psoriasis inversa', 'plaque psoriasis',
        'psoriasis guttata', 'psoriasis inversa', 'psoriasis palmoplantaris', 'psoriasis pustulosa',
        'psoriasis pustulosa palmoplantaris', 'generalized pustular psoriasis',
        'nan']

    adata.obs['sdiag'] = adata.obs['sdiag'].astype(str).astype('category')
    adata.obs['sdiag'] = adata.obs['sdiag'].cat.reorder_categories(sdiag_order)
    adata.uns['sdiag_colors'] = [
        'sandybrown', 'peru', 'bisque', 'peachpuff',
        'lightcoral', 'indianred', 'orangered', 'brown', 'salmon', 'tomato', 'red', 'maroon',
        'cornflowerblue', 'royalblue', 'midnightblue', 'deepskyblue', 'blue', 'dodgerblue', 'steelblue',
        'darkviolet', 'mediumorchid',
        'grey']

    return adata, sdiag_order


def pattern_order_lesion_nonlesion(adata):
    # Replace undefined with UD and non-lesional with NL
    adata.obs = adata.obs.replace({'Pattern': {"undefined": "UD", "non-lesional": "NL"}})

    # Reorder categories
    pattern_order = ['1', '2a', '2b', '3', '4a', '4b', '5', 'UD', 'NL']
    adata.obs['Pattern'] = adata.obs['Pattern'].cat.reorder_categories(pattern_order)
    adata.uns['Pattern_colors'] = ['darkorange', 'darkred', 'tomato', 'royalblue', 'mediumvioletred',
                                   'magenta', 'darkviolet', 'grey', 'pink']  # 'aquamarine'

    return adata, pattern_order


def endotype_order_lesion_nonlesion(adata, obs_name):
    endotype_order = list(np.sort(adata.obs.loc[adata.obs[obs_name] != 'non lesional', obs_name].astype(
        str).map(lambda x: x.lstrip('E')).astype(int).unique()))
    endotype_order = np.asarray(endotype_order + ['non lesional'])
    len_obsname = len(endotype_order)

    if obs_name == 'Endotypes':
        endotype_order = [s if s == 'non lesional' else 'E' + s for s in list(endotype_order.astype(str))]
        adata.obs[obs_name] = adata.obs[obs_name].cat.reorder_categories(endotype_order)
    else:
        adata.obs[obs_name] = adata.obs[obs_name].cat.reorder_categories(endotype_order)

    try:
        adata.uns["{}_colors".format(obs_name)] = ["#006ddb", "#b6dbff", "#004949", "#009292", "#ff6db6", "#490092",
                                                   "#b66dff", "#000000", "#920000", "#E69F00", "#D55E00", "#8B4513",
                                                   "#999999", 'aquamarine'][:len_obsname]
    except IndexError:
        adata.uns["{}_colors".format(obs_name)] = sc.pl.palettes.default_102[:len_obsname]

    return adata, endotype_order


def diag_order_lesion(adata):
    diag_order = [
        'lichen planus', 'lupus erythematosus', 'lichenoid drug reaction', 'eczema', 'prurigo simplex subacuta',
        'bullous pemphigoid', 'psoriasis', 'pityriasis rubra pilaris', 'morphea', 'venous ulcer', 'systemic sclerosis',
        'granuloma annulare', 'sarcoidosis', 'psoriasis pustulosa', 'pyoderma gangrenosum', 'cutaneous lymphoma',
        'cutaneous side effects of biologics', 'darier disease', 'keratosis lichenoides chronica', 'erythrodermia',
        'parapsoriasis']
    adata.obs['diag'] = adata.obs['diag'].cat.reorder_categories(diag_order)
    adata.uns['diag_colors'] = [
        'orange', 'goldenrod', 'gold', 'maroon', 'firebrick', 'salmon', 'blue',
        'dodgerblue', 'palevioletred', 'mediumvioletred', 'deeppink', 'fuchsia',
        'violet', 'darkviolet', 'mediumorchid', 'dimgrey',
        'darkgray', 'silver', 'grey', 'slategrey',
        'lightgrey']

    return adata, diag_order


def sdiag_order_lesion(adata):
    # 'eczemaundefined', 'nummular eczemaundefined', 'lupus erythematosusundefined',
    # 'psoriasis pustulosa palmoplantarisundefined', 'psoriasis pustulosaundefined'
    sdiag_order = [
        'chilblain lupus', 'chronic discoid lupus erythematosus', 'subacute cutaneous lupus erythematosus',
        'lupus erythematosus',
        'erythrodermia', 'asteatotic eczema', 'atopic dermatitis', 'hyperkeratotic-rhagadiform eczema of the hands',
        'nummular eczema', 'rosacea', 'seborrheic eczema', 'eczema',
        'plaque psoriasis and psoriasis arthritis', 'plaque psoriasis and psoriasis inversa', 'plaque psoriasis',
        'psoriasis guttata', 'psoriasis inversa', 'psoriasis palmoplantaris', 'psoriasis pustulosa',
        'psoriasis pustulosa palmoplantaris', 'generalized pustular psoriasis']

    adata.obs['sdiag'] = adata.obs['sdiag'].astype('category')
    adata.obs['sdiag'] = adata.obs['sdiag'].cat.reorder_categories(sdiag_order)
    adata.uns['sdiag_colors'] = [
        'sandybrown', 'peru', 'bisque', 'peachpuff',
        'lightcoral', 'indianred', 'orangered', 'brown', 'salmon', 'tomato', 'red', 'maroon',
        'cornflowerblue', 'royalblue', 'midnightblue', 'deepskyblue', 'blue', 'dodgerblue', 'steelblue',
        'darkviolet', 'mediumorchid']

    return adata, sdiag_order


def pattern_order_lesion(adata):
    adata.obs['Pattern'] = adata.obs['Pattern'].astype('category')
    # Replace undefined with UD and non-lesional with NL
    adata.obs = adata.obs.replace({'Pattern': {"undefined": "UD"}})
    # Reorder categories
    pattern_order = ['1', '2a', '2b', '3', '4a', '4b', '5', 'UD']
    adata.obs['Pattern'] = adata.obs['Pattern'].cat.reorder_categories(pattern_order)
    adata.uns['Pattern_colors'] = ['darkorange', 'darkred', 'tomato', 'royalblue', 'mediumvioletred',
                                   'magenta', 'darkviolet', 'grey']

    return adata, pattern_order


def endotype_order_lesion(adata, obs_name):
    endotype_order = np.sort(adata.obs[obs_name].astype(str).map(lambda x: x.lstrip('E')).astype(int).unique())
    len_obsname = len(endotype_order)

    if obs_name == 'Endotypes':
        endotype_order = ['E' + s for s in list(endotype_order.astype(str))]
        adata.obs[obs_name] = adata.obs[obs_name].cat.reorder_categories(endotype_order)
    else:
        adata.obs[obs_name] = adata.obs[obs_name].cat.reorder_categories(endotype_order)

    try:
        adata.uns["{}_colors".format(obs_name)] = ["#006ddb", "#b6dbff", "#004949", "#009292", "#ff6db6", "#490092",
                                                   "#b66dff", "#000000", "#920000", "#E69F00", "#D55E00", "#8B4513",
                                                   "#999999"][:len_obsname]
    except IndexError:
        adata.uns["{}_colors".format(obs_name)] = sc.pl.palettes.default_102[:len_obsname]

    return adata, endotype_order
