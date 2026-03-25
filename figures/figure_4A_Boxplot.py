import matplotlib.pyplot as plt
import matplotlib as mpl
import statsmodels.stats.multitest
from statannotations.Annotator import Annotator

import seaborn as sns
import scanpy as sc

from datetime import date
import os
import pandas as pd
import numpy as np
import scipy.stats as stats

plt.rcParams["font.family"] = "Arial"
fileformat = '.pdf'


def cohens_d(d1, d2):
    """ function to calculate Cohen's d for independent samples

    Parameters
    ----------
    d1
    d2

    Returns
    -------

    """
    # calculate the size of samples
    n1, n2 = len(d1), len(d2)
    # calculate the variance of the samples
    s1, s2 = np.var(d1, ddof=1), np.var(d2, ddof=1)
    # calculate the pooled standard deviation
    s = np.sqrt(((n1 - 1) * s1 + (n2 - 1) * s2) / (n1 + n2 - 2))
    # calculate the means of the samples
    u1, u2 = np.mean(d1), np.mean(d2)
    # calculate the effect size
    return (u1 - u2) / s


def main():
    save_folder = '/Volumes/CH__data/Thesis/writing/figures/Chapter_4/therapy_response'
    df = pd.read_excel('/Volumes/CH__data/Projects/Eyerich_AG_projects/BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer/coworkers/Martin/Therapy_response_week_12_20_with_umap_E8_11_12_13.xlsx')
    # df = pd.read_excel('/Volumes/CH__data/Thesis/writing/figures/Chapter_4/files/Therapy_response_supplements_with_umap.xlsx')

    df_boxplot = df[['Delta PGA Rate 12', 'Endotypes', 'Therapeutic target']]

    # Drop Endotypes not beein associated with psoriasis
    df_boxplot = df_boxplot.loc[df_boxplot['Endotypes'].isin(['E8', 'E11', 'E12', 'E13']), :]

    PROPS = {
        'boxprops': {'edgecolor': 'black'},
        'medianprops': {'color': 'black'},
        'whiskerprops': {'color': 'black'},
        'capprops': {'color': 'black'}
    }

    # hatches must equal the number of hues (3 in this case)
    hatches = ['//', '..', 'xx']

    df_boxplot['Endotypes'] = df_boxplot['Endotypes'].astype('category')

    pairs = [
        (("E11", "IL-17"), ("E11", "IL-23")), (("E11", "IL-17"), ("E11", "TNF")), (("E11", "IL-23"), ("E11", "TNF")),
        (("E12", "IL-17"), ("E12", "IL-23")), (("E12", "IL-17"), ("E12", "TNF")), (("E12", "IL-23"), ("E12", "TNF")),
        (("E13", "IL-17"), ("E13", "IL-23")), (("E13", "IL-17"), ("E13", "TNF")), (("E13", "IL-23"), ("E13", "TNF"))]

    fig, ax = plt.subplots(figsize=(6, 6))
    sns.boxplot(df_boxplot, x='Endotypes', y='Delta PGA Rate 12', hue='Therapeutic target',
                ax=ax, showfliers=False, order=['E8', 'E11', 'E12', 'E13'], **PROPS)
    sns.swarmplot(data=df_boxplot, x='Endotypes', y='Delta PGA Rate 12', hue='Therapeutic target',
                  palette='dark:black', ax=ax, dodge=True, legend=False, size=5, order=['E8', 'E11', 'E12', 'E13'])

    annotator = Annotator(ax, pairs, data=df_boxplot, x='Endotypes', y='Delta PGA Rate 12', hue='Therapeutic target',
                          order=['E8', 'E11', 'E12', 'E13'])
    annotator.configure(test='Mann-Whitney', text_format='star', loc='outside')  # simple
    # annotator.set_custom_annotations(["first pair", "second pair", "third pair"])
    # annotator.annotate()
    annotator.apply_and_annotate()

    # select the correct patches
    patches = [patch for patch in ax.patches if type(patch) == mpl.patches.PathPatch]
    # the number of patches should be evenly divisible by the number of hatches
    h = hatches * (len(patches) // len(hatches))
    # iterate through the patches for each subplot
    for patch, hatch in zip(patches, h):
        patch.set_hatch(hatch)
        fc = patch.get_facecolor()
        patch.set_edgecolor(fc)
        patch.set_facecolor('none')

    ax.set_ylabel(r'$\Delta$PGA score', fontsize=18)
    ax.set_xlabel('', fontsize=18)
    plt.setp(ax.get_xticklabels(), rotation=0, fontsize=18)
    plt.setp(ax.get_yticklabels(), fontsize=14)
    sns.despine(fig=fig, ax=ax, trim=False)
    ax.axhline(y=0.5, ls='dashed', c='black', label="Success")

    labels = ["IL-17", "TNF", "IL-23", "Success"]
    handles, _ = ax.get_legend_handles_labels()
    l = ax.legend(title="Drug target", loc='center left', bbox_to_anchor=(1, 0.5),
                  ncol=1, prop={'size': 14}, frameon=False, title_fontsize=16, handles=handles, labels=labels)
    for lp, hatch in zip(l.get_patches(), hatches):
        lp.set_hatch(hatch)
        fc = lp.get_facecolor()
        lp.set_edgecolor(fc)
        lp.set_facecolor('none')

    plt.tight_layout()
    plt.savefig(os.path.join(save_folder, 'Boxplot_drug_targets.pdf'), bbox_inches='tight')
    plt.close('all')

    # For thesis stats
    # df = df[df['Diag'] == 'psoriasis']
    df = df[df['Therapeutic target'].isin(['IL-17', 'TNF', 'IL-23'])]
    adata = sc.read(os.path.join(
        '/Volumes', 'CH__data', 'Projects', 'Eyerich_AG_projects', 'BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer',
        'analysis', 'Molecular_subtypes', 'input', 'h5_files', 'LESION',
        'Lesion_RNAseq_20210720_patient_meta_data_v04__CH_20220209_TherapyResponse_v04_v220715_TherapyResponseCorrected__Endotypes_230620.h5'))
    adata.obs['MUC ID'] = adata.obs.index
    adata = adata[adata.obs['MUC ID'].isin(df['MUC ID'].values)]
    df_adata = adata.obs[['MUC ID', 'age', 'Sex_x', 'batchID', 'Endotypes', 'Pseudo ID']]
    # adata = adata[adata.obs['Pseudo ID'].isin(df['Pseudo ID'].values)]
    # df_adata = adata.obs[['Pseudo ID', 'age', 'Sex_x', 'batchID', 'Endotypes']]

    df = df[df['MUC ID'].isin(df_adata['MUC ID'].values)]
    assert np.all(df['MUC ID'].unique() == df_adata['MUC ID'].values), 'Nope'

    df['Sex_x'] = 'unknown'
    df['age'] = 0
    for patient_id in df['MUC ID'].unique():
        mask = df['MUC ID'].isin([patient_id])
        mask_adata = df_adata['MUC ID'].isin([patient_id])
        if np.count_nonzero(mask) > 1:
            df.loc[mask, 'Sex_x'] = [df_adata.loc[mask_adata, :]['Sex_x'].astype(str).values[0]] * np.count_nonzero(mask)
            df.loc[mask, 'age'] = [df_adata.loc[mask_adata, :]['age'].astype(float).values[0]] * np.count_nonzero(mask)
        else:
            df.loc[mask, 'Sex_x'] = df_adata.loc[mask_adata, :]['Sex_x'].astype(str).values[0]
            df.loc[mask, 'age'] = df_adata.loc[mask_adata, :]['age'].astype(float).values[0]

    # Read out patients where we have GEx information from
    df = df[df['Pseudo ID'].isin(df_adata['Pseudo ID'].values)]
    df['Sex_x'] = 'unknown'
    df['age'] = 0
    # assign sex and age to patients
    for patient in df['Pseudo ID'].unique():
        mask = df['Pseudo ID'].isin([patient])
        mask_adata = df_adata['Pseudo ID'].isin([patient])
        if np.count_nonzero(mask) > 1:
            df.loc[mask, 'Sex_x'] = [df_adata.loc[mask_adata, :]['Sex_x'].astype(str).values[0]] * np.count_nonzero(mask)
            df.loc[mask, 'age'] = [df_adata.loc[mask_adata, :]['age'].astype(float).values[0]] * np.count_nonzero(mask)
        else:
            df.loc[mask, 'Sex_x'] = df_adata.loc[mask_adata, :]['Sex_x'].astype(str).values[0]
            df.loc[mask, 'age'] = df_adata.loc[mask_adata, :]['age'].astype(float).values[0]

    # Subset to psoriasis
    df = df[df['Diag'].isin(['psoriasis'])]
    df['Therapeutic target'].value_counts()

    # for drug_target in ['IL-23', 'TNF', 'IL-17']:
    #     drug_id = df[['Pseudo ID', 'Therapeutic target']].loc[df['Therapeutic target'] == drug_target, 'Pseudo ID'].values
    #     for sex in ['m', 'f']:
    #         val_mean = df_adata.loc[df_adata['Pseudo ID'].isin(drug_id) & (df_adata['Sex_x'] == sex), 'age'].astype(float).mean()
    #         std_mean = df_adata.loc[df_adata['Pseudo ID'].isin(drug_id) & (df_adata['Sex_x'] == sex), 'age'].astype(float).std()
    #         print("{} {}: {:.2f} pm {:.2f}".format(drug_target, sex, val_mean, std_mean))
    #
    # for drug_target in ['IL-23', 'TNF', 'IL-17']:
    #     drug_id = df[['Pseudo ID', 'Therapeutic target']].loc[df['Therapeutic target'] == drug_target, 'Pseudo ID'].values
    #     for sex in ['m', 'f']:
    #         donors = df_adata.loc[df_adata['Pseudo ID'].isin(drug_id) & (df_adata['Sex_x'] == sex), :].shape
    #         print("{} {}: {}".format(drug_target, sex, donors))

    for drug_target in ['IL-23', 'TNF', 'IL-17']:
        drug_id = df['Therapeutic target'] == drug_target
        for sex in ['m', 'f']:
            donors = len(df.loc[drug_id & (df['Sex_x'] == sex), 'Pseudo ID'].unique())
            print("Patient {} {}: {}".format(drug_target, sex, donors))

    for drug_target in ['IL-23', 'TNF', 'IL-17']:
        drug_id = df['Therapeutic target'] == drug_target
        for sex in ['m', 'f']:
            # week 0: therapy start
            mean_val = df.loc[drug_id & (df['Sex_x'] == sex), 'week 0'].mean()
            std_val = df.loc[drug_id & (df['Sex_x'] == sex), 'week 0'].std()
            print("DeltaPGA score {} {}: {:.2f} \pm {:.2f}".format(drug_target, sex, mean_val, std_val))

    for drug_target in ['IL-23', 'TNF', 'IL-17']:
        drug_id = df['Therapeutic target'] == drug_target
        for sex in ['m', 'f']:
            # 12 weeks after initiation
            mean_val = df.loc[drug_id & (df['Sex_x'] == sex), '12 weeks after initiation '].mean()
            std_val = df.loc[drug_id & (df['Sex_x'] == sex), '12 weeks after initiation '].std()
            print("DeltaPGA score  {} {}: {:.2f} \pm {:.2f}".format(drug_target, sex, mean_val, std_val))

    for drug_target in ['IL-23', 'TNF', 'IL-17']:
        drug_id = df['Therapeutic target'] == drug_target
        for sex in ['m', 'f']:
            # 20 weeks after initiation
            mean_val = df.loc[drug_id & (df['Sex_x'] == sex), '20 weeks after initiation '].mean()
            std_val = df.loc[drug_id & (df['Sex_x'] == sex), '20 weeks after initiation '].std()
            print("DeltaPGA score  {} {}: {:.2f} \pm {:.2f}".format(drug_target, sex, mean_val, std_val))

    for drug_target in ['IL-23', 'TNF', 'IL-17']:
        drug_id = df['Therapeutic target'] == drug_target
        for sex in ['m', 'f']:
            mean_val = df.loc[drug_id & (df['Sex_x'] == sex), 'age'].mean()
            std_val = df.loc[drug_id & (df['Sex_x'] == sex), 'age'].std()
            print("Age of {} {}: {:.2f} \pm {:.2f}".format(drug_target, sex, mean_val, std_val))

    for drug_target in ['IL-23', 'TNF', 'IL-17']:
        drug_id = df['Therapeutic target'] == drug_target
        for diag in ['psoriasis', 'pityriasis rubra pilaris', 'keratosis lichenoides chronica']:
            donors = len(df.loc[drug_id & (df['Diag'] == diag), 'Pseudo ID'].unique())
            print("Patient {} {}: {}".format(drug_target, diag, donors))

    # # Get duplicates
    # df.loc[df['Pseudo ID'].isin(df.loc[df['Pseudo ID'].duplicated(), 'Pseudo ID'].values), :][
    #     ['MUC ID', 'Therapeutic target', 'Diag']]

    boxplot = df.loc[df['Pseudo ID'].isin(df.loc[df['Pseudo ID'].duplicated(), 'Pseudo ID'].values), :][
        ['Pseudo ID', 'Delta PGA Rate 12', 'Therapeutic target', 'Endotypes']]
    sns.boxplot(data=boxplot, x='Endotypes', y='Delta PGA Rate 12', hue='Therapeutic target')

    # composition drug targets by endotypes
    pd.crosstab(df['Therapeutic target'], df['Endotypes'], normalize='index').plot(kind="bar", stacked=True, rot=0)
    plt.close('all')

    # pd.melt(df[['Therapeutic target', 'Endotypes', 'Delta PGA Rate 12']], id_vars=['Therapeutic target', 'Endotypes'])

    # Get mean and std per endotype and drug target
    for endotype in ['E8', 'E11', 'E12', 'E13']:
        for drug_target in ['IL-23', 'IL-17', 'TNF']:
            mean_val = df.loc[(df['Endotypes'] == endotype) & (df['Therapeutic target'] == drug_target), 'Delta PGA Rate 12'].mean()
            std_val = df.loc[(df['Endotypes'] == endotype) & (df['Therapeutic target'] == drug_target), 'Delta PGA Rate 12'].std()
            median_val = df.loc[(df['Endotypes'] == endotype) & (df['Therapeutic target'] == drug_target), 'Delta PGA Rate 12'].median()
            print("{} {}: {:.2f} \pm {:.2f}".format(endotype, drug_target, mean_val, std_val))
            print('{} {}: {} \n'.format(endotype, drug_target, median_val))

    # Log2FC of IL-23 vs IL-17 and TNF in E12
    pvals = []
    for endotype in ['E8', 'E11', 'E12', 'E13']:
        df_endotype = df.loc[df['Endotypes'] == endotype, :][['Therapeutic target', 'Delta PGA Rate 12']]
        log2fc_IL23_IL17 = np.log2((df_endotype.loc[df_endotype['Therapeutic target'] == 'IL-23', 'Delta PGA Rate 12'].mean()/
                                   df_endotype.loc[df_endotype['Therapeutic target'] == 'IL-17', 'Delta PGA Rate 12'].mean()))
        log2fc_IL23_TNF = np.log2((df_endotype.loc[df_endotype['Therapeutic target'] == 'IL-23', 'Delta PGA Rate 12'].mean()/
                                   df_endotype.loc[df_endotype['Therapeutic target'] == 'TNF', 'Delta PGA Rate 12'].mean()))
        log2fc_IL17_TNF = np.log2((df_endotype.loc[df_endotype['Therapeutic target'] == 'IL-17', 'Delta PGA Rate 12'].mean()/
                                   df_endotype.loc[df_endotype['Therapeutic target'] == 'TNF', 'Delta PGA Rate 12'].mean()))

        cd_val_IL23_IL17 = cohens_d(
            d1=df_endotype.loc[df_endotype['Therapeutic target'] == 'IL-23', 'Delta PGA Rate 12'],
            d2=df_endotype.loc[df_endotype['Therapeutic target'] == 'IL-17', 'Delta PGA Rate 12'])
        cd_val_IL23_TNF = cohens_d(
            d1=df_endotype.loc[df_endotype['Therapeutic target'] == 'IL-23', 'Delta PGA Rate 12'],
            d2=df_endotype.loc[df_endotype['Therapeutic target'] == 'TNF', 'Delta PGA Rate 12'])
        cd_val_IL17_TNF = cohens_d(
            d1=df_endotype.loc[df_endotype['Therapeutic target'] == 'IL-17', 'Delta PGA Rate 12'],
            d2=df_endotype.loc[df_endotype['Therapeutic target'] == 'TNF', 'Delta PGA Rate 12'])

        _, pval_il23_vs_il17 = stats.mannwhitneyu(
            df_endotype.loc[df_endotype['Therapeutic target'] == 'IL-23', 'Delta PGA Rate 12'],
            df_endotype.loc[df_endotype['Therapeutic target'] == 'IL-17', 'Delta PGA Rate 12'])
        _, pval_il23_vs_tnf = stats.mannwhitneyu(
            df_endotype.loc[df_endotype['Therapeutic target'] == 'IL-23', 'Delta PGA Rate 12'],
            df_endotype.loc[df_endotype['Therapeutic target'] == 'TNF', 'Delta PGA Rate 12'])
        _, pval_il17_vs_tnf = stats.mannwhitneyu(
            df_endotype.loc[df_endotype['Therapeutic target'] == 'IL-17', 'Delta PGA Rate 12'],
            df_endotype.loc[df_endotype['Therapeutic target'] == 'TNF', 'Delta PGA Rate 12'])
        pvals.extend([pval_il23_vs_il17, pval_il23_vs_tnf, pval_il17_vs_tnf])

        n_IL23_IL17 = df_endotype.loc[df_endotype['Therapeutic target'].isin(['IL-23', 'IL-17']), 'Delta PGA Rate 12'].shape[0]
        n_IL23_TNF = df_endotype.loc[df_endotype['Therapeutic target'].isin(['IL-23', 'TNF']), 'Delta PGA Rate 12'].shape[0]
        n_IL17_TNF = df_endotype.loc[df_endotype['Therapeutic target'].isin(['IL-17', 'TNF']), 'Delta PGA Rate 12'].shape[0]

        print('{} n={} IL-23 vs IL-17: p:{:.2f} log2FC: {:.2f} Cohens D: {:.2f}'.format(
            endotype, n_IL23_IL17, pval_il23_vs_il17, log2fc_IL23_IL17, cd_val_IL23_IL17))
        print('{} n={} IL-23 vs TNF: p:{:.2f} log2FC: {:.2f} Cohens D: {:.2f}'.format(
            endotype, n_IL23_TNF, pval_il23_vs_tnf, log2fc_IL23_TNF, cd_val_IL23_TNF))
        print('{} n={} IL-17 vs TNF: p:{:.2f} log2FC: {:.2f} Cohens D: {:.2f}'.format(
            endotype, n_IL17_TNF, pval_il17_vs_tnf, log2fc_IL17_TNF, cd_val_IL17_TNF))

    _, padj_in_endotypes, _, _ = statsmodels.stats.multitest.multipletests(pvals=pvals[3:], method='fdr_bh')

    # Compare stats for therapy response hypothesis
    df_E12E13_il17 = df.loc[df['Endotypes'].isin(['E12', 'E13']) & (df['Therapeutic target'] == 'IL-17'), 'Delta PGA Rate 12']
    df_E8E11_il17 = df.loc[df['Endotypes'].isin(['E8', 'E11']) & (df['Therapeutic target'] == 'IL-17'), 'Delta PGA Rate 12']
    df_E12E13_il23 = df.loc[df['Endotypes'].isin(['E12', 'E13']) & (df['Therapeutic target'] == 'IL-23'), 'Delta PGA Rate 12']
    df_E8E11_il23 = df.loc[df['Endotypes'].isin(['E8', 'E11']) & (df['Therapeutic target'] == 'IL-23'), 'Delta PGA Rate 12']
    df_E12 = df.loc[df['Endotypes'].isin(['E12']) & (df['Therapeutic target'] == 'TNF'), 'Delta PGA Rate 12']
    df_E13 = df.loc[df['Endotypes'].isin(['E13']) & (df['Therapeutic target'] == 'TNF'), 'Delta PGA Rate 12']
    _, pval_il17 = stats.mannwhitneyu(df_E12E13_il17, df_E8E11_il17)
    _, pval_il23 = stats.mannwhitneyu(df_E12E13_il23, df_E8E11_il23)
    _, pval_tnf = stats.mannwhitneyu(df_E12, df_E13)

    log2fc_IL17 = np.log2((df_E8E11_il17.mean()/df_E12E13_il17.mean()))
    log2fc_IL23 = np.log2((df_E12E13_il23.mean()/df_E8E11_il23.mean()))
    log2fc_TNF = np.log2((df_E13.mean()/df_E12.mean()))

    cd_val_IL17 = cohens_d(d1=df_E8E11_il17, d2=df_E12E13_il17)
    cd_val_IL23 = cohens_d(d1=df_E12E13_il23, d2=df_E8E11_il23)
    cd_val_TNF = cohens_d(d1=df_E13, d2=df_E12)

    # _, padj, _, _ = statsmodels.stats.multitest.multipletests(pvals=[pval_il17, pval_il23, pval_tnf], method='fdr_bh')

    print('{}: p:{:.2f} log2FC: {:.2f} Cohens D: {:.2f}'.format('IL-17', pval_il17, log2fc_IL17, cd_val_IL17))
    print('{}: p:{:.2f} log2FC: {:.2f} Cohens D: {:.2f}'.format('IL-23', pval_il23, log2fc_IL23, cd_val_IL23))
    print('{}: p:{:.2f} log2FC: {:.2f} Cohens D: {:.2f}'.format('TNF', pval_tnf, log2fc_TNF, cd_val_TNF))

    adata = adata[:, adata.var_names.isin(['IL23A', 'IL17A', 'TNF'])]
    df_counts = adata.to_df()
    # df_counts['sampleType'] = adata.obs['sampleType']
    # df_counts['sampleType'] = adata.obs['sampleType']
    df_counts['Endotypes'] = adata.obs['Endotypes']
    df_counts = df_counts[df_counts['Endotypes'].isin(['E8', 'E11', 'E12', 'E13'])]
    df_counts['Endotypes'] = df_counts['Endotypes'].cat.remove_unused_categories()
    df_counts_melted = pd.melt(df_counts, id_vars='Endotypes')

    pairs = [
        (("TNF", "E8"), ("TNF", "E11")), (("TNF", "E8"), ("TNF", "E12")), (("TNF", "E8"), ("TNF", "E13")),
        (("TNF", "E11"), ("TNF", "E12")), (("TNF", "E12"), ("TNF", "E13")), (("TNF", "E12"), ("TNF", "E13")),
        (("IL17A", "E8"), ("IL17A", "E11")), (("IL17A", "E8"), ("IL17A", "E12")), (("IL17A", "E8"), ("IL17A", "E13")),
        (("IL17A", "E11"), ("IL17A", "E12")), (("IL17A", "E12"), ("IL17A", "E13")), (("IL17A", "E12"), ("IL17A", "E13")),
        (("IL23A", "E8"), ("IL23A", "E11")), (("IL23A", "E8"), ("IL23A", "E12")), (("IL23A", "E8"), ("IL23A", "E13")),
        (("IL23A", "E11"), ("IL23A", "E12")), (("IL23A", "E12"), ("IL23A", "E13")), (("IL23A", "E12"), ("IL23A", "E13"))]
    fig, ax = plt.subplots(figsize=(4, 6))
    sns.boxplot(df_counts_melted, x='hgnc_symbol', y='value',  hue='Endotypes',
                palette=['#000000', '#D55E00', '#8B4513', '#999999'], ax=ax)
    annotator = Annotator(ax, pairs, data=df_counts_melted, x='hgnc_symbol', y='value', hue='Endotypes')
    annotator.configure(test='Mann-Whitney', text_format='star', loc='outside')  # simple
    annotator.apply_and_annotate()
    ax.set_ylabel('normed counts')
    ax.set_xlabel('')
    # Put a legend below current axis
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1),
              fancybox=False, shadow=False, ncol=5, frameon=False)
    sns.despine()
    plt.tight_layout()
    plt.savefig(os.path.join(save_folder, 'Boxplot_Lesion_TNF_IL23A_IL17A.pdf'), bbox_inches='tight')
    plt.close('all')


if __name__ == '__main__':
    main()