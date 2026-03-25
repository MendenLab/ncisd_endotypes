import os
import scanpy as sc
import pandas as pd
import numpy as np

from datetime import date


def sanity_check(df_metadata, result):
    # Sanity check
    age_mean = df_metadata.loc[df_metadata['Disease'] == 'eczema', 'Age'].mean()
    age_std = df_metadata.loc[df_metadata['Disease'] == 'eczema', 'Age'].std()
    age_mean_std = "{:.1f} ± {:.1f}".format(age_mean, age_std)
    sex_percent = ((df_metadata.loc[df_metadata['Disease'] == 'eczema', 'Sex'].str.lower() == 'm').astype(
        int).mean() * 100).round(1)
    assert np.all(result.loc[result['Disease'] == 'eczema', 'Age (mean ± SD)'] == age_mean_std) & np.all(
        result.loc[result['Disease'] == 'eczema', 'Sex (% of male)'] == sex_percent), 'Eczema: Not the same'
    age_mean = df_metadata.loc[df_metadata['Disease'] == 'psoriasis', 'Age'].mean()
    age_std = df_metadata.loc[df_metadata['Disease'] == 'psoriasis', 'Age'].std()
    age_mean_std = "{:.1f} ± {:.1f}".format(age_mean, age_std)
    sex_percent = ((df_metadata.loc[df_metadata['Disease'] == 'psoriasis', 'Sex'].str.lower() == 'm').astype(
        int).mean() * 100).round(1)
    assert np.all(result.loc[result['Disease'] == 'psoriasis', 'Age (mean ± SD)'] == age_mean_std) & np.all(
        result.loc[result['Disease'] == 'psoriasis', 'Sex (% of male)'] == sex_percent), \
        'Psoriasis: Not the same'


def get_infos(df_metadata):
    # # Replace string with integer in Sex
    # mapping = {'m': 0, 'f': 1}
    # df_metadata = df_metadata.replace({'Sex': mapping})

    df_metadata['Age'] = df_metadata['Age'].astype(int)

    # Regroup by Pattern and diag
    # Gruppieren nach Pattern und diag
    df_grouped = df_metadata.groupby(['Pattern', 'Disease']).agg(
        Patients_n=('Patient', lambda x: x.nunique()),
        Age_mean=('Age', 'mean'),
        Age_std=('Age', 'std'),
        Male_percentage=('Sex', lambda x: (x.str.lower() == 'm').mean() * 100)
    ).reset_index()

    # Alter als "mean ± SD" formatieren
    df_grouped['Age (mean ± SD)'] = df_grouped.apply(
        lambda row: f"{row['Age_mean']:.1f} ± {row['Age_std']:.1f}", axis=1
    )

    # Prozentsatz männlich runden und formatieren
    df_grouped['Sex (% of male)'] = df_grouped['Male_percentage'].round(1)

    # Endauswahl der Spalten
    result = df_grouped[['Pattern', 'Disease', 'Patients_n', 'Age (mean ± SD)', 'Sex (% of male)']]

    # Optional: Spalte umbenennen
    result = result.rename(columns={'Patients_n': 'Patients (n)'})

    mask_nan = ~result.isna().sum(axis=1).astype(bool)
    result = result.loc[mask_nan]

    return result


def get_metadata(df):
    # Get Pattern, diag, Patients (n), Age (mean ± SD), Sex (% of male)
    df_metadata = df[['healthysamp_Pattern_v2', 'healthysamp_diag_v2', 'Patient ID', 'age', 'Sex_x']]
    # rename diag to Disease, Pseudo ID to Patient, age to Age, and Sex_x to Sex
    df_metadata.rename(columns={'healthysamp_Pattern_v2': 'Pattern', 'healthysamp_diag_v2': 'Disease',
                                'Patient ID': 'Patient', 'age': 'Age', 'Sex_x': 'Sex'}, inplace=True)

    return df_metadata


def rename_nonlesion_samples_to_paired_lesion_name(
        obs, ref_obsname='diag', obsname='healthysamp_diag_v2', group_name='diag'):

    obs[ref_obsname] = obs[ref_obsname].astype(str)
    obs[ref_obsname] = obs[ref_obsname].astype('category')

    levels_ref_obsname = list(obs[ref_obsname].cat.categories)
    # remove non-lesional to set it on first place
    levels_ref_obsname.remove("non-lesional")
    # Reorder levels such that non-lesional level comes first
    obs[ref_obsname] = obs[ref_obsname].cat.reorder_categories(["non-lesional"] + levels_ref_obsname)

    print(list(obs[ref_obsname].cat.categories))

    obs[obsname] = 'Unknown'
    # 40 only non-lesion samples (after filtering and everything)
    for observable in list(obs[ref_obsname].cat.categories):
        mask_observable = obs[ref_obsname] == observable
        # Find healthy partner via Patient ID == Pseudo.ID
        patient_id = obs.loc[mask_observable, 'Patient ID']
        # multiple matching using which and %in%, as match only return first encounter
        ind_patients = obs.index[obs['Patient ID'].isin(patient_id)]
        ind_patients = ind_patients[
            (obs.loc[ind_patients, group_name] == observable) | (
                    obs.loc[ind_patients, group_name] == 'non-lesional')]
        obs.loc[ind_patients, obsname] = observable

    obs[obsname] = obs[obsname].astype("str")

    return obs


def main():
    # Load L bulk samples
    adata = sc.read(os.path.join(
        '/Volumes', 'CH__data', 'Projects', 'Eyerich_AG_projects', 'BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer',
        'analysis', 'Molecular_subtypes', 'input', 'h5_files', 'LESION',
        'Lesion_RNAseq_20210720_patient_meta_data_v04__CH_20220209_TherapyResponse_v04_v220715_TherapyResponseCorrected__Endotypes_230620.h5'))

    # PatientID sometimes contains HG16.711_H and HG16.711_Lym -> remove rubstring to make same patient
    adata.obs['Patient ID'] = adata.obs['Pseudo ID'].str.split('_').str[0]

    df_metadata = get_metadata(df=adata.obs)
    print(df_metadata.columns.to_list())
    LS_result = get_infos(df_metadata=df_metadata)
    # Ausgabe
    print(LS_result)
    print('Only L samples: ', LS_result['Patients (n)'].sum())

    sanity_check(df_metadata=df_metadata, result=LS_result)

    # Save result but these are only L samples .. Not enough
    LS_result.to_excel(os.path.join(output_dir, 'Table_S1.xlsx'), index=False, sheet_name='Table_S1_LS_samples')

    # L and NL: 632 samples nach preprocessing
    adata = sc.read(os.path.join(
        '/Volumes', 'CH__data', 'Projects', 'Eyerich_AG_projects', 'BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer',
        'analysis', 'Molecular_subtypes', 'input', 'h5_files', 'NON_LESION',
        'LNL_RNAseq_20210720_patient_meta_data_v04__CH__Endotypes_230620.h5'))
    # PatientID sometimes contains HG16.711_H and HG16.711_Lym -> remove rubstring to make same patient
    adata.obs['Patient ID'] = adata.obs['Pseudo ID'].str.split('_').str[0]

    # rename non lesional to non-lesional
    adata.obs['diag'] = adata.obs['diag'].astype(str)
    adata.obs = adata.obs.replace({'diag': r'non lesional'}, {'diag': 'non-lesional'}, regex=True)
    adata.obs = rename_nonlesion_samples_to_paired_lesion_name(
        obs=adata.obs, ref_obsname='diag', obsname='healthysamp_diag_v2', group_name='diag')
    # rename non lesional to non-lesional
    adata.obs['Pattern'] = adata.obs['Pattern'].astype(str)
    adata.obs = adata.obs.replace({'Pattern': r'NL'}, {'Pattern': 'non-lesional'}, regex=True)
    adata.obs = rename_nonlesion_samples_to_paired_lesion_name(
        obs=adata.obs, ref_obsname='Pattern', obsname='healthysamp_Pattern_v2', group_name='Pattern')

    # read out unique patients (patient ID),
    df = adata.obs.drop_duplicates(subset=['Patient ID'], keep='last')

    df_metadata = get_metadata(df=df)
    print(df_metadata.columns.to_list())
    result = get_infos(df_metadata=df_metadata)
    # Ausgabe
    print(result)
    print('LS and NL patients: ', result['Patients (n)'].sum())

    sanity_check(df_metadata=df_metadata, result=result)

    # Save both to Table
    result.to_excel(os.path.join(output_dir, 'Table_S1.xlsx'), index=False, sheet_name='Table_S1_LSNL_samples')

    # Number of patients after preprocessing:
    print('Number of samples: ', adata.shape[0])  # 632
    print('Number of patients: ', len(np.unique(adata.obs['Patient ID'])))  # 368
    df_crosstab = pd.crosstab(adata.obs['Patient ID'], adata.obs['sampleType'])
    print('Number of patients with 1 L & 2 NL samples: ',
          df_crosstab.loc[(df_crosstab['d'] == 1) & (df_crosstab['h'] == 2), :].shape[0])  # 1
    print('Number of patients with 2 L & 1 NL samples: ',
          df_crosstab.loc[(df_crosstab['d'] == 2) & (df_crosstab['h'] == 1), :].shape[0])  # 3
    print('Number of patients with 1 L & 1 NL samples: ',
          df_crosstab.loc[(df_crosstab['d'] == 1) & (df_crosstab['h'] == 1), :].shape[0])  # 255
    print('Number of patients with only L samples: ',
          df_crosstab.loc[(df_crosstab['d'] > 0) & (df_crosstab['h'] == 0), :].shape[0])  # 82 (1 patient has 2 LS)
    print('Number of patients with only H samples: ',
          df_crosstab.loc[(df_crosstab['d'] == 0) & (df_crosstab['h'] == 1), :].shape[0])  # 27

    df_overview = pd.read_excel('/Volumes/CH__data/Projects/data/Summary/Data_Overview.xlsx', sheet_name='BRAIN')
    df_overview['SampleID'] = df_overview['SampleID'].str.replace('Sample_', '')
    sampleIDs = df_overview.loc[~df_overview['SampleID'].isna(), 'SampleID'].to_list()
    df_reindexed = adata.obs.reindex(sampleIDs)
    df_reindexed[['Patient ID', 'diag', 'Pattern', 'sampleType', 'batchID', 'age', 'Sex_x']].to_excel(
        '/Volumes/CH__data/Projects/data/Summary/MetaData_Overview.xlsx')


if __name__ == '__main__':
    output_dir = os.path.join('/Volumes', 'CH__data', 'Projects', 'Eyerich_AG_projects',
                              'BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer', 'analysis', 'Molecular_subtypes',
                              'output', 'Table_S1_data_overview', str(date.today()))
    os.makedirs(output_dir, exist_ok=True)

    main()
