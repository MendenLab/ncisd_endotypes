import pickle
import scanpy as sc
import os
import platform

if platform.system() == 'Linux':
    data_root = '/home/christina/2022_BRAIN_Eyerich'
else:
    data_root = '/Volumes/CH__data/Projects/Eyerich_AG_projects/BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer/analysis'


filename = '/Volumes/CH__data/Projects/Eyerich_AG_projects/BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer/analysis/Molecular_subtypes/output/Leiden/2024-07-02/resolution_[0.9]__featureselection_HIG/seed0__percentage7__resolution0.9/Scores_GridSearch_1_minPercentage7_maxPercentage7_minResolution0.9_maxResolution0.9_.pkl'


with open(filename, 'rb') as ff:
    gs_hig = pickle.load(ff)

# train and test sample IDs
train_ids = gs_hig['train samples']
test_ids = gs_hig['test samples']

# Load data
data_filename = 'Lesion_RNAseq_20210720_patient_meta_data_v04__CH_20220209_TherapyResponse_v04_v220715_TherapyResponseCorrected__Endotypes_230620.h5'
data_path = os.path.join(data_root, 'Molecular_subtypes', 'input', 'h5_files', 'LESION')

adata = sc.read(os.path.join(data_path, data_filename))

adata.obs['Train/Test'] = 'not used'
adata.obs.loc[adata.obs.index.isin(list(train_ids)), 'Train/Test'] = 'train'
adata.obs.loc[adata.obs.index.isin(list(test_ids)), 'Train/Test'] = 'test'


adata.obs.to_excel(
    os.path.join('/Volumes/CH__data/Projects/Eyerich_AG_projects/BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer',
                 'analysis/Molecular_subtypes/output/Leiden/2024-07-02', 'Train_test_annotation.xlsx'))

