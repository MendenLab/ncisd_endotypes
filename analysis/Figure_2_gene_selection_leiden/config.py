from easydict import EasyDict as edict
import os
from datetime import date
import numpy as np
import platform

if platform.system() == 'Linux':
    data_root = '/home/christina/2022_BRAIN_Eyerich'
else:
    data_root = '/Volumes/CH__data/Projects/Eyerich_AG_projects/BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer/analysis'


class config :
    def __init__(self):
        self.dataset = edict()
        self.dataset.path = os.path.join(data_root, 'Molecular_subtypes', 'input', 'h5_files', 'LESION')
        self.dataset.file_path = os.path.join(
            self.dataset.path,
            'Lesion_RNAseq_20210720_patient_meta_data_v04__CH_20220209_TherapyResponse_v04_v220715_TherapyResponseCorrected__Endotypes_230620.h5')
        # for application
        self.dataset.patterns_to_keep = ['1', '2a', '2b', '3', '4a', '4b', '5', 'UD']
        # for train and test
        # self.dataset.patterns_to_keep = ['1', '2a', '3', '5']
        self.dataset.normalise = 'edgeR'
        self.dataset.batch_keys = ['batchID', 'Sex_x']  # can be also None

        self.part = edict()
        self.part.ratio = 0.80  # Train ratio

        self.feature_selection = edict()
        self.feature_selection.do_feature_selection = True
        self.feature_selection.normalization = 'mean'
        self.feature_selection.nclusters = 50
        self.feature_selection.genepercentage_keep = 74  # 7 opt value for HRG
        self.feature_selection.method = 'std'  # 'std' # 'HRG' # 'HVG'
        self.feature_selection.sort_rows = True

        self.cluster = edict()
        self.cluster.method = 'leiden'
        self.cluster.resolution = 0.9
        self.cluster.nearest_neighbors = 4  # TOTO optimise value

        self.results = edict()
        path_output = os.path.join(os.path.join(
            data_root, 'Molecular_subtypes', 'output', self.cluster.method.capitalize(), str(date.today())))
        if self.feature_selection.do_feature_selection is None:
            save_folder = os.path.join(
                os.path.join(
                    path_output, "resolution_{}__dofeatureselection_{}__percentage{}".format(
                        str(self.cluster.resolution), str(self.feature_selection.do_feature_selection),
                        str(self.feature_selection.genepercentage_keep))))
        else:
            save_folder = os.path.join(
                path_output, "resolution_{}__percentage{}__featureselection_{}".format(
                    str(self.cluster.resolution), str(self.feature_selection.genepercentage_keep),
                    str(self.feature_selection.method)))

        self.results.save_folder = save_folder

        # Set seed
        self.seeds = 0


class GridSearchParameters(config):
    def __init__(self):
        # In order tp overwrite the variables from parent class you have to first call super().__init__()
        super().__init__()

        # Overwrite variables for GirdSearch
        # start_percentage = [0.5]
        # start_percentage.extend(list(np.round(np.linspace(1, 100, 100))))
        # self.feature_selection.genepercentage_keep = start_percentage
        # For application
        # self.feature_selection.genepercentage_keep = np.asarray([7])
        self.feature_selection.genepercentage_keep = np.asarray([74])

        # self.seeds = [0, 42, 50, 100, 3210, 500, 300, 5000, 123]
        # For application
        self.seeds = [0]

        # self.cluster.resolution = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
        # self.cluster.resolution = np.array([0.5, 0.8, 0.9])
        # self.cluster.resolution = np.array([0.5, 0.6, 0.7, 0.8, 0.9])
        # For application
        self.cluster.resolution = np.array([0.9])

        self.njobs = 1  # 8

        self.results = edict()
        path_output = os.path.join(os.path.join(
            data_root, 'Molecular_subtypes', 'output', self.cluster.method.capitalize(), str(date.today())))

        if self.feature_selection.do_feature_selection is None:
            save_folder = os.path.join(
                os.path.join(
                    path_output, "resolution_{}__dofeatureselection_{}".format(
                        str(self.cluster.resolution), str(self.feature_selection.do_feature_selection))))
        else:
            save_folder = os.path.join(
                path_output, "resolution_{}__featureselection_{}".format(str(self.cluster.resolution),
                                                                         str(self.feature_selection.method)))

        self.results.save_folder = save_folder


