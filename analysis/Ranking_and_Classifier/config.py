import os
import platform 
from easydict import EasyDict as edict

if platform.system() == "Darwin" :
    dataset_root = "/Users/heydar/Work/void/data/bio/brain"
else :
    dataset_root = "/ssd/data/bio/brain"

cache_root = './'

#dataset_name = "Lesion_RNAseq_20210720_patient_meta_data_v04__CH.h5"
#dataset_name = "Lesion_RNAseq_20210720_patient_meta_data_v04__CH_20220209_TherapyResponse_v04_v220715.h5" 
dataset_name = "Lesion_RNAseq_20210720_patient_meta_data_v04__CH_20220209_TherapyResponse_v04_v220715_TherapyResponseCorrected__Endotypes_230620.h5"


class Config :
    def make_directories( self ):
        def check_create_dir( p ):
            if not os.path.exists( p ):
                os.mkdir( p )


        plots_dir = os.path.join(cache_root,'plots')
        data_dir = os.path.join(cache_root,'data')

        check_create_dir( plots_dir )
        check_create_dir( data_dir )

    def __init__( self, args ):
        self.data_path = os.path.join(dataset_root, dataset_name)
        self.seed = 12345

        self.rank_model = args.rank_model 
        self.feat_model = args.feat_model

        self.parts_ratio = 0.7
        self.parts_num = 20 

        self.sampling_ratio = 0.7
        self.sampling_num = 10

        self.genes_per_endotype =  0

        self.dnn = edict()
        self.dnn.batch_size = 16
        self.dnn.niter = 5000

        self.make_directories()

        self.path = edict()
        self.path.plots = f'plots/plot_endotype_{args.rank_model}_%s.pdf'
        self.path.ranking = f'data/ranked_genes_endotype_{args.rank_model}.pkl'
        self.path.feats = f'data/feats_endotype_{args.rank_model}_{args.feat_model}.pkl'

        self.endotypes = ['E1','E2','E3','E4','E5','E6','E7','E8','E9','E10','E11','E12','E13'] 
        self.patterns = []
