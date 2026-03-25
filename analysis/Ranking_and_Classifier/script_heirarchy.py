import os
import pickle
import scanpy
import pandas as pd
import numpy as np
import sys
import yaml
import time
import argparse

from sklearn.linear_model import LinearRegression
from sklearn.svm import SVR
from easydict import EasyDict as edict
from tqdm import tqdm 

from config import Config
from dataframe_tools import partition 
from modeling import utils
from modeling import FeatureModels
from report_generator import ReportGenerator, ReportGeneratorHeirarchy

from matplotlib import pyplot as pp 
pp.ion()

def load_brain( cfg ):
    adata = scanpy.read( cfg.data_path )

    index = adata.to_df().index
    df = adata.to_df()
    pattern = adata.obs["Pattern"].to_numpy()
    endotypes = adata.obs["Endotypes"].to_numpy()

    aa = pd.DataFrame( index=index )
    aa['pattern'] = pattern 
    aa['endotypes'] = endotypes

    return df, aa

endotypes = ['E1','E2','E3','E4','E5','E6','E7','E8','E9','E10','E11','E12','E13']

class Module :
    def __init__( self, args ):
        self.cfg = Config( args )

        with open('cfg_endotype.yaml','r') as ff :
            self.experiments = yaml.safe_load(ff)['experiments']

        self.data, self.annots = load_brain( self.cfg )

        self.parts = partition( self.annots, "endotypes", endotypes, 
                                self.cfg.parts_ratio, self.cfg.parts_num,
                                self.cfg.sampling_ratio, self.cfg.sampling_num,
                                numpy_seed=self.cfg.seed )

    def create_binary_labels( self, annots ):
        labels = []

        columns = []

        for idx in range(len(endotypes)) :
            for idy in range(idx) :
                e0 = endotypes[idx]
                e1 = endotypes[idy]

                columns.append('%s-%s' % (e0,e1))

                ll = utils.create_binary_labels_np( annots, [e0], [e1] )
                ll = ll.reshape([-1,1])

                labels.append(ll)

        labels = np.concatenate(labels,axis=1)

        return pd.DataFrame( labels, index=annots.index, columns=columns )

    def create_endotype_labels( self, annots ):
        labels = []

        for endotype in endotypes :
            ll = utils.create_binary_labels_geo( annots, endotype ).transpose()
            labels.append(ll)

        labels = np.concatenate( labels, axis=1 )

        return pd.DataFrame( labels, index=annots.index, columns=endotypes )

    def prepare_data( self, part, data, labels ):
        dd_train = data.loc[ part['train'] ]
        ll_train = labels.loc[ part['train'] ]
        dd_test = data.loc[ part['test'] ]
        ll_test = labels.loc[ part['test'] ]

        return dd_train, ll_train, dd_test, ll_test

    def create_crossval_data( self ):
        labels = self.create_binary_labels( self.annots['endotypes'] )

        classifier = eval(self.cfg.rank_model)

        cross_val = []

        for idx, part in enumerate(self.parts) :
            print( idx )
            X_tr, y_tr, X_te, y_te = self.prepare_data( part, self.data, labels ) 

            model = FeatureModels( classifier, ignore_idx=0 )
            model.fit( X_tr, y_tr )

            val = model.eval(X_te).astype( np.float32 )
            cross_val.append( val )

        cross_val = np.array( cross_val, dtype=np.float32 )
        print( cross_val.shape )


        #with open('cross_val.pkl','wb') as ff :
        #    pickle.dump( cross_val, ff )

        #print('Saved to cross_val.pkl')

        self.cross_val = cross_val

    def calculate_gene_score( self, values, labels ):

        def sigmoid( z ):
            return 1/(1 + np.exp(-z))

        scores_all = []

        for part_idx, part in enumerate(self.parts) :
            labels_part = labels.loc[ part['test'] ].to_numpy()

            values_part = values[part_idx]

            scores = []

            for l_idx in range( labels_part.shape[1] ):
                ll = labels_part[:,l_idx].ravel()
                keep = np.where( ll == 1 )[0]

                vv = values_part[keep,:]
                vv = sigmoid( vv )
                vv = np.mean( vv, axis=0 ).reshape([1,-1])

                scores.append( vv )

            scores = np.concatenate( scores, axis=0 )

            scores_all.append( scores )

        scores_all = np.mean( np.array( scores_all ), axis=0 )

        return scores_all 

        #scores_all = np.matrix( scores_all )
        #ss = scores_all * scores_all.T
        #pp.matshow( ss )

        #print( values.shape )
        #print( labels )

    def calculate_gene_scores( self ):
        cross_val = self.cross_val.transpose(1,0,2,3) 
        labels = self.create_endotype_labels( self.annots['endotypes'] )

        genes = np.array(list(self.data.columns))
        index = np.array(list(self.data.index))
        data = self.data.to_numpy().transpose()
        
        stds = np.std( data, axis=1 )
        average_std = np.mean( stds )

        keep = np.where( stds >= average_std )[0]

        cross_val = cross_val[keep]
        genes = genes[keep]

        gene_idx = 150

        gene_scores = {}

        for gene_idx, gene in tqdm(enumerate(genes)):
            ss = self.calculate_gene_score( cross_val[gene_idx], labels )
            gene_scores[gene] = ss

        self.gene_scores = gene_scores 

        with open('./data/gene_scores.pkl','wb') as ff :
            pickle.dump( self.gene_scores, ff )


        return
        ss = self.calculate_gene_score( cross_val[gene_idx], labels )

        print( ss.shape )


        return

        for part_idx, part in enumerate(self.parts):
            cross_val_part = cross_val[part_idx]
            labels_part = labels.loc[ part['test'] ].to_numpy()

            for labels_idx in range(labels_part.shape[1]):
                l = labels_part[:,labels_idx]
                keep = np.where( l == 1 )[0]
                cv = cross_val_part[:,keep,:]

                
                min_values = np.min( cv, axis=1 )
                max_values = np.max( cv, axis=1 )

                print( min_values.shape, max_values.shape )

def get_module():
    args = edict()
    args.rank_model = 'LinearRegression'
    args.feat_model = None 
    return Module( args )

if __name__=="__main__" :
    m = get_module()
    m.do_stuff()
