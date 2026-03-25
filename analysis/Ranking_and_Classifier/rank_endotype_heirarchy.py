from io import UnsupportedOperation
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

class1 = ['E5']
class2 = ['E6']

classifier = SVR

def eval_function( x,y ):
    return np.mean( (x-y) ** 2, axis=1 )


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

    def create_labels( self, c1, c2, all ):
        negative = []
        for a in all :
            if a not in c1 and a not in c2 :
                negative.append( a )

        blabels_c1 = utils.create_binary_labels_new( self.annots['endotypes'], c1, negative )
        blabels_c2 = utils.create_binary_labels_new( self.annots['endotypes'], c2, negative )

        return blabels_c1, blabels_c2

    def prepare_data( self, part, data, labels ):
        dd_train = data.loc[ part['train'] ]
        ll_train = labels.loc[ part['train'] ]
        dd_test = data.loc[ part['test'] ]
        ll_test = labels.loc[ part['test'] ]

        keep_train = ll_train['labels'] != 0
        keep_test = ll_test['labels'] != 0
            
        dd_train = dd_train[ keep_train ]
        ll_train = ll_train[ keep_train ]

        dd_test = dd_test[ keep_test ]
        ll_test = ll_test[ keep_test ]

        return dd_train, ll_train, dd_test, ll_test


    def do_stuff( self ):

        labels_c1, labels_c2 = self.create_labels( class1, class2, endotypes )

        cross_val_scores = []

        for idx, part in enumerate(self.parts) :
            print(f' {idx+1}/{len(self.parts)}')

            X_train_c1, y_train_c1, X_test_c1, y_test_c1 = self.prepare_data( part, self.data, labels_c1 )
            X_train_c2, y_train_c2, X_test_c2, y_test_c2 = self.prepare_data( part, self.data, labels_c2 )

            y_train_c1 = y_train_c1['labels'].to_numpy().reshape((1,-1))
            y_train_c2 = y_train_c2['labels'].to_numpy().reshape((1,-1))

            y_test_c1 = y_test_c1['labels'].to_numpy().reshape((1,-1))
            y_test_c2 = y_test_c2['labels'].to_numpy().reshape((1,-1))

            fmodels_c1 = FeatureModels( classifier )
            fmodels_c1.fit( X_train_c1, y_train_c1 )

            fmodels_c2 = FeatureModels( classifier )
            fmodels_c2.fit( X_train_c2, y_train_c2 )

            eval_c2 = fmodels_c1.eval( X_test_c2 )
            eval_c1 = fmodels_c2.eval( X_test_c1 )

            y_test_c1 = y_test_c1.T
            y_test_c1 = np.tile(y_test_c1, (len(eval_c1), 1, 1))
            
            y_test_c2 = y_test_c2.T
            y_test_c2 = np.tile(y_test_c2, (len(eval_c2), 1, 1))

            scores_c1 = eval_function( eval_c1, y_test_c1 )
            scores_c2 = eval_function( eval_c2, y_test_c2 )

            scores = (scores_c1 + scores_c2) * 0.5  
            scores = np.expand_dims( scores, axis=0 )

            cross_val_scores.append( scores )

        cross_val_scores = np.concatenate( cross_val_scores, axis=0 )
        cross_val_scores = np.mean( cross_val_scores, axis=0 )

        self.cross_val_scores = cross_val_scores 

    def do_report( self ):
        cross_val_scores = self.cross_val_scores

        scores = cross_val_scores[:,0].ravel()
        inds = np.argsort(scores)

        feats = np.array(list(self.data.columns))[inds[:50]]

        inds = inds[:50]

        groups = {}
        groups['all'] = endotypes 
        groups['positives'] = class1 + class2

        report = ReportGeneratorHeirarchy( self.data, self.annots['endotypes'], inds, scores[inds], groups )
        report.plot()


        #for idx, endotype in enumerate( endotypes ):
        #    scores = cross_val_scores[:,idx].ravel()
        #    inds = np.argsort( scores )
        #    feats = np.array(list(self.data.columns))[inds[:1000]]

        #    inds = inds[:50]

        #    report = ReportGenerator( self.data, self.annots['endotypes'], inds, scores[inds], endotypes )
        #    report.plot( "Brain - Endotype", 'plot.pdf'  )

        #pass

def get_module():
    args = edict()
    args.rank_model = 'LinearRegression'
    args.feat_model = None 
    return Module( args )

def parse_commandline( rank_model_required=False, feat_model_required=False ):

    parser = argparse.ArgumentParser()
    parser.add_argument( '--rank-model', 
                         help="Modeling architecture for ranking", 
                         type=str,
                         required=rank_model_required )
    parser.add_argument( '--feat-model', 
                         help="Modeling architecture for feature extraction", 
                         type=str,
                         required=feat_model_required )
    args = parser.parse_args() 

    return args

if __name__=="__main__" :
    args = parse_commandline( rank_model_required=True, 
                              feat_model_required=False )

    m = Module( args )
    #m.do_stuff()
