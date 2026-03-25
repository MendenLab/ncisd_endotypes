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

class Module :
    def __init__( self ):
        args = edict()
        args.rank_model = 'LinearRegression'
        args.feat_model = None 

        self.cfg = Config(args)

        self.adata = scanpy.read(self.cfg.data_path)

    def get_counts( self ):
        df = self.adata.to_df() 
        counts = pd.DataFrame( self.adata.layers['counts'], index=df.index, columns=df.columns )
        return counts

    def get_norm( self ):
        return self.adata.to_df()

    def do_stuff( self ):
        counts = self.get_counts()
        #norm = self.get_norm().to_numpy()

        norm = counts.to_numpy().T

        means = np.mean( norm, axis=1 )
        stds = np.std( norm, axis=1 )

        sinds = np.argsort(stds)[::-1]

        print( sinds )

        pp.figure()

        for ii, idx in enumerate(sinds):
            print( norm.shape, means.shape, idx )
            vals = norm[idx]

            #m = means[idx] 
            
            #vals = vals - m 
            pp.plot(np.ones_like(vals)*ii, vals, 'b.')

        pp.grid()
        pp.xlabel("Gene index")
        pp.ylabel("Normalized Gene Expression")
        pp.savefig("expressions_raw.png")
