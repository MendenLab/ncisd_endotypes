import numpy as np
import multiprocessing
from joblib import Parallel, delayed
from tqdm import tqdm

from .feature_model import FeatureModel

class FeatureModels :
    def __init__( self, classifier ):
        self.classifier = classifier

        self.models = None 
        self.columns = None

    def fit( self, X, y_binary ):
        # X : rows are samples and columns are genes 
        columns = X.columns

        def fit_single( col ):
            Xcol = X.loc[:,col]

            model = FeatureModel( self.classifier )
            model.fit( Xcol, y_binary )

            return model

        self.columns = columns
        self.models = Parallel(n_jobs=multiprocessing.cpu_count())( delayed(fit_single)(col) for col in columns )

    def eval( self, X ):
        def eval_single(  idx, col ):
            Xcol = X.loc[:,col]
            model = self.models[idx]
            res = model.eval( Xcol )
            return res 

        out = []
        for idx, col in enumerate(self.columns):
            r = eval_single( idx, col )
            out.append( r )

        return np.array( out )
