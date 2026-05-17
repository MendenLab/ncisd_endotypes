import pandas as pd
import numpy as np
from sklearn.preprocessing import MinMaxScaler
import rpy2
import rpy2.robjects as ro
import rpy2.robjects.packages as rpackages
from rpy2.robjects import numpy2ri, pandas2ri
from rpy2.robjects.conversion import localconverter

import platform

numpy2ri.activate()
pandas2ri.activate()

if platform.system() == 'Linux':
    ro.r('library(DESeq2, lib="/home/christina/R_packages")')
else:
    deseq2 = rpackages.importr('DESeq2')

edger = rpackages.importr('edgeR')
sva = rpackages.importr('sva')
summarizedExperiment = rpackages.importr("SummarizedExperiment")
base = rpackages.importr("base")


def edger_normalise(X: np.ndarray, coldata: pd.DataFrame, batch_keys: [str, list], logcpm: bool = True):
    # 1. Create DGEList object
    # 1.1 First read out save sample ids to dataframe
    sampledata = pd.DataFrame({'sample': list(coldata.index)})
    # 1.2 Use design 1 to not make use of any design information
    y = edger.DGEList(counts=X.transpose(), samples=sampledata)

    # Library Size Normalization: Calculate normalization factors
    y = edger.calcNormFactors(y, method="TMM")

    # Estimate dispersion - robust against outliers
    ro.r.assign('coldata', sampledata)
    design = ro.r('model.matrix(~1, data=coldata)')
    y = edger.estimateDisp(y, design, robust=True)

    # Get normalized log counts
    X_normed = edger.cpm(y, normalized_lib_sizes=True, log=logcpm)

    X_normed_batchcorrected = apply_batch_correction(X_normed=X_normed, coldata=coldata, batch_keys=batch_keys)

    return X_normed_batchcorrected.T


def minmax_normalise(X: np.ndarray, coldata: pd.DataFrame, batch_keys: [str, list]):
    # Min Max scaler for normalisation
    norm_factor = MinMaxScaler()

    # Apply transformation
    X_normed = norm_factor.fit_transform(X=X)

    # Apply batch corrrection
    X_normed_batchcorrected = apply_batch_correction(X_normed=X_normed.T, coldata=coldata, batch_keys=batch_keys)

    return X_normed_batchcorrected.T


def pandasToR(df: pd.DataFrame):
    """
    Transforms a pandas counts to R
    Parameters
    ----------
    df: Pandas counts

    Returns: R counts corresponding to given pandas Dataframe
    -------

    """
    with localconverter(ro.default_converter + pandas2ri.converter):
        return ro.conversion.py2rpy(df)


def apply_batch_correction(X_normed, coldata, batch_keys=None):
    mat = X_normed

    if batch_keys is not None:
        if isinstance(batch_keys, list):
            for key in batch_keys:
                if key in coldata.columns:
                    # Convert to categories
                    coldata[key] = coldata[key].astype('category')
                    try:
                        print("Correcting for: ", key)
                        mat = sva.ComBat(dat=mat, batch=coldata[key])
                    except rpy2.rinterface_lib.embedded.RRuntimeError:
                        # System ist genau singulär: U[1,1] = 0
                        print("Singular Matrix")
        else:
            if batch_keys in coldata.columns:
                # Convert to categories
                coldata[batch_keys] = coldata[batch_keys].astype('category')
                try:
                    print("Correcting for: ", batch_keys)
                    mat = sva.ComBat(dat=mat, batch=coldata[batch_keys])
                except rpy2.rinterface_lib.embedded.RRuntimeError:
                    # System ist genau singulär: U[1,1] = 0
                    print("Singular Matrix")

    return mat