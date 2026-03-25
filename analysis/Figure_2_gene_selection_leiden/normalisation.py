# conda install -c bioconda bioconductor-deseq2

import pandas as pd
import numpy as np
from sklearn.preprocessing import MinMaxScaler
import rpy2
import rpy2.robjects as ro
import rpy2.robjects.packages as rpackages
from rpy2.robjects import numpy2ri, pandas2ri, Formula
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


class Normalisation:
    def __init__(self):
        self.dispersion = None
        self.sizefactor = None
        self.norm_factor = None

    def normalise(self, X: np.array, normalisation_method, parts, coldata: pd.DataFrame, batch_keys: [str, list]):
        if normalisation_method == 'DESeq2':
            if platform.system() == 'Linux':
                X_normed_bc = self.deseq2_normalise_linux(X=X, parts=parts, coldata=coldata, batch_keys=batch_keys)
            else:
                X_normed_bc = self.deseq2_normalise(X=X, parts=parts, coldata=coldata, batch_keys=batch_keys)
        elif normalisation_method == 'MinMax':
            X_normed_bc = self.minmax_normalise(X=X, parts=parts, coldata=coldata, batch_keys=batch_keys)
        elif normalisation_method == 'edgeR':
            X_normed_bc = self.edger_normalise(X=X, parts=parts, coldata=coldata, batch_keys=batch_keys)
        else:
            return X
        return X_normed_bc

    def deseq2_normalise_linux(self, X: np.ndarray, parts, coldata: pd.DataFrame, batch_keys: [str, list]):
        """ Normalise data using DESeq2

        Parameters
        ----------
        X : Filtered, raw count matrix

        Returns
        -------

        """
        # 1. Create DESeq2 object
        # 1.1 First read out save sample ids to dataframe
        # coldata = self._parts.sample_names()
        sample_id = pd.DataFrame({'sample': coldata})
        ro.r.assign('ncols', sample_id)
        # Transpose the count matrix
        ro.r.assign('counts', X.T)
        # 1.2 Use design 1 to not make use of any design information
        ro.r('dds <- DESeq2::DESeqDataSetFromMatrix(countData=counts, colData=ncols, design=~1)')

        if parts._model_status != 'test':
            # Calculate sizeFactors on whole train or application set
            # remove sizeFactors on train set, otherwise dispersion will be calculated using
            # the existing sizeFactors which were calculate on whole dataset
            # dds.sizeFactor = None
            ro.r('dds <- DESeq2::estimateSizeFactors(dds, type="ratio", quiet=TRUE)')

            #  Calculate dispersion using precalculated sizefactors
            if self.dispersion is None:
                # Calculate dispersion only once
                ro.r('dds <- DESeq2::estimateDispersions(dds)')
                self.dispersion = ro.r('DESeq2::dispersionFunction(dds)')
                ro.r('dds <- DESeq2::DESeq(dds, quiet=TRUE)')
            else:
                ro.r.assign('dispersion', self.dispersion)
                ro.r('DESeq2::dispersionFunction(dds) <- dispersion')

            # Apply Variance stabilisation
            vst = ro.r('DESeq2::varianceStabilizingTransformation(dds, blind=FALSE)')

            X_normed = summarizedExperiment.assay(vst)

            self.sizefactor = np.asarray(ro.r('dds$sizeFactor'))
        else:
            assert self.dispersion is not None, "Please run clustering on train set first"
            ro.r('dds <- DESeq2::DESeq(dds, quiet=TRUE)')
            # Reuse dispersion estimation derived from train set
            ro.r.assign('dispersion', self.dispersion)
            ro.r('DESeq2::dispersionFunction(dds) <- dispersion')
            # Apply Variance stabilisation
            vst = ro.r('DESeq2::varianceStabilizingTransformation(dds, blind=FALSE)')
            X_normed = summarizedExperiment.assay(vst)

        X_normed_batchcorrected = self.apply_batch_correction(X_normed=X_normed, coldata=coldata, batch_keys=batch_keys)

        return X_normed_batchcorrected.T

    def deseq2_normalise(self, X: np.ndarray, parts, coldata: pd.DataFrame, batch_keys: [str, list]):
        """ Normalise data using DESeq2

        Parameters
        ----------
        X : Filtered, raw count matrix

        Returns
        -------

        """
        # 1. Create deseq2 object
        # 1.1 First read out save sample ids to dataframe
        # coldata = self._parts.sample_names()
        sample_id = pd.DataFrame({'sample': coldata})
        # 1.2 Use design 1 to not make use of any design information
        dds = deseq2.DESeqDataSetFromMatrix(
            countData=X.T, colData=self.pandasToR(df=sample_id), design=Formula("~1"))

        if parts._model_status != 'test':
            # Calculate sizeFactors on whole train or application set
            # remove sizeFactors on train set, otherwise dispersion will be calculated using
            # the existing sizeFactors which were calculate on whole dataset
            # dds.sizeFactor = None
            dds = deseq2.estimateSizeFactors_DESeqDataSet(dds, type='ratio', quiet=True)

            #  Calculate dispersion using precalculated sizefactors
            if self.dispersion is None:
                # Calculate dispersion only once
                self.dispersion = deseq2.dispersionFunction(dds)
                dds = deseq2.DESeq(dds, quiet=True)
            else:
                dds.dispersionFunction = self.dispersion

            # Apply Variance stabilisation
            vst = deseq2.varianceStabilizingTransformation(dds, blind=False)

            X_normed = summarizedExperiment.assay(vst)

            dollar = base.__dict__["$"]
            self.sizefactor = np.asarray(dollar(dds, "sizeFactor"))
        else:
            assert self.dispersion is not None, "Please run clustering on train set first"
            dds = deseq2.DESeq(dds, quiet=True)
            dds.dispersionFunction = self.dispersion
            vst = deseq2.varianceStabilizingTransformation(dds, blind=False)
            X_normed = summarizedExperiment.assay(vst)

        X_normed_batchcorrected = self.apply_batch_correction(X_normed=X_normed, coldata=coldata, batch_keys=batch_keys)

        return X_normed_batchcorrected.T

    def edger_normalise(self, X: np.ndarray, parts, coldata: pd.DataFrame, batch_keys: [str, list]):
        # 1. Create DGEList object
        # 1.1 First read out save sample ids to dataframe
        # coldata = self._parts.sample_names()
        sample_id = pd.DataFrame({'sample': coldata.index})
        # 1.2 Use design 1 to not make use of any design information
        X = X.transpose()
        y = edger.DGEList(counts=X, samples=sample_id)
        if parts._model_status != 'test':
            # Library Size Normalization: Calculate normalization factors on training set
            y = edger.calcNormFactors(y, method="TMM")

            # Estimate dispersion - robust against outliers
            ro.r.assign('coldata', sample_id)
            design = ro.r('model.matrix(~1, data=coldata)')
            y = edger.estimateDisp(y, design, robust=True)
            self.dispersion = y

            # Get normalized counts of train set
            X_normed = edger.cpm(y, normalized_lib_sizes=True, log=True)
        else:
            assert self.dispersion is not None, "Please run clustering on train set first"
            # Get normalized counts of test set
            dollar = base.__dict__["$"]
            X_normed = edger.cpm(
                y, normalized_lib_sizes=True, dispersion=np.asarray(dollar(self.dispersion, "tagwise.dispersion")),
                log=True)

        X_normed_batchcorrected = self.apply_batch_correction(X_normed=X_normed, coldata=coldata, batch_keys=batch_keys)

        return X_normed_batchcorrected.T

    def minmax_normalise(self, X: np.ndarray, parts, coldata: pd.DataFrame, batch_keys: [str, list]):
        if parts._model_status != 'test':
            # Min Max scaler for normalisation
            self.norm_factor = MinMaxScaler()

            # Apply transformation
            X_normed = self.norm_factor.fit_transform(X=X)
        else:
            X_normed = self.norm_factor.transform(X=X)

        # Apply batch corrrection
        X_normed_batchcorrected = self.apply_batch_correction(
            X_normed=X_normed.T, coldata=coldata, batch_keys=batch_keys)

        return X_normed_batchcorrected.T

    def pandasToR(self, df: pd.DataFrame):
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

    def apply_batch_correction(self, X_normed, coldata, batch_keys=None):
        # Get col data stored in df
        # coldata = self._parts.observations()
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
