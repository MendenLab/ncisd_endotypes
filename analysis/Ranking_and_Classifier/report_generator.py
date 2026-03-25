import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats import multitest
from matplotlib import pyplot as pp

class ReportGenerator :
    def __init__( self, data, annot, feats, scores, patterns ):
        self.data = data.T
        self.annot = annot
        self.feats = feats
        self.scores = scores
        self.patterns = patterns

    def calculate_pvalues( self, feats, X, Y ):
        interest_patterns = self.patterns
        pattern_indices = []

        for p in interest_patterns :
            pinds = np.where( Y == p )[0]
            ninds = np.where( Y != p )[0]
            pattern_indices.append([ pinds, ninds ])

        p_values = []

        print( X.shape )

        for idx, [ pinds, ninds ] in enumerate( pattern_indices ):
            pX = X[:,feats][pinds]
            nX = X[:,feats][ninds]

            print( pX.shape, nX.shape )

            U1, pvals = mannwhitneyu( pX, nX )
            _ ,pvals = multitest.fdrcorrection( pvals )

            p_values.append( pvals )

        return p_values

    def build_gene_list( self, fname, scores, feats, X, Y, pattern_set, p_values ):
        interest_patterns = self.patterns
        pattern_indices = []

        for p in interest_patterns :
            idx = pattern_set.index(p)
            pinds = np.where( Y == idx )[0]
            ninds = np.where( Y != idx )[0]
            pattern_indices.append([ pinds, ninds ])

        headers = []
        headers.append("Scores")
        headers.append("Gene Name")
        headers.append("External Gene Name")

        for p in interest_patterns :
            headers.append(f'Pattern {p} (Mean)')
            headers.append(f'Pattern {p} (Std)')
            headers.append(f'Pattern {p} (p-value)')

        with open( fname, 'w' ) as ff :

            ff.write( ','.join(headers) + '\n' )

            for idx, (s,f) in enumerate(zip( scores, feats )):
                gene = self.data.gene_names()[f]

                row = []
                row.append( str(s) )
                row.append( gene['gene_name'] )
                row.append( gene['external_gene_name'] )

                vals = X[:,f].ravel()

                for idy, (pinds, ninds) in enumerate(pattern_indices) :
                    row.append( str( np.mean(vals[pinds]) ) )
                    row.append( str( np.std(vals[pinds]) ) )
                    row.append(str(p_values[idy][idx]))

                ff.write( ','.join(row) + '\n' )


        #with open( fname, 'w' ) as ff :

    def build_pd( self, scores, feats, X, Y, p_values, gene_names ):
        interest_patterns = self.patterns
        pattern_indices = []

        for p in interest_patterns :
            pinds = np.where( Y == p )[0]
            ninds = np.where( Y != p )[0]
            pattern_indices.append([ pinds, ninds ])

        headers = []
        headers.append("Scores")
        headers.append("Gene Name")
        headers.append("External Gene Name")

        for p in interest_patterns :
            headers.append(f'Pattern {p} (Mean)')
            headers.append(f'Pattern {p} (Std)')
            headers.append(f'Pattern {p} (p-value)')

        mat = []
        for idx, (s,f) in enumerate(zip( scores, feats )):
            gene = gene_names[f]

            row = []
            row.append( str(s) )
            row.append( gene )
            row.append( gene )
            vals = X[:,f].ravel()

            for idy, (pinds, ninds) in enumerate(pattern_indices) :
                row.append( str( np.mean(vals[pinds]) ) )
                row.append( str( np.std(vals[pinds]) ) )
                row.append(str(p_values[idy][idx]))

            mat.append( row )


        return pd.DataFrame( mat, columns=headers )

    def __call__( self, fname=None ):
        X, Y, pattern_set = self.data.get_norm()

        p_values = self.calculate_pvalues( self.feats, X, Y, pattern_set )
        #self.save_gene_list( fname, self.scores, self.feats, X, Y, pattern_set, p_values)

    def plot( self, title=None, fname=None ):
        
        X = self.data.to_numpy().T
        Y = self.annot.to_numpy().ravel()

        gene_names = self.data.index.to_numpy()
        p_values = self.calculate_pvalues( self.feats, X, Y )
        data = self.build_pd( self.scores, self.feats, X, Y, p_values, gene_names )

        means = []
        stds = []

        for p in self.patterns :
            means.append( data[[f'Pattern {p} (Mean)']].to_numpy().ravel().astype(np.float32) )
            stds.append( data[[f'Pattern {p} (Std)']].to_numpy().ravel().astype(np.float32) )

        pp.close('all')

        fig, ax = pp.subplots()

        inds = np.arange( len(means[0]))

        colors = pp.cm.get_cmap(name="prism",lut=len(self.patterns))

        plots = []

        for idx, p in enumerate(self.patterns):
            pt = ax.errorbar( inds, means[idx], yerr=stds[idx],
                              label=self.patterns[idx], fmt='.', color=colors(idx) )
            plots.append(pt)

        ax.legend( loc='upper right', ncol=3 )
        ax.set_xticks(inds)
        ax.set_xticklabels(gene_names[self.feats],rotation=90,fontsize=5)
        ax.grid()

        if title != None :
            pp.title( title )
        if fname != None :
            pp.savefig( fname )
        else :
            pp.show()

class ReportGeneratorHeirarchy :

    def _group_indices( self, annot, groups ):

        positive_indices = {}
        negative_indices = []

        index = list(annot.index)
        for i, a in  zip( index, annot ):

            if a in groups['positives'] :
                if not a in positive_indices :
                    positive_indices[a] = []

                positive_indices[a].append(i)
            else :
                negative_indices.append(i)

        return positive_indices, negative_indices


    def __init__( self, data, annot, feats, scores, groups ):
        self.data = data
        self.annot = annot
        self.feats = feats
        self.scores = scores
        self.groups = groups

    def plot( self ):

        positives, negatives = self._group_indices( self.annot, self.groups )
        gene_names = self.data.columns.to_numpy()[ self.feats ]

        means = []
        stds = []

        data = self.data.loc[:,gene_names]
        
        for name, indices in positives.items() :
            dd = data.loc[indices,:].to_numpy()

            mean = dd.mean( axis=0 )
            std = dd.std( axis=0 )

            means.append( mean )
            stds.append( std )

        dd = data.loc[negatives,:].to_numpy()

        mean = dd.mean( axis=0 )
        std = dd.std( axis=0 )

        means.append( mean )
        stds.append( std )

        names = list( positives.keys() ) + ['Others']

        print( names )

        fig, ax = pp.subplots()

        for idx in range(3):
            pt = ax.errorbar( gene_names, means[idx], yerr=stds[idx],
                              label=names[idx], fmt='.' )

        ax.legend( loc='upper right', ncol=3 )
        ax.set_xticks(gene_names)
        ax.set_xticklabels(gene_names,rotation=90,fontsize=5)
        ax.grid()

        pp.show()

        #pp.savefig('plot.pdf')



    
        

        #for gene in gene_names :
        #    data = self.data.loc[:,gene]

        #    for name, indices in positives.items() :
        #        dd = data[indices].to_numpy()

        #    dd = data[negatives].to_numpy()    

        #    print( dd.shape )
        

