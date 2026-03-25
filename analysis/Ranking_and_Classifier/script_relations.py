import numpy as np 
import pickle
import pandas as pd
import scanpy
from easydict import EasyDict as edict

from matplotlib import pyplot as pp 
pp.ion()

from config import Config


endotypes = ['E1','E2','E3','E4','E5','E6','E7','E8','E9','E10','E11','E12','E13']

groups = [
    [['E5'],['E6']],
    [['E9'],['E10']],
    [['E12'],['E13']],
    [['E1'],['E2']],
    [['E3'],['E4']],
    [['E8'],['E9','E10']],
    [['E7'],['E8','E9','E10']],
    [['E11'],['E12','E13']],
    [['E1'],['E2']],
    [['E3'],['E4']],
    [['E5','E6'],['E7','E8','E9','E10']],
    [['E1','E2'],['E3','E4']],
    [['E11','E12','E13'],['E1','E2','E3','E4']],
    [['E5','E6','E7','E8','E9','E10'],['E11','E12','E13','E1','E2','E3','E4']]
]

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

        with open('data/gene_scores.pkl','rb') as ff :
            self.data = pickle.load(ff)

        self.brain_data, self.brain_annots = load_brain(self.cfg)

    def _to_index( self, values ):
        indices = []

        for v in values :
            indices.append( endotypes.index(v) )

        return np.array(indices)

    def _get_descriptor( self, values, indices ):
        vv = values[indices].mean(axis=0)
        return vv

    def _build_name( self, group ):

        n0 = '_'.join(group[0])
        n1 = '_'.join(group[1])

        return '%s_vs_%s' % (n0,n1)

    def _get_samples( self, endotypes, annots ):

        ll = []

        for e in endotypes :
            aa = annots[ annots == e ]

            for idx in list(aa.index) :
                ll.append(idx)

        return np.array(ll)

    def get_expressions( self, g ):
        indices = self._get_samples( g, self.brain_annots['endotypes'] )
        values = self.brain_data.loc[indices]

        return values

    def gene_map( self, g0, g1, num_genes=50, mode="similar" ):

        inds0 = self._to_index(g0)
        inds1 = self._to_index(g1)

        genes = []
        scores = []


        for gene, values in self.data.items() :
            v0 = self._get_descriptor( values, inds0 )
            v1 = self._get_descriptor( values, inds1 )

            if mode == "disimilar" :
                v1 = 1 - v1

            diff = np.abs(v0 - v1)

            genes.append( gene )
            scores.append( diff )

        genes = np.array( genes )
        scores = np.array( scores )

        pp.figure()
        pp.matshow( np.abs(scores), aspect="auto" ) 
        pp.title("Descriptor difference between the two groups for each gene")

        scores = np.mean( scores, axis=1 )
        sinds = np.argsort( scores )

        pp.figure()

        g0_expressions = self.get_expressions( g0 )
        g1_expressions = self.get_expressions( g1 ) 

        for gene in genes[sinds[:num_genes]] :
            g0_values = g0_expressions[gene].to_numpy()
            g1_values = g1_expressions[gene].to_numpy()

            print( gene, np.mean( g0_values ), np.std( g0_values ), np.mean( g1_values ), np.std( g1_values ) )


        #scores = np.linalg.norm( scores,axis=1 )
        #sinds = np.argsort( scores )

    def similar_genes( self, g0, g1 ):

        inds0 = self._to_index(g0)
        inds1 = self._to_index(g1)

        genes = []
        scores = []


        for gene, values in self.data.items() :
            v0 = self._get_descriptor( values, inds0 )
            v1 = self._get_descriptor( values, inds1 )

            diff = np.linalg.norm( v0 - v1 )

            genes.append( gene )
            scores.append( diff )

        genes = np.array( genes )
        scores = np.array( scores )

        inds = np.argsort( scores )

        return genes[inds[:20]]

    def disimilar_genes( self, g0, g1 ):
        inds0 = self._to_index(g0)
        inds1 = self._to_index(g1)

        genes = []
        scores = []


        for gene, values in self.data.items() :
            v0 = self._get_descriptor( values, inds0 )
            v1 = 1-self._get_descriptor( values, inds1 )

            diff = np.linalg.norm( v0 - v1 )

            genes.append( gene )
            scores.append( diff )

        genes = np.array( genes )
        scores = np.array( scores )

        inds = np.argsort( scores )

        return genes[inds[:20]]

    def create_gene_pool( self ):
        
        gene_pool = set()

        for group in groups :
            gg = self.similar_genes( group[0], group[1] )

            for g in gg:
                gene_pool.add(g)

            gg = self.dismilar_genes( group[0], group[1] )

            for g in gg :
                gene_pool.add(g)

        self.gene_pool = np.array(list(gene_pool))
        

    def do_stuff( self ):

        columns = []

        values = []


        for group in groups :

            columns.append( self._build_name(group) )

            inds0 = self._to_index(group[0])
            inds1 = self._to_index(group[1])

            mat0 = []
            mat1 = []


            for gene in self.gene_pool :

                v0 = self._get_descriptor( self.data[gene], inds0 )
                v1 = self._get_descriptor( self.data[gene], inds1 )

                mat0.append(v0)
                mat1.append(v1)

            mat0 = np.array( mat0 )
            mat1 = np.array( mat1 )

            diff = np.abs( mat0 - mat1 ).mean(axis=1)

            diff = diff.reshape([-1,1])

            values.append( diff )

        values = np.concatenate( values, axis=1 )

        df = pd.DataFrame( values, index=self.gene_pool, columns=columns )
        df.to_csv('data.csv')

    def do_stuff2( self ):


        gene = 'RGS6'


        d = 1-self.data[gene]

        pp.matshow(d)

    def do_stuff3( self ):

        def get_samples( endotypes, annots ):

            ll = []

            for e in endotypes :
                aa = annots[ annots == e ]

                for idx in list(aa.index) :
                    ll.append(idx)

            return np.array(ll)

        g0, g1 = groups[12]

        g0_indices = get_samples( g0, self.brain_annots['endotypes'] )
        g1_indices = get_samples( g1, self.brain_annots['endotypes'] )

        g0_values = self.brain_data.loc[ g0_indices ]
        g1_values = self.brain_data.loc[ g1_indices ]

        print( g0_values.shape )
        print( g1_values.shape )

