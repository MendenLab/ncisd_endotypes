import os
import numpy as np
import pandas as pd
from easydict import EasyDict as edict
from tqdm import tqdm
import pickle

from config import Config
from tools import load_brain
from dataframe_tools import partition

from modeling import utils
from modeling import FeatureModels
from sklearn.linear_model import LinearRegression
from sklearn.svm import SVR, SVC
from sklearn.metrics import confusion_matrix


import matplotlib.pyplot as pp
pp.ion()
import seaborn as sn

def eval_function( x,y ):
    return np.mean( (x-y) ** 2, axis=1 )

endotypes = ['E1','E2','E3','E4','E5','E6','E7','E8','E9','E10','E11','E12','E13']
#endotypes = ['E1','E2','E3','E4','E5','E6','E7','E9','E10','E11','E12','E13']


class Module :
    def __init__( self ):
        args = edict()
        args.rank_model = 'LinearRegression'
        args.feat_model = None 

        self.cfg = Config(args)
        self.brain_data, self.brain_annots = load_brain( self.cfg )

        self.parts = partition( self.brain_annots, "endotypes", endotypes, 
                                self.cfg.parts_ratio, self.cfg.parts_num,
                                self.cfg.sampling_ratio, self.cfg.sampling_num,
                                numpy_seed=self.cfg.seed )

    def create_gene_scores( self ):
        if os.path.exists('data/gene_scores.pkl') :
            with open('data/gene_scores.pkl','rb') as ff :
                self.gene_scores = pickle.load(ff)
            return


        data = self.brain_data
        classifier = 'LinearRegression'

        labels = []

        for endotype in endotypes :
            annots = self.brain_annots['endotypes']

            annots = (annots == endotype).astype(float).to_numpy()
            annots = annots * 2 - 1

            annots = annots.reshape([-1,1])
            labels.append( annots )

        labels = np.concatenate( labels, axis=1 )
        labels = pd.DataFrame( labels, index=self.brain_annots.index, columns=endotypes )

        data = self.brain_data# 

        cross_val = []

        for part in tqdm(self.parts) :
            X_tr = data.loc[ part['train'] ]
            y_tr = labels.loc[ part['train'] ]
            X_te = data.loc[ part['test'] ]
            y_te = labels.loc[ part['test'] ]

            model = FeatureModels( LinearRegression )
            model.fit( X_tr, y_tr )

            val = model.eval(X_te).astype( np.float32 )

            y_te = np.tile( y_te, ( len(val), 1, 1))

            scores = eval_function( val, y_te )
            cross_val.append( scores )

        cross_val = np.array( cross_val )
        self.gene_scores = np.mean( cross_val, axis=0 )

        with open('data/gene_scores.pkl','wb') as ff :
            pickle.dump( self.gene_scores, ff )

    def create_gene_scores_endotype( self ):
        if os.path.exists('data/gene_scores_endotype.pkl') :
            with open('data/gene_scores_endotype.pkl','rb') as ff :
                self.gene_scores_endotype = pickle.load(ff)
            return

        data = self.brain_data 
        annots = self.brain_annots['endotypes'].to_frame(name="endotypes")

        gene_scores_endotype = {}

        for idx, e0 in enumerate(endotypes):
            for idy, e1 in enumerate(endotypes[:idx]):

                cross_val = []

                print(e0, e1)

                for part in tqdm(self.parts) :
                    X_tr = data.loc[ part['train'] ]
                    y_tr = annots.loc[ part['train'] ]
                    X_te = data.loc[ part['test'] ]
                    y_te = annots.loc[ part['test'] ]

                    tr_index = y_tr[ (y_tr['endotypes'] == e0) | (y_tr['endotypes'] == e1) ].index  
                    te_index = y_te[ (y_te['endotypes'] == e0) | (y_te['endotypes'] == e1) ].index

                    X_tr = X_tr.loc[ tr_index ]
                    y_tr = y_tr.loc[ tr_index ]
                    X_te = X_te.loc[ te_index ]
                    y_te = y_te.loc[ te_index ]

                    ll_tr = (y_tr == e0).astype(float)
                    ll_te = (y_te == e1).astype(float)

                    ll_tr = ll_tr * 2 - 1 
                    ll_te = ll_te * 2 - 1

                    ll_tr = pd.DataFrame( ll_tr.to_numpy(), index=y_tr.index, columns=y_tr.columns )
                    ll_te = pd.DataFrame( ll_te.to_numpy(), index=y_te.index, columns=y_tr.columns )

                    model = FeatureModels( LinearRegression )
                    model.fit( X_tr, ll_tr )

                    ll_te = ll_te.to_numpy()

                    val = model.eval(X_te).astype( np.float32 )
                    ll_te = np.tile( ll_te, ( len(val), 1, 1))

                    scores = eval_function( val, ll_te )

                    cross_val.append( scores )

                cross_val = np.array( cross_val )
                cross_val = np.mean( cross_val, axis=0 )

                gene_scores_endotype[f'{e0}_{e1}'] = cross_val

        with open('data/gene_scores_endotype.pkl','wb') as ff :
            pickle.dump( gene_scores_endotype, ff )

        self.gene_scores_endotype = gene_scores_endotype

    def gene_scores_df( self ):
        genes = np.array(self.brain_data.columns) 
        scores = self.gene_scores

        df = pd.DataFrame( scores, index=genes, columns = endotypes )

        df.to_csv('gene_endotype_scores.csv')



    def export( self ):
        genes = np.array(self.brain_data.columns) 


        scores = self.gene_scores.T 

        ranked = {}

        for e, s in zip( endotypes, scores ):
            sinds = np.argsort(s)
            ranked[e] = genes[sinds]


        out = {}
        out['ranked'] = ranked 
        out['partitions'] = self.parts

        with open('data.pkl','wb') as ff :
            pickle.dump( out, ff )

        print( self.brain_annots.loc[self.parts[0]['test'],'endotypes'].to_numpy() )

    def stats( self ):
        aa = self.brain_annots['endotypes'].to_numpy()
        for e in endotypes :
            num = len( np.where( aa == e )[0] )
            print( e, num )

    def binary_labels( self, labels, endotype ):
        
        ll = np.zeros( len(labels) )

        pos_inds = np.where( labels == endotype )[0]
        neg_inds = np.where( labels != endotype )[0]

        ll[pos_inds] = 1.0 
        ll[neg_inds] = -1.0

        return ll

    def category_labels( self, labels ):
        out = []

        for l in labels :
            out.append( endotypes.index(l) )

        return np.array( out )


    def select_genes( self, num_per_endotype ):

        genes = list(self.brain_data.columns)
        gene_set = set()
        for scores in self.gene_scores.T :
            sinds = np.argsort( scores )

            for i in sinds[:num_per_endotype] :
                gene_set.add( genes[i] )

        selected_genes = np.array(list(gene_set))

        return selected_genes

    def map_features( self, X_tr, y_tr, X_te ):

        blabels = []
        for e in endotypes :
            blabels.append( self.binary_labels(y_tr,e) )

        blabels = np.array(blabels).T
        blabels = pd.DataFrame( blabels, index=y_tr.index, columns=endotypes )

        model = FeatureModels( LinearRegression )
        model.fit( X_tr, blabels )

        X_tr_mapped = model.eval( X_tr ).transpose(1,0,2).reshape(( len(X_tr),-1 ) )
        X_te_mapped = model.eval( X_te ).transpose(1,0,2).reshape(( len(X_te),-1 ) )
    
        return X_tr_mapped, X_te_mapped 

    def do_stuff_part( self, data, annots, part ):
        X_tr = data.loc[ part['train'] ]
        y_tr = annots.loc[ part['train'] ] 
        X_te = data.loc[ part['test'] ]
        y_te = annots.loc[ part['test'] ]

        map_features = False

        if map_features :
            X_tr, X_te = self.map_features( X_tr, y_tr, X_te )
        else :
            X_tr = X_tr.to_numpy()
            X_te = X_te.to_numpy()
        y_tr = y_tr.to_numpy()
        y_te = y_te.to_numpy();

        use_svc = True
        use_linear = False

        if use_svc :
            l_tr = self.category_labels( y_tr )
            m = SVC( kernel='linear' )
            m.fit( X_tr, l_tr )
            pred = m.predict(X_te)

        if use_linear : 
            scores = [] 

            for idx, endotype in enumerate( endotypes ):
                l_tr = self.binary_labels( y_tr, endotype )

                #m = LinearRegression()
                m = SVR()
                m.fit( X_tr, l_tr )
        
                s = m.predict( X_te )
                scores.append( s )

            scores = np.array( scores ).T
            pred = np.argmax( scores, axis=1 )

        cmat = np.zeros(( len(endotypes), len(endotypes)) )

        for l,p in zip( y_te, pred ):
            l = endotypes.index(l)

            cmat[l,p] += 1

        return cmat 

    def normalize_cmat( self, cmat ):
        for idx, r in enumerate(cmat):
            r = r / r.sum()
            cmat[idx] = r
        return cmat


    def do_stuff( self, num_genes ):
        genes = self.select_genes(num_genes)

        data = self.brain_data.loc[:,genes]
        annots = self.brain_annots.loc[:,'endotypes']

        print( data.shape )

        cmats = []
        #accs = []

        for part in tqdm(self.parts) :
            cmat = self.do_stuff_part( data, annots, part )
            cmats.append(cmat)

            #ncmat = cmat.copy()
            #ncmat = self.normalize_cmat( ncmat )
            #accs.append( np.diag(ncmat).mean() )

        cmats = np.array( cmats )
        cmats = np.mean( cmats, axis=0 )

        cmats = self.normalize_cmat( cmats )

        diag = np.diag(cmats)
        df_conf_mat = pd.DataFrame( cmats, index=endotypes, columns=endotypes )

        pp.figure(figsize = (10,10))
        sn.heatmap(df_conf_mat, annot=True)

        pp.title(f'Mean acc : {np.mean(diag)} - Num genes {len(genes)}')
        pp.savefig(f'conf_mat_{num_genes}.pdf')

        #print( accs, np.mean(accs) )
        #print( np.mean(accs) )

        return np.mean( diag ), len(genes)

    def do_stuff2( self ):
        num_genes = [1,5,10,20,50,100,200,300,400, 500]

        pp.close('all')

        scores = []
        real_num_genes = []
        for ng in num_genes :
            s,n = self.do_stuff( ng )
            scores.append( s )
            real_num_genes.append(n)

        pp.figure()
        pp.plot( real_num_genes, scores, '.-' )
        pp.grid()

        pp.xlabel('Number of genes')
        pp.ylabel('Accuracy')

        pp.savefig('num_genes.pdf')





        
if __name__=="__main__" :
    m = Module()
    m.create_gene_scores()
    m.create_gene_scores_endotype()
