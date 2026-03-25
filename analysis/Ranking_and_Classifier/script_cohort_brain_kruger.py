import os
import pickle
import numpy as np
import pandas as pd
import scanpy
from easydict import EasyDict as edict
from tqdm import tqdm

from config import Config
from dataframe_tools import partition

from modeling import utils
from modeling import FeatureModels
from sklearn.linear_model import LinearRegression
from sklearn.svm import SVR, SVC
from sklearn.metrics import confusion_matrix

import matplotlib.pyplot as pp
pp.ion()
import seaborn as sn

endotypes = ['E1','E2','E3','E4','E5','E6','E7','E8','E9','E10','E11','E12','E13']

def eval_function( x,y ):
    return np.mean( (x-y) ** 2, axis=1 )

class ExportModule :
    def __init__( self, mode='il23' ):
        args = edict()
        args.rank_model = 'LinearRegression'
        args.feat_model = None 

        self.mode = mode

        self.cfg = Config( args )

        adata = scanpy.read( self.cfg.data_path )

        var = adata.var 
        ensamble_gene_id = {}
        valid_genes = []

        for e,n in zip( var.gene_name, var.external_gene_name ):
            valid_genes.append(n)
            ensamble_gene_id[e] = n

        self.ensamble_gene_id = ensamble_gene_id 
        self.valid_genes = valid_genes

        self.gene_map = None

    def load_brain( self ):
        adata = scanpy.read( self.cfg.data_path )

        df = adata.to_df() 
        df = df.T         

        self.brain = {}
        self.brain['expressions'] = df

    def load_kruger( self ):
        expressions = pd.read_csv(f'./data-cohort/kruger_{self.mode}/expression_data.csv', index_col=0)

        self.kruger = {}
        self.kruger['expressions'] = expressions

    def create_gene_map( self, brain_data, kruger_data ):
        all_genes = list(brain_data.index)
        self.gene_map = {}
        for gene in tqdm(all_genes) :
            self.gene_map[gene] = self.ext_genes_and_ids( kruger_data, gene )

     
    def ext_genes_and_ids( self, data, key ):
        genes = list(data.index)

        selected_genes = []
        for gene in genes :

            if '...' in gene :
                continue

            if key in gene :
                selected_genes.append(gene)

        return selected_genes


    def common_genes( self, data ):

        genes = None 
        for d in data :
            index = set(d.index)

            if genes is None :
                genes = index 
            else :
                genes = genes.intersection(index)

        return list(genes) 

    def get_data( self, data, genes, gene_map ):
        
        expressions = [] 

        valid_genes = []

        for gene in genes:
            gg = gene_map[gene] 

            if len(gg) == 0 :
                continue

            d = data.loc[ gg ]
            d = d.to_numpy()
            d = np.mean(d,axis=0)

            valid_genes.append( gene )
            expressions.append( d )

        expressions = np.array( expressions )
        out = pd.DataFrame( expressions, index=valid_genes, columns=data.columns )

        return out
        index = self.s

        expressions = []

        for gene in genes :
            d = self.expression_data.loc[ gene_map[gene] ]
            d = d.to_numpy()
            d = np.mean(d,axis=0)

            expressions.append(d)

        expressions = np.array( expressions ).T

        data = pd.DataFrame( expressions, index=index, columns=genes )

        ids = np.array( list(self.meta_data.index) )

        labels = []
        for id in ids :
            t = self.meta.loc[id,'title']
            r = self.meta.loc[id,'pasi75resp']

            if 'NL Week 0' in t and r == 'Yes' :
                labels.append( annot_count + 0)
            elif 'LS Week 0' in t and r == 'Yes' :
                labels.append( annot_count + 1)
            elif 'LS Week 12' in t and r == 'Yes' :
                labels.append( annot_count + 2)
            elif 'NL Week 0' in t and r == 'No' :
                labels.append( annot_count + 3)
            elif 'LS Week 0' in t and r == 'No' :
                labels.append( annot_count + 4)
            elif 'LS Week 12' in t and r == 'No' :
                labels.append( annot_count + 5)
            else :
                print('Ooops')

        annots = np.array( labels ).reshape([-1,1])
        annots = pd.DataFrame( annots, index=index )
        return data, annots.to_numpy().ravel()



    def merge( self ):
        brain_exp = self.brain['expressions']
        kruger_exp = self.kruger['expressions']

        if self.gene_map is None :
            self.create_gene_map( brain_exp, kruger_exp )

        kruger_exp = self.get_data( kruger_exp, list(brain_exp.index), self.gene_map  )

        brain_exp.sort_index(inplace=True)

        brain_meta = pd.DataFrame( np.ones( len(brain_exp.columns), dtype=int ), 
                                   index=brain_exp.columns, columns=['batch'] )

        kruger_exp.sort_index(inplace=True)

        kruger_meta = pd.DataFrame( np.ones( len(kruger_exp.columns), dtype=int )*2, 
                                      index=kruger_exp.columns, columns=['batch'] )

        genes = self.common_genes( [ brain_exp, kruger_exp ] )

        brain_exp = brain_exp.loc[genes]
        kruger_exp = kruger_exp.loc[genes]

        counts = pd.concat( [ brain_exp, kruger_exp ], axis=1 ) 
        meta = pd.concat( [ brain_meta, kruger_meta ], axis=0 )

        counts.to_csv(f'./data-cohort/brain_kruger_{self.mode}/counts.csv')
        meta.to_csv(f'./data-cohort/brain_kruger_{self.mode}/meta.csv')

class AnalyticsModule :
    def __init__( self, mode ):
        args = edict()
        args.rank_model = 'LinearRegression'
        args.feat_model = None 

        self.cfg = Config( args ) 

        self.mode = mode

    def load_data( self ):

        adata = scanpy.read( self.cfg.data_path )

        index = adata.to_df().index
        pattern_ = adata.obs["Pattern"].to_numpy()
        endotypes_ = adata.obs["Endotypes"].to_numpy()

        aa = pd.DataFrame( index=index )
        aa['pattern'] = pattern_ 
        aa['endotypes'] = endotypes_

        print( f'./data-cohort/brain_kruger_{self.mode}/counts_log2_DESeq2_corrected.csv' )

        meta = pd.read_csv(f'./data-cohort/brain_kruger_{self.mode}/meta.csv', index_col=0)
        data = pd.read_csv(f'./data-cohort/brain_kruger_{self.mode}/counts_log2_DESeq2_corrected.csv', index_col=0)

        brain_indices = np.array(list(meta[ meta['batch'] == 1 ].index))
        kruger_indices = np.array(list(meta[ meta['batch'] == 2 ].index))

        data = data.T 

        brain_data = data.loc[ brain_indices ]
        brain_annots = aa
        kruger_data = data.loc[ kruger_indices ]

        kruger_annots = pd.read_csv(f'data-cohort/kruger_{self.mode}/metadata.csv', index_col=0)

        title = kruger_annots['title'].copy()

        arr = []

        for index in title.index :
            v = title.loc[index]
            v = v[11:]

            title.loc[index] = v

        kruger_annots['check'] = title

        keep = kruger_annots[ kruger_annots['check'] == 'LS Week 0' ].index 

        print( len(keep) )

        kruger_data = kruger_data.loc[keep]
        kruger_annots = kruger_annots.loc[keep]

        #print( kruger_annots )
        #michigan_annots = pd.read_excel('data-cohort/michigan/Michigan_Data_Original_Summarized.xlsx', index_col=0)

        self.data = {}
        self.data['brain_data'] = brain_data
        self.data['brain_annots'] = brain_annots 
        self.data['kruger_data'] = kruger_data 
        self.data['kruger_annots'] = kruger_annots

        self.parts = partition( brain_annots, "endotypes", endotypes, 
                                self.cfg.parts_ratio, self.cfg.parts_num,
                                self.cfg.sampling_ratio, self.cfg.sampling_num,
                                numpy_seed=self.cfg.seed )


    def create_gene_scores( self, recompute=False ):

        cache_file = f'./data-cohort/brain_kruger_{self.mode}/gene_scores.pkl'

        if os.path.exists( cache_file ) and not recompute:
            with open( cache_file, 'rb' ) as ff :
                self.gene_scores = pickle.load(ff)
                return

        data = self.data['brain_data']
        brain_annots = self.data['brain_annots']
        classifier = 'LinearRegression'

        labels = []

        for endotype in endotypes :
            annots = brain_annots['endotypes']

            annots = (annots == endotype).astype(float).to_numpy()
            annots = annots * 2 - 1

            annots = annots.reshape([-1,1])
            labels.append( annots )

        labels = np.concatenate( labels, axis=1 )
        labels = pd.DataFrame( labels, index=brain_annots.index, columns=endotypes )

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
        gene_scores = np.mean( cross_val, axis=0 )

        with open( cache_file, 'wb' ) as ff :
            pickle.dump( gene_scores, ff )

        self.gene_scores = gene_scores

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

        genes = list(self.data['brain_data'].columns)
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

        print( genes )

        data = self.data['brain_data'].loc[:,genes]
        annots = self.data['brain_annots'].loc[:,'endotypes']

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

        return np.mean( diag ), len(genes)

    def michigan_meta_data( self ):
        annots = pd.read_csv('data-cohort/michigan/sampleinfo_0.12.txt', delimiter="\t", index_col=0 )

        #annots.to_excel('data-cohort/michigan/sampleinfo_0.12.xlsx')

        print( annots.shape )
        print( self.data['michigan_data'].shape )

        return

        annots = pd.read_excel('data-cohort/kiel/Pso_patients_clin_data_short_20240311.mapping.blinded.xlsx')
        kiel_data = self.data['kiel_data']

        indices = list(kiel_data.index)

        for idx in indices :
            sub_annots = annots[ annots['KEYERICH_SAMPLE_ID'] == idx ]
            print( sub_annots['pasi_score_result'] )

 
        #print( data['KEYERICH_SAMPLE_ID'] )


    def do_stuff2( self, num_genes ):
        genes = self.select_genes(num_genes)

        data = self.data['brain_data'].loc[:,genes]
        annots = self.data['brain_annots'].loc[:,'endotypes']

        data = data.to_numpy()
        annots = annots.to_numpy()

        l_annots = self.category_labels( annots )
        m = SVC( kernel='linear' )
        m.fit( data, l_annots )


        kruger_data = self.data['kruger_data'].loc[:,genes]
        kruger_annots = self.data['kruger_annots']

        aa = []

        index = kruger_annots[ kruger_annots['pasi75resp'] == 'Yes' ].index
        dd = kruger_data.loc[ index ]

        pred = m.predict( dd.to_numpy() )

        a = np.zeros( len(endotypes) )
        for p in pred :
            a[p] += 1

        aa.append( a )

        index = kruger_annots[ kruger_annots['pasi75resp'] == 'No' ].index
        dd = kruger_data.loc[ index ]

        a = np.zeros( len(endotypes) )
        if len( index ) > 0 :
            pred = m.predict( dd.to_numpy() )
            for p in pred :
                a[p] += 1

        aa.append( a )


        aa = np.array(aa, dtype=int )

        aa = pd.DataFrame( aa, index=[f'{self.mode.upper()} R', f'{self.mode.upper()} NR'], columns=endotypes )

        print( aa )

        aa.to_excel(f'data-cohort/brain_kruger_{self.mode}/results_corrected.xlsx')

    def export_endotypes( self, num_genes ):
        genes = self.select_genes(num_genes)

        data = self.data['brain_data'].loc[:,genes]
        annots = self.data['brain_annots'].loc[:,'endotypes']

        data = data.to_numpy()
        annots = annots.to_numpy()

        l_annots = self.category_labels( annots )
        m = SVC( kernel='linear' )
        m.fit( data, l_annots )

        kruger_data = self.data['kruger_data'].loc[:,genes]

        pred = m.predict( kruger_data.to_numpy() )

        pred_endotypes = []

        for p in pred :
            pred_endotypes.append( endotypes[p] )

        pred_endotypes = np.array( pred_endotypes ).reshape([-1,1])

        pred_endotypes = pd.DataFrame( pred_endotypes, index=kruger_data.index, columns=['Endotype'])
        
        pred_endotypes.to_excel(f'./data-cohort/brain_kruger_{self.mode}/kruger_{self.mode.upper()}_endotype_assignment.xlsx')


