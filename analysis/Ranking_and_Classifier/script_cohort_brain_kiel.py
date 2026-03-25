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
    def __init__( self ):
        args = edict()
        args.rank_model = 'LinearRegression'
        args.feat_model = None 

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

    def load_brain( self ):
        adata = scanpy.read( self.cfg.data_path )

        df = adata.to_df() 

        counts = pd.DataFrame( adata.layers['counts'], index=df.index, columns=df.columns )

        counts = counts.T

        self.brain = {}
        self.brain['counts'] = counts

    def load_kiel( self ):
        counts = pd.read_csv( './data-cohort/kiel/Pso_patients_clin_data_short_20240311.counts.csv', 
                          index_col='ENSEMBL_GENE_ID' )

        index = list(counts.index)

        keep = []
        name = []

        for idx in index :
            if idx in self.ensamble_gene_id :
                keep.append(idx)
                name.append( self.ensamble_gene_id[idx] )

        counts = counts.loc[ keep, : ]
        counts.index = name 

        self.kiel = {}
        self.kiel['counts'] = counts

    def load_michigan( self ):
        pass

    def common_genes( self, data ):

        genes = None 
        for d in data :
            index = set(d.index)

            if genes is None :
                genes = index 
            else :
                genes = genes.intersection(index)

        return list(genes) 
    def merge( self ):
        metadata = pd.read_csv('./data-cohort/brain_kiel/Combine_MetaData_BRAIN_Kiel.csv',
                               index_col=0 )

        brain_counts = self.brain['counts']
        brain_counts.sort_index(inplace=True)

        brain_meta = pd.DataFrame( np.ones( len(brain_counts.columns), dtype=int ), 
                                   index=brain_counts.columns, columns=['batch'] )

        kiel_counts = self.kiel['counts']
        kiel_counts.sort_index(inplace=True)

        kiel_meta = pd.DataFrame( np.ones( len(kiel_counts.columns), dtype=int )*2, 
                                  index=kiel_counts.columns, columns=['batch'] )

        genes = self.common_genes( [ brain_counts, kiel_counts ] )

        brain_counts = brain_counts.loc[genes]
        kiel_counts = kiel_counts.loc[genes]

        counts = pd.concat( [ brain_counts, kiel_counts ], axis=1 ) 
        meta = pd.concat( [ brain_meta, kiel_meta ], axis=0 )

        counts = counts.T 
        counts = counts.loc[ metadata.index ]
        counts = counts.T 

        meta = meta.loc[ metadata.index ]

        brain_metadata = metadata.loc[ brain_counts.columns ]
        kiel_metadata = metadata.loc[kiel_counts.columns]

        brain_counts.to_csv('./data-cohort/brain_kiel/counts_brain.csv')
        brain_metadata.to_csv('./data-cohort/brain_kiel/metadata_brain.csv')

        kiel_counts.to_csv('./data-cohort/brain_kiel/counts_kiel.csv')
        kiel_metadata.to_csv('./data-cohort/brain_kiel/metadata_kiel.csv')

        counts.to_csv('./data-cohort/brain_kiel/counts.csv')
        meta.to_csv('./data-cohort/brain_kiel/meta.csv')

class AnalyticsModule :
    def __init__( self ):
        args = edict()
        args.rank_model = 'LinearRegression'
        args.feat_model = None 

        self.cfg = Config( args ) 

    def load_data( self ):

        def load_data_deseq2():
            print('Loading DESeq2')
            data = pd.read_csv('./data-cohort/brain_kiel/counts_log2_DESeq2.csv', index_col=0)
            return data 

        def load_data_tmm():
            print('Loading TMM') 
            data = pd.read_csv('./data-cohort/brain_kiel/Normed_counts_BRAIN_Kiel_TMM.csv', index_col=0)
            data = data.T
            return data

        adata = scanpy.read( self.cfg.data_path )

        index = adata.to_df().index
        pattern_ = adata.obs["Pattern"].to_numpy()
        endotypes_ = adata.obs["Endotypes"].to_numpy()

        aa = pd.DataFrame( index=index )
        aa['pattern'] = pattern_ 
        aa['endotypes'] = endotypes_

        meta = pd.read_csv('./data-cohort/brain_kiel/meta.csv', index_col=0)

        #data = load_data_deseq2()
        #data = load_data_tmm()


        brain_indices = np.array(list(meta[ meta['batch'] == 1 ].index))
        kiel_indices = np.array(list(meta[ meta['batch'] == 2 ].index))

        #data = data.T 

        #brain_data = data.loc[ brain_indices ]
        brain_data = pd.read_csv('data-cohort/brain_kiel/counts_log2_DESeq2_brain.csv',index_col=0).T
        brain_annots = aa
        #kiel_data = data.loc[ kiel_indices ]
        kiel_data = pd.read_csv('data-cohort/brain_kiel/counts_log2_DESeq2_kiel.csv',index_col=0).T

        kiel_annots = pd.read_excel('data-cohort/kiel/2025-01-26_Summary_KielData_RED_MODIFIED_O.xlsx')
        kiel_annots.index = kiel_annots['KEYERICH_SAMPLE_ID']

        self.data = {}
        self.data['brain_data'] = brain_data
        self.data['brain_annots'] = brain_annots 
        self.data['kiel_data'] = kiel_data 
        self.data['kiel_annots'] = kiel_annots

        self.parts = partition( brain_annots, "endotypes", endotypes, 
                                self.cfg.parts_ratio, self.cfg.parts_num,
                                self.cfg.sampling_ratio, self.cfg.sampling_num,
                                numpy_seed=self.cfg.seed )

    def create_gene_scores( self, recompute=False ):
        cache_file = './data-cohort/brain_kiel/gene_scores.pkl'

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

        #with open( cache_file, 'wb' ) as ff :
        #    pickle.dump( gene_scores, ff )

        self.gene_scores = gene_scores

    def export_gene_list( self ):
        genes = np.array(self.data['brain_data'].columns)
        scores = self.gene_scores.T 

        sorted_gene_list = {}

        for e, s in zip( endotypes, scores ):
            sinds = np.argsort(s)
            genes_sorted = genes[sinds]
            print( genes_sorted )
            sorted_gene_list[e] = genes_sorted

        with open('./data-cohort/brain_kiel/brain_kiel_sorted_genes.pkl','wb') as ff :
            pickle.dump( sorted_gene_list, ff )
        
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

        map_features = True

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
            m = SVC( kernel='rbf', C=1.5, class_weight='balanced' )
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

    def do_stuff2( self, num_genes ):
        genes = self.select_genes(num_genes)

        data = self.data['brain_data'].loc[:,genes]
        annots = self.data['brain_annots'].loc[:,'endotypes']

        data = data.to_numpy()
        annots = annots.to_numpy()

        l_annots = self.category_labels( annots )

        preds = pd.read_excel('./data-cohort/kiel/kiel_endotypes.xlsx', index_col=0)

        #m = SVC( kernel='linear' )
        #m.fit( data, l_annots )

        kiel_data = self.data['kiel_data'].loc[:,genes]
        kiel_annots = self.data['kiel_annots']

        def get_hist( key ):
            index = kiel_annots[ key ].index 

            a = np.zeros( len(endotypes) )

            if len( index ) == 0 :
                return a 

            pp = preds.loc[index].to_numpy()

            for p in pp :
                idx = endotypes.index(p)
                a[idx] += 1




            #dd = kiel_data.loc[index]
            
            #pred = m.predict( dd.to_numpy() )
            #for p in pred :
            #    a[p] += 1
            return a

        
        aa = []
        rows = []


        # IL23
        aa.append( get_hist( kiel_annots['IL23 SR M'] == 1 ) )
        rows.append('IL23 SR')
        aa.append( get_hist( kiel_annots['IL23 R M'] == 1 ) )
        rows.append('IL23 R')
        aa.append( get_hist( kiel_annots['IL23 NR M'] == 1 ) )
        rows.append('IL23 NR')
        aa.append( get_hist( kiel_annots['IL23 SNR M'] == 1 ) )
        rows.append('IL23 SNR')

        # TNF 
        aa.append( get_hist( kiel_annots['TNF SR M'] == 1 ) )
        rows.append('TNF SR')
        aa.append( get_hist( kiel_annots['TNF R M'] == 1 ) )
        rows.append('TNF R')
        aa.append( get_hist( kiel_annots['TNF R M'] == 1 ) )
        rows.append('TNF NR')
        aa.append( get_hist( kiel_annots['TNF SNR M'] == 1 ) )
        rows.append('TNF SNR')

        # IL-17
        aa.append( get_hist( kiel_annots['IL-17 SR'] == 1 ) )
        rows.append('IL-17 SR')
        aa.append( get_hist( kiel_annots['IL-17 R'] == 1 ) )
        rows.append('IL-17 R')
        aa.append( get_hist( kiel_annots['IL-17 R'] == 2 ) )
        rows.append('IL-17 NR')
        aa.append( get_hist( kiel_annots['IL-17 SNR'] == 1 ) )
        rows.append('IL-17 SNR')
 
        aa = np.array(aa, dtype=int )

        print( aa.sum() )

        aa = pd.DataFrame( aa, index=rows, columns=endotypes )

        print( aa )

        aa.to_excel('data-cohort/brain_kiel/results_corrected.xlsx')

    def export_endotypes( self, num_genes ):
        genes = self.select_genes(num_genes)

        data = self.data['brain_data'].loc[:,genes]
        annots = self.data['brain_annots'].loc[:,'endotypes']

        data = data.to_numpy()
        annots = annots.to_numpy()

        l_annots = self.category_labels( annots )
        m = SVC( kernel='linear' )
        m.fit( data, l_annots )

        kiel_data = self.data['kiel_data'].loc[:,genes]
        kiel_annots = self.data['kiel_annots']

        pred = m.predict( kiel_data.to_numpy() )

        print( pred )

        pred_endotypes = []

        for p in pred :
            pred_endotypes.append( endotypes[p] )

        pred_endotypes = np.array( pred_endotypes ).reshape([-1,1])

        pred_endotypes = pd.DataFrame( pred_endotypes, index=kiel_data.index, columns=['Endotype'])

        print( pred_endotypes )
        
        #pred_endotypes.to_excel('./data-cohort/brain_kiel/endotype_assignment.xlsx')

    def do_stuff3( self ):

        def eval( genes ):
            data = self.data['brain_data'].loc[:,genes]
            annots = self.data['brain_annots'].loc[:,'endotypes']

            data = data.to_numpy()
            annots = annots.to_numpy()

            l_annots = self.category_labels( annots )

            m = SVC( kernel='linear' )
            m.fit( data, l_annots )

            kiel_data = self.data['kiel_data'].loc[:,genes]
            kiel_annots = self.data['kiel_annots']

            return m.predict( kiel_data.to_numpy() )

        genes = self.select_genes(100)
        pred_base = eval(genes)


        num_genes = [1,2,3,4,5,10,20,50,100,200,300]

        total_genes = []
        scores = [] 

        reduced_num_genes = []
        reduced_scores = []

        for n in num_genes :
            genes = self.select_genes(n)

            pred = eval(genes)
            s = len(np.where( pred == pred_base )[0])/len(pred_base)

            total_genes.append( len(genes))
            scores.append(s)

            np.random.shuffle( genes )

            genes = genes[:int(0.9 * len(genes))]

            pred = eval(genes)
            s = len(np.where( pred == pred_base )[0])/len(pred_base)

            reduced_num_genes.append( len(genes) )
            reduced_scores.append( s )

        reduced_num_genes = np.array( reduced_num_genes )
        reduced_scores = np.array( reduced_scores )

        sinds = np.argsort( reduced_num_genes )

        reduced_num_genes = reduced_num_genes[sinds]
        reduced_scores = reduced_scores[sinds]

        pp.figure()
        pp.plot( total_genes, scores, '.-', label="No reduction" )
        pp.plot( reduced_num_genes, reduced_scores, '.-',label="Rand reduction")

        pp.xlabel('Num genes')
        pp.ylabel('Acc')

        pp.legend()
        pp.grid()

    def do_stuff4( self ):
        def eval( genes ):
            data = self.data['brain_data'].loc[:,genes]
            annots = self.data['brain_annots'].loc[:,'endotypes']

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
            
            return np.mean(diag)

        num_genes = [1,2,3,4,5,10,20,50,100,200,300]
        
        total_genes = []
        scores = [] 

        reduced_num_genes = []
        reduced_scores = []

        for n in num_genes :
            genes = self.select_genes(n)

            p = eval(genes)

            total_genes.append( len(genes))
            scores.append(p)

            np.random.shuffle( genes )

            genes = genes[:int(0.5 * len(genes))]

            p = eval(genes)

            reduced_num_genes.append( len(genes) )
            reduced_scores.append( p )

        
        pp.figure()
        pp.subplot()
        pp.plot( total_genes, scores, '.-', label="No reduction" )
        pp.subplot()
        pp.plot( reduced_num_genes, reduced_scores, '.-',label="Rand reduction")

        pp.xlabel('Num genes')
        pp.ylabel('Acc')

        pp.legend()
        pp.grid()

