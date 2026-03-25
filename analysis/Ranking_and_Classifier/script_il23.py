import numpy as np 
import pandas as pd 
import scanpy
from easydict import EasyDict as edict
from tqdm import tqdm

from config import Config
import umap
from dataframe_tools import partition

from modeling import utils
from modeling import FeatureModels
from sklearn.linear_model import LinearRegression

from matplotlib import pyplot as pp 
import matplotlib.colors as mcolors
pp.ion()

endotypes = ['E1','E2','E3','E4','E5','E6','E7','E8','E9','E10','E11','E12','E13']

def eval_function( x,y ):
    return np.mean( (x-y) ** 2, axis=1 )

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

class CSVData :
    def __init__( self, data_name ):
        self.expression_data = pd.read_csv(f'./{data_name}/expression_data.csv', index_col=0)
        self.meta_data = pd.read_csv(f'./{data_name}/metadata.csv', index_col=0)
        self.meta = self.meta_data.loc[ :, ['title', 'subject id','pasi75resp'] ]

        self.expression_data = self.expression_data
        self.ids = list(self.expression_data.columns)

    def get_data( self, genes, gene_map, annot_count ):
        index = self.expression_data.columns

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

class SampleClassifier :
    def __init__( self ):
        self.models = {}

    def create_binary_labels( self, labels, target ):

        ll = np.zeros_like( labels )

        pos_inds = np.where( labels == target )[0]
        neg_inds = np.where( labels != target )[0] 

        ll[pos_inds] = 1 
        ll[neg_inds] = -1

        return ll

    def train( self, data, labels ):
        labels = labels.astype( int )

        num_labels = np.max( labels )+1 

        self.models = {}

        #keep_labels = [ 7, 10, 11, 13 ]

        for i in tqdm( range( num_labels ) ):

            #if not i in keep_labels :
            #    continue

            b_labels = self.create_binary_labels( labels, i ) 

            m = LinearRegression()
            m.fit( data, b_labels )

            self.models[i] = m

    def eval( self, s ):
        
        indices = []
        scores = []
        for idx , m in self.models.items() :

            scores.append(m.predict(s)[0])
            indices.append(idx)

        scores = np.array( scores )

        midx = np.argmax(scores)

        return indices[midx]

class Module :
    def __init__( self ):
        args = edict()
        args.rank_model = 'LinearRegression'
        args.feat_model = None 

        self.cfg = Config(args)

        self.brain_data, self.brain_annots = load_brain( self.cfg )

        self.data_il23 = CSVData('data-il-23')
        self.data_tnf = CSVData('data-tnf')

        self.parts = partition( self.brain_annots, "endotypes", endotypes, 
                                self.cfg.parts_ratio, self.cfg.parts_num,
                                self.cfg.sampling_ratio, self.cfg.sampling_num,
                                numpy_seed=self.cfg.seed )

    def create_gene_scores( self ):

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

    def create_gene_map( self ):
        all_genes = list(self.brain_data.columns)
        self.gene_map = {}
        for gene in tqdm(all_genes) :
            self.gene_map[gene] = self.ext_genes_and_ids( self.data_il23, gene )

     
    def ext_genes_and_ids( self, data, key ):
        genes = list(data.expression_data.index)

        selected_genes = []
        for gene in genes :

            if '...' in gene :
                continue

            if key in gene :
                selected_genes.append(gene)

        return selected_genes
        #self.selected_genes = selected_genes

    def get_ids( self, respond ):
        ids_nl_week0 = []
        ids_ls_week0 = []
        ids_ls_week12 = []

        ids = np.array( list(self.meta_data.index) )
        titles = self.meta_data.title

        self.meta = self.meta_data.loc[ :, ['title', 'subject id','pasi75resp'] ]

        for id in ids :
            t = self.meta.loc[id,'title']
            r = self.meta.loc[id,'pasi75resp']

            if 'NL Week 0' in t and r == respond :
                ids_nl_week0.append( id )
            if 'LS Week 0' in t and r == respond :
                ids_ls_week0.append( id )
            if 'LS Week 12' in t and r == respond :
                ids_ls_week12.append( id )

        ids_nl_week0 = np.array( ids_nl_week0 )
        ids_ls_week0 = np.array( ids_ls_week0 )
        ids_ls_week12 = np.array( ids_ls_week12 ) 

        return ids_nl_week0, ids_ls_week0, ids_ls_week12

    def get_genes( self, num, valid_endotypes=None ):

        genes = []

        all_genes = list(self.brain_data.columns)
        for idx, endotype in enumerate(endotypes):

            if valid_endotypes is not None :
                if endotype not in valid_endotypes :
                    continue

            scores = self.gene_scores[:,idx].ravel()
            sinds = np.argsort( scores )

            e_gene_idx = 0
            e_genes = []

            for s_idx in sinds :
                g = all_genes[s_idx]
                
                if len(self.gene_map[g]) > 0 :
                    e_genes.append(g)

                if len(e_genes) == num :
                    break

            genes += e_genes 

        return genes            

    def normalize_data( self, data ):

        data = data.to_numpy().T 

        for idx, r in enumerate(data):
            m = np.mean(r)
            s = np.std(r)

            data[idx,:] = (r - m)/s


        return data.T

    def normalize_data_pd( self, data ):
        index = data.index 
        columns = data.columns

        data = data.to_numpy().T 

        for idx, r in enumerate(data):
            m = np.mean(r)
            s = np.std(r)
        
            data[idx,:] = (r - m)/s

        return pd.DataFrame( data.T, index=index, columns=columns )

    def normalize_data_2( self, data, annot, vals ):

        inds = []

        for v in vals :
            ii = np.where( annot == v )[0]
            inds = np.array( ii )

        inds = np.concatenate(inds)

        sub_data = data[inds,:]

        sub_data = self.normalize_data( sub_data )
        data = data.to_numpy().T 

        for idx, r in enumerate(data):
            m = np.mean(r)
            s = np.std(r)

            data[idx,:] = (r - m)/s

        return data.T

    def normalize_data3( self, data ):
        data = data.to_numpy()

        for idx, r in enumerate(data):
            data[idx] = r / np.linalg.norm(r)

        return data 

    def do_stuff( self ):
        # E11 E13 E8

        legend_arr = []

        for e in endotypes :
            legend_arr.append( e )

        legend_arr.append('NL Week0 - Resp : Yes (IL23)')
        legend_arr.append('LS Week0 - Resp : Yes (IL23)')
        legend_arr.append('LS Week12 - Resp : Yes (IL23)')
        legend_arr.append('NL Week0 - Resp : No (IL23)')
        legend_arr.append('LS Week0 - Resp : No (IL23)')
        legend_arr.append('LS Week12 - Resp : No (IL23)')

        legend_arr.append('NL Week0 - Resp : Yes (TNF)')
        legend_arr.append('LS Week0 - Resp : Yes (TNF)')
        legend_arr.append('LS Week12 - Resp : Yes (TNF)')
        legend_arr.append('NL Week0 - Resp : No (TNF)')
        legend_arr.append('LS Week0 - Resp : No (TNF)')
        legend_arr.append('LS Week12 - Resp : No (TNF)')

        markers = ['.', '*', '+', '1', '2' ]
        
        genes = self.get_genes( 10 )
        brain_data = self.brain_data.loc[:,genes]

        annot_count = 0
        ba = self.brain_annots.loc[:,'endotypes'].to_numpy()

        labels = np.zeros( len(ba) )

        for idx, e in enumerate(endotypes) :
            inds = np.where( ba == e )[0]
            labels[ inds ] = idx 
            annot_count += 1

        brain_annots = labels

        il23_data, il23_annots = self.data_il23.get_data( genes, self.gene_map, annot_count )

        tnf_data, tnf_annots = self.data_tnf.get_data( genes, self.gene_map, annot_count+6 )

        other_annots = tnf_annots

        brain_data = self.normalize_data( brain_data )
        il23_data = self.normalize_data( il23_data )
        tnf_data = self.normalize_data( tnf_data )

        #other_data = self.normalize_data_pd( pd.concat([tnf_data,il23_data]) ) 
        #il23_data = other_data.loc[ self.data_il23.ids ].to_numpy()
        #tnf_data = other_data.loc[ self.data_tnf.ids ].to_numpy()

        data = np.concatenate( [brain_data, il23_data, tnf_data ], axis=0 )
        annots = np.concatenate( [brain_annots, il23_annots, tnf_annots ] ).astype( int )


        colormap = pp.cm.autumn  # You can choose from various colormaps
        colors = [colormap(i) for i in np.linspace(0, 1,len(legend_arr)+2 )]
        

        sub_data_names = ['E8','E11','E12','E13','LS Week0 - Resp : Yes', 'LS Week0 - Resp : No' ]

        reducer = umap.UMAP()
        embedding = reducer.fit_transform( data )

        pp.figure()

        c = 0 

        ignore_list =  [ 'NL Week0 - Resp : Yes (IL23)', 'LS Week12 - Resp : Yes (IL23)',
                         'NL Week0 - Resp : No (IL23)', 'LS Week12 - Resp : No (IL23)',
                         'NL Week0 - Resp : Yes (TNF)', 'LS Week12 - Resp : Yes (TNF)',
                         'NL Week0 - Resp : No (TNF)', 'LS Week12 - Resp : No (TNF)' ]

        target_endotypes = ['E8', 'E11', 'E12','E13', 'LS Week0 - Resp : Yes (IL23)', 'LS Week0 - Resp : No (IL23)', 'LS Week0 - Resp : Yes (TNF)', 'LS Week0 - Resp : No (TNF)']

        for idx in range( np.max(annots)+1 ):
            inds = np.where( annots == idx )[0]
                
            if len(inds) == 0 :
                continue

            d = embedding[inds,:]
            a = annots[inds]

            annot_idx = a[0]
            if legend_arr[annot_idx] in ignore_list :
                continue
  
            
            m = markers[ idx % len(markers)]
            
            if legend_arr[annot_idx] in target_endotypes : 
                pp.plot( d[:,0], d[:,1], m, label=legend_arr[annot_idx] )#, color=colors[c] )
            else :
                pp.plot( d[:,0], d[:,1], m, label=legend_arr[annot_idx], color='black' )


            c += 1

        pp.legend( fontsize=4 )
        pp.grid()

        pp.savefig('plot_il23_tnf.pdf')

    def do_stuff_classifier( self, retrain_classifier=False ):
        legend_arr = []

        for e in endotypes :
            legend_arr.append( e )

        legend_arr.append('NL Week0 - Resp : Yes (IL23)')
        legend_arr.append('LS Week0 - Resp : Yes (IL23)')
        legend_arr.append('LS Week12 - Resp : Yes (IL23)')
        legend_arr.append('NL Week0 - Resp : No (IL23)')
        legend_arr.append('LS Week0 - Resp : No (IL23)')
        legend_arr.append('LS Week12 - Resp : No (IL23)')

        legend_arr.append('NL Week0 - Resp : Yes (TNF)')
        legend_arr.append('LS Week0 - Resp : Yes (TNF)')
        legend_arr.append('LS Week12 - Resp : Yes (TNF)')
        legend_arr.append('NL Week0 - Resp : No (TNF)')
        legend_arr.append('LS Week0 - Resp : No (TNF)')
        legend_arr.append('LS Week12 - Resp : No (TNF)')

        markers = ['.', '*', '+', '1', '2' ]
        
        genes = self.get_genes( 5 )
        brain_data = self.brain_data.loc[:,genes]

        annot_count = 0
        ba = self.brain_annots.loc[:,'endotypes'].to_numpy()

        labels = np.zeros( len(ba) )

        for idx, e in enumerate(endotypes) :
            inds = np.where( ba == e )[0]
            labels[ inds ] = idx 
            annot_count += 1

        brain_annots = labels

        il23_data, il23_annots = self.data_il23.get_data( genes, self.gene_map, annot_count )

        tnf_data, tnf_annots = self.data_tnf.get_data( genes, self.gene_map, annot_count+6 )

        other_annots = tnf_annots

        brain_data = self.normalize_data( brain_data )
        il23_data = self.normalize_data( il23_data )
        tnf_data = self.normalize_data( tnf_data )

        print('Training model')

        ignore_list =  [ 'NL Week0 - Resp : Yes (IL23)', 'LS Week12 - Resp : Yes (IL23)',
                         'NL Week0 - Resp : No (IL23)', 'LS Week12 - Resp : No (IL23)',
                         'NL Week0 - Resp : Yes (TNF)', 'LS Week12 - Resp : Yes (TNF)',
                         'NL Week0 - Resp : No (TNF)', 'LS Week12 - Resp : No (TNF)' ]

        print( brain_data.shape )

        if retrain_classifier :
            self.classifier = SampleClassifier()
            self.classifier.train( brain_data, brain_annots ) 

        keys = [ 14, 17, 20, 23 ]

        res = np.zeros([len(keys), 13])

        for s,l in zip( il23_data, il23_annots ):
            l = int(l)

            if legend_arr[l] in ignore_list :
                continue 

            c = self.classifier.eval( s.reshape([1,-1]) ) 

            res[keys.index(l), c] += 1

        for s,l in zip( tnf_data, tnf_annots ):
            l = int(l)

            if legend_arr[l] in ignore_list :
                continue 

            c = self.classifier.eval( s.reshape([1,-1]) ) 

            res[keys.index(l), c] += 1


        res_dp = pd.DataFrame( res, 
                               index=[ legend_arr[idx] for idx in keys ],
                               columns=endotypes
                              )

        print( res_dp )
        print( res_dp.sum( axis=1 ))
        res_dp.to_csv('il23_tnf.csv')

        


