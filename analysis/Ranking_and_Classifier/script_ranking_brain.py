import numpy as np 
import pandas as pd
from easydict import EasyDict as edict 
from tqdm import tqdm

from sklearn.linear_model import LinearRegression
from sklearn.svm import SVC,SVR

from matplotlib.colors import LinearSegmentedColormap

from config import Config
from tools import load_brain
from dataframe_tools import partition
from modeling import FeatureModels

import matplotlib.pyplot as pp
pp.ion()
import seaborn as sn

features = {}

features['endotypes'] = ['E1','E2','E3','E4','E5','E6','E7','E8','E9','E10','E11','E12','E13']
features['pattern'] = ['1','2a', '3']

def eval_function( x,y ):
    return np.mean( (x-y) ** 2, axis=1 )

class Module :
    def __init__( self, norm ):
        args = edict()
        args.rank_model = 'LinearRegression'
        args.feat_model = None 

        self.norm = norm

        self.cfg = Config(args)
        self.brain_data, self.brain_annots = load_brain(self.cfg, norm=self.norm)

    def create_parts( self, label ):
        self.label = label 

        self.parts = partition( self.brain_annots, label, features[label], 
                                self.cfg.parts_ratio, self.cfg.parts_num,
                                self.cfg.sampling_ratio, self.cfg.sampling_num,
                                numpy_seed=self.cfg.seed )

    def calc_gene_scores( self ):
        data = self.brain_data
        annots = self.brain_annots[self.label]
        classifier = 'LinearRegression'
        ff = features[self.label]

        labels = []

        for f in ff :
            aa = (annots == f).astype(float).to_numpy()
            aa = aa * 2 - 1
            aa = aa.reshape([-1,1])
            labels.append( aa )

        labels = np.concatenate( labels, axis=1 )
        labels = pd.DataFrame( labels, index=annots.index, columns=ff )

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

        self.gene_scores = gene_scores

    def category_labels( self, labels ):
        out = []

        for l in labels :
            out.append( features[self.label].index(l) )

        return np.array( out )

    def binary_labels( self, labels, endotype ):
        
        ll = np.zeros( len(labels) )

        pos_inds = np.where( labels == endotype )[0]
        neg_inds = np.where( labels != endotype )[0]

        ll[pos_inds] = 1.0 
        ll[neg_inds] = -1.0

        return ll

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


    def select_genes( self, num, save=False ):

        ff = features[self.label]

        genes = list(self.brain_data.columns)
        gene_set = set()

        genes_df = []
        for f, scores in zip( ff, self.gene_scores.T ) :
            sinds = np.argsort( scores )
    
            f_genes = []
            for i in sinds[:num] :
                gene_set.add( genes[i] )
                f_genes.append( genes[i] )

            genes_df.append( f_genes )

        if save :
            genes_df = np.array(genes_df).T
            genes_df = pd.DataFrame( genes_df, columns=ff )
            genes_df.to_excel(f'genes_{num}.xlsx')

        selected_genes = np.array(list(gene_set))

        return selected_genes

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
            
            m = SVC( kernel='linear', C=1.5, class_weight='balanced' )
            m.fit( X_tr, l_tr )
            pred = m.predict(X_te)

            #np.random.shuffle(pred)

        if use_linear : 
            scores = [] 
    
            for idx, ff in enumerate( features[self.label] ):
                l_tr = self.binary_labels( y_tr, ff )

                #m = LinearRegression()
                m = SVR()
                m.fit( X_tr, l_tr )

                #r = l_tr.copy()
                #np.random.shuffle(r)
        
                s = m.predict( X_te )
                scores.append( r )

            scores = np.array( scores ).T
            pred = np.argmax( scores, axis=1 )

        ff = features[self.label]

        cmat = np.zeros(( len(ff), len(ff)) )

        np.random.shuffle(pred)

        for l,p in zip( y_te, pred ):
            l = ff.index(l)

            cmat[l,p] += 1

        return cmat 

    def normalize_cmat( self, cmat ):
        for idx, r in enumerate(cmat):
            r = r / r.sum()
            cmat[idx] = r
        return cmat

    def eval( self, num_genes, show_figure=False, show_genes=False ):
        sn.set(font_scale=1.0)

        genes = self.select_genes(num_genes)

        if show_genes :
            for gene in genes :
                print( gene )

        data = self.brain_data.loc[:,genes]
        annots = self.brain_annots.loc[:,self.label]

        cmats = []

        for part in tqdm(self.parts) :
            cmat = self.do_stuff_part( data, annots, part )
            cmats.append(cmat)

        cmats = np.array( cmats )
        cmats = np.mean( cmats, axis=0 )

        cmats = self.normalize_cmat( cmats )

        ff = features[self.label]

        diag = np.diag(cmats)
        df_conf_mat = pd.DataFrame( cmats, index=ff, columns=ff )

        if show_figure :
            #colors = [(0, 0.3, 0), (0, 0, 0), (0.5, 0, 0)]  # Green → Black → Red
            #colors = [(0, 0.7, 0), (0, 0, 0), (1.0, 0, 0)]  # Green → Black → Red

            #cmap_name = ''
            #cm = LinearSegmentedColormap.from_list(cmap_name, colors, N=256)

            pp.figure(figsize = (10,10))
            sn.heatmap(df_conf_mat, annot=True,cmap="coolwarm")

            pp.title(f'Mean acc : {np.round(np.mean(diag),3)} $\\pm$ {np.round(np.std(diag),3)} - Num genes : {len(genes)} - Norm : {self.norm}')
            pp.savefig(f'conf_mat_{self.label}_{num_genes}_rand.pdf')

        return len(genes), np.mean(diag), np.std(diag)

    def eval_multiple( self ):

        num_genes = []
        acc = []

        for i in [ 1,5,10,20,30,40,50,100,150,200,300,400 ] :
            n,a = self.eval(i)

            num_genes.append(n)
            acc.append(a)

        pp.figure()

        pp.plot( num_genes, acc, '.-' )
        pp.xlabel('Number of genes')
        pp.ylabel('Classifiation Accuracy')
        pp.grid()
        pp.title(f'Classifiction Acc - Norm : {self.norm}')

        pp.plot([num_genes[7]], [acc[7]],'r*', label='100 genes per endotype')

        pp.legend()
        pp.savefig(f'acc_{self.label}.pdf')

    def cluster( self, num_genes ):
        genes = self.select_genes(num_genes)
        data = self.brain_data.loc[:,genes].T
        annots = self.brain_annots.loc[:,self.label]

        ff = features[self.label]

        pooled_data = []

        for f in ff :
            a = annots[ annots == f ]
            index = a.index 

            df = data[index]

            df = df.mean(axis=1).to_numpy()
            pooled_data.append(df)

        pooled_data = np.array( pooled_data ).T

        df_pooled = pd.DataFrame( pooled_data, index=genes, columns=ff)

        sn.set(font_scale=1.0)

        #colors = [(0, 0.7, 0), (0, 0, 0), (1.0, 0, 0)]  # Green → Black → Red
        #cmap_name = 'green_black_red'
        #cm = LinearSegmentedColormap.from_list(cmap_name, colors, N=256)

        print( df_pooled.shape )


        cluster = sn.clustermap(df_pooled, standard_scale=0, cmap="coolwarm", 
                      method="complete", metric="euclidean",
                      xticklabels=1, yticklabels=1,
                      figsize=(13 * 0.3,65 * 0.3), col_cluster=False)

        cluster.cax.set_ylabel("Normalized Expression over endotypes", fontsize=5, labelpad=30, rotation=270)

        xtick_colors = {"Sample1": "red", "Sample2": "blue", "Sample3": "green", "Sample4": "purple"}

        colormap = pp.cm.hsv  # You can choose from various colormaps
        colors = [colormap(i) for i in np.linspace(0, 1,len(ff) )]
        np.random.shuffle(colors)

        print( colors )


        #for label in cluster.ax_heatmap.get_xticklabels():
        #    t = label.get_text()
        #    f = ff.index(annots.loc[t])
        #    label.set_color( colors[f] ) 

        pp.savefig('cluster.pdf')


        

if __name__=="__main__" :

    m = Module()
    m.create_parts('pattern')
    m.calc_gene_scores()
    
    m.eval(1)
    m.eval(5)
    m.eval(10)
