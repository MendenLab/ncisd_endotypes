import numpy as np 
import pandas as pd

def create_binary_labels( y, classes ):

    y = y.to_numpy()

    labels = []
    for c in classes :
        l = np.ones_like( y )
        neg = np.where( y != c )[0]
        l[neg] = -1
        labels.append( l.reshape([1,-1]))

    labels = np.concatenate( labels, axis=0 )
    return labels

def create_binary_labels_geo( y, pos_idx ):
    y = y.to_numpy()

    labels = []

    l = np.ones_like( y )
    neg = np.where( y != pos_idx )[0]
    l[neg] = -1

    l = l.reshape([1,-1])
    return l

def create_binary_labels_new( y, positives, negatives ):
    index = y.index
    y = y.to_numpy()
    labels = np.zeros_like( y )

    for p in positives :
        inds = np.where( y == p )[0]
        labels[inds] = 1

    for n in negatives :
        inds = np.where( y == n )[0]
        labels[inds] = -1

    return pd.DataFrame( labels, index=index, columns=['labels'] ) 

def create_binary_labels_np( y, positives, negatives ):
    index = y.index
    y = y.to_numpy()
    labels = np.zeros_like( y )

    for p in positives :
        inds = np.where( y == p )[0]
        labels[inds] = 1

    for n in negatives :
        inds = np.where( y == n )[0]
        labels[inds] = -1

    return labels
    #return pd.DataFrame( labels, index=index, columns=['labels'] ) 




