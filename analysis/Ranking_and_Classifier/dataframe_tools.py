import pandas as pd 
import numpy as np

from numpy_tools import NumpyTemporarySeed 

def partition_df( data, column, groups, ratio, num ):

    parts = []

    for part_idx in range( num ): 

        train = []
        test = []

        for g in groups : 
            indices = data[ data.loc[:,column] == g ].index
            indices = indices.to_numpy()

            np.random.shuffle( indices )

            nindices = len(indices)

            train_group = indices[:int(nindices*ratio)]
            test_group = indices[int(nindices*ratio):]
                
            train.append( train_group )
            test.append( test_group )

        train = np.concatenate( train )
        test = np.concatenate( test )

        parts.append( {'train':train,'test':test} )

    return parts


def partition( data, column, groups, 
               parts_ratio, parts_num,
               sampling_ratio, sampling_num,
               numpy_seed=1234 ):

    with NumpyTemporarySeed( numpy_seed ):
        parts = partition_df( data, column, groups, parts_ratio, parts_num )

        for part in parts :
            data_tr = data.loc[ part['train'] ]
            sampling = partition_df( data_tr, column, groups, sampling_ratio, sampling_num )

            sampling_arr = []
            for ss in sampling :
                sampling_arr.append( ss['train'] )

            part['sampling'] = sampling_arr

    return parts
