from sklearn import preprocessing
import pandas as pd


def scale_data(data_train, sparse_data=True):
    """
    1. Only positive values in data
    2. For many features we only have a few data points => Sparse data

    Parameters
    ----------
    data_train
    sparse_data

    Returns
    -------

    """
    # TODO: Should we apply scaling on whole data set?
    #  -> No only on training data, use scaler to scale valid and test set
    if sparse_data:
        scaler = preprocessing.MaxAbsScaler().fit(data_train)
    else:
        scaler = preprocessing.StandardScaler().fit(data_train)

    data_scaled = scaler.transform(data_train)

    # convert to df
    df_data_scaled = pd.DataFrame(data=data_scaled, index=data_train.index, columns=data_train.columns)

    return df_data_scaled, scaler


def scale_gex(data_train, sparse_data=True):
    """
    1. Only positive values in data
    2. For many features we only have a few data points => Sparse data

    Parameters
    ----------
    data_train
    sparse_data

    Returns
    -------

    """
    # TODO: Should we apply scaling on whole data set?
    #  -> No only on training data, use scaler to scale valid and test set
    if sparse_data:
        scaler = preprocessing.MaxAbsScaler().fit(data_train)
    else:
        scaler = preprocessing.StandardScaler().fit(data_train)

    data_scaled = scaler.transform(data_train)

    # convert to df
    df_data_scaled = pd.DataFrame(data=data_scaled, index=data_train.index, columns=data_train.columns)

    return df_data_scaled, scaler

