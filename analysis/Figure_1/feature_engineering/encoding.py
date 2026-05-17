
from sklearn.preprocessing import OrdinalEncoder, OneHotEncoder, LabelEncoder

import pandas as pd
import numpy as np


def encode_label(label):
    """Works with numeric data and strings
    Note: LabelEncode does not encode the order, this needs to be done manually

    Parameters
    ----------
    label : numpy.array, pandas.Series

    Returns
    -------

    """
    label_encoder = LabelEncoder()
    label_encoder.fit(label)
    return label_encoder


def encode_ordinal_categories(data):
    """Before imputing categorical variables using fancyimpute you have to encode the strings to numerical values.
        -> then call function `impute_numerical(func_data_numeric)`
    OrdinalEncoder converts also string categories to numerical categories

    Note: OrdinalEncoder does not encode the order, this needs to be done manually

    see: https://journalofbigdata.springeropen.com/articles/10.1186/s40537-020-00305-w and Python Scikit-learn library

    Parameters
    ----------
    data : pandas.Dataframe

    Returns
    -------

    """
    # Factorize -> does not contain the order, e.g. 3 could be encoded as 0
    # bulk_data.obs[ordinal_cols]
    # bulk_data.obs[ordinal_cols[0]].mask(bulk_data.obs[ordinal_cols[0]] == ' ').factorize()[0]

    # Encode Ordinal features with OrdinalEncoder
    # convert dataframe to type object
    data = data.astype('object')

    # Build encoder
    encoder = OrdinalEncoder()

    for ordinal_features in data.columns:
        # data[ordinal_features] = data[ordinal_features].cat.reorder_categories(data[ordinal_features].cat.categories.sort_values())
        nonans = data[ordinal_features].dropna().values
        impute_reshape = np.asarray(nonans).reshape(-1, 1)
        data.loc[data[ordinal_features].notna(), [ordinal_features]] = encoder.fit_transform(impute_reshape)

        # fill up NaN with np.nan
        data.loc[:, ordinal_features] = data.loc[:, ordinal_features].fillna(np.nan)
        # convert to type category
        data.loc[:, ordinal_features] = data.loc[:, ordinal_features].astype('category')

    return data


def encode_nominal_feature(data, dummy_variable_encoding=False):
    """Apply One-Hot encoding on nominal categories
    https://machinelearningmastery.com/one-hot-encoding-for-categorical-data/
    https://www.ritchieng.com/machine-learning-evaluate-linear-regression-model/

    Parameters
    ----------
    data : pandas.Dataframe
    dummy_variable_encoding : bool
        If True, data will by encoded using pandas.get_dummies or OneHotEncoder(drop='first, sparse=False) which drops
        the first column of each feature
        If False: data will be encoded with OneHotEncoder(sparse=False)

    Returns
    -------

    """
    if dummy_variable_encoding:
        # define one hot encoding
        encoder = OneHotEncoder(drop='first', sparse=False)
        # or pd.get_dummies(data, columns=col).head()
    else:
        # define one hot encoding
        encoder = OneHotEncoder(sparse=False)

    # transform data
    onehot_data = encoder.fit_transform(data)
    encoded_feature_names = encoder.get_feature_names_out(data.columns)

    # convert back to dataframe
    df_onehot_data = pd.DataFrame(onehot_data, columns=encoded_feature_names)

    return df_onehot_data


def imputation_replace_zero_nan(df_data, original_nominal_feature_names):
    """
    Step 4:
    For each row and each original categorical feature, when the 1 is in the "missing" column, replace the 0's with
    np.nan; then delete the missing indicator column.

    Parameters
    ----------
    df_data : pandas.Datafrane
    original_nominal_feature_names : list of strings

    Returns
    -------

    """
    for original_feature_name in original_nominal_feature_names:
        # get all columns with original_feature_name
        encoded_cols = [for_col for for_col in df_data.columns if original_feature_name in for_col]
        # get column of nan category
        matching_nan = [s for s in encoded_cols if "nan" in s]
        # Only do that step if you have missing values in feature
        if len(matching_nan) > 0:
            # get other columns which don't encode nan category
            encoded_cols_notnan = encoded_cols.copy()
            encoded_cols_notnan.remove(matching_nan[0])
            # find index where 1 is in nan category
            ind_nancolumn = np.where(df_data[matching_nan] == 1)[0]
            # replace 0 with NaN in other categories
            df_data.loc[ind_nancolumn, encoded_cols_notnan] = np.nan
            # delete the nan category column
            df_data = df_data.drop(matching_nan, axis=1)

    df_data[df_data.columns] = df_data[df_data.columns].astype('category')

    return df_data
