from sklearn.impute import KNNImputer

import pandas as pd


def impute_knn(func_data_numeric, ordinal_categories: list, k: int = 1):
    """ Impute missing values using KNN
    https://krrai77.medium.com/using-fancyimpute-in-python-eadcffece782

    Notes:
    - Can become computationally difficult with more predictor variables and more instances (doesn’t scale well)
    - CAN work for categorical variables, but will require transformation to dummy variable(s),
     in the case of nominal categorical data OR numeric conversion in the case of ordinal data.

    Parameters
    ----------
    func_data_numeric : pandas.Dataframe
    ordinal_categories:
    k : int

    Returns
    -------

    """
    # Instantiate KNN imputer
    knn_imputer = KNNImputer(n_neighbors=k, weights='uniform')
    # imputing the missing value with knn imputer
    array_imputed = knn_imputer.fit_transform(func_data_numeric)
    # convert to dataframe:
    df_func_data_fillna = pd.DataFrame(array_imputed, index=func_data_numeric.index, columns=func_data_numeric.columns)

    # convert objects back to categories
    for cat in ordinal_categories:
        df_func_data_fillna = df_func_data_fillna.astype({cat: 'category'})

    return df_func_data_fillna
