import os
import glob
import pandas as pd
import numpy as np
from datetime import date


def get_unique_elements_for_each_key(input_dict):
    # Initialize a dictionary to store the unique elements for each key
    unique_dict = {}

    # Get all the values in the dictionary
    all_values = [value for values in input_dict.values() for value in values]

    # Iterate through the dictionary
    for key, values in input_dict.items():
        # Initialize a set to store the unique values for the current key
        unique_values = set(values)

        # Iterate through other values to remove non-unique items
        for other_key, other_values in input_dict.items():
            if key != other_key:
                unique_values -= set(other_values)

        # Convert the set to a list and add to the result dictionary
        unique_dict[key] = list(unique_values)

    return unique_dict


def main(save_folder, work_dir):
    path_degs = os.path.join(work_dir, '2024-01-17')
    # Get all .xlsx files in directory
    files = glob.glob(os.path.join(path_degs, '*.xlsx'))

    dict_degs = {}
    for file in files:
        names = file.split(os.sep)[-1].split('__')
        name_pos = names[0].split('DEGs_')[-1]
        name_neg = names[2]

        df_tmp = pd.read_excel(file)

        # Apply cut-off parameters padj < 0.05 and abs(log2FC) > 1
        mask = (np.abs(df_tmp['log2FoldChange']) > 1) & (df_tmp['padj'] < 0.05)
        degs_tmp = df_tmp.loc[mask, 'hgnc_symbol'].tolist()

        key = "{}__ vs __{}".format(name_pos, name_neg)
        dict_degs[key] = degs_tmp

    unique_dict_degs = get_unique_elements_for_each_key(input_dict=dict_degs)

    # create a dataframe for each endotype adding log2FoldChange and padj again
    # writer = pd.ExcelWriter(os.path.join(save_folder, 'Unique_DEGs_Dendrogram.xlsx'))
    for key in unique_dict_degs.keys():
        xlsx_filename = [s for s in files if key in s][0]
        df_init_tmp = pd.read_excel(xlsx_filename)

        df_new_tmp = pd.DataFrame(unique_dict_degs[key], columns=['hgnc_symbol'])
        df_result = pd.merge(df_new_tmp, df_init_tmp, on='hgnc_symbol', how='left')

        # save to .xlsx file
        df_result.to_excel(os.path.join(save_folder, 'Unique_DEGs_{}.xlsx'.format(key)))


if __name__ == '__main__':
    input_dir = os.path.join('/Volumes', 'CH__data', 'Projects', 'Eyerich_AG_projects',
                              'BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer', 'analysis', 'Molecular_subtypes',
                              'output', 'DGE_analysis', 'Endotypes_Dendrogram')
    output_dir = os.path.join(input_dir, 'unique_DEGs', str(date.today()))
    os.makedirs(output_dir, exist_ok=True)

    main(save_folder=output_dir, work_dir=input_dir)
