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
    path_pas = os.path.join(work_dir, '2023-09-13')
    # Get all .xlsx files in directory
    files = glob.glob(os.path.join(path_pas, 'GSEA*.xlsx'))

    dict_pas = {}
    for file in files:
        names = file.split(os.sep)[-1].split('__')
        name_pos = names[1]
        name_neg = names[3]

        df_tmp = pd.read_excel(file)

        # Apply cut-off parameters padj < 0.05 and abs(log2FC) > 1
        mask = df_tmp['p.adjust'] < 0.05
        pas_tmp = df_tmp.loc[mask, 'Description'].tolist()

        key = "{}__ vs __{}".format(name_pos, name_neg)
        dict_pas[key] = pas_tmp

    unique_dict_pas = get_unique_elements_for_each_key(input_dict=dict_pas)

    # create a dataframe for each endotype adding log2FoldChange and padj again
    # writer = pd.ExcelWriter(os.path.join(save_folder, 'Unique_DEGs_Dendrogram.xlsx'))
    for key in unique_dict_pas.keys():
        xlsx_filename = [s for s in files if key in s][0]
        df_init_tmp = pd.read_excel(xlsx_filename)

        df_new_tmp = pd.DataFrame(unique_dict_pas[key], columns=['Description'])
        df_result = pd.merge(df_new_tmp, df_init_tmp, on='Description', how='left')

        # save to .xlsx file
        df_result.to_excel(os.path.join(save_folder, 'Unique_GSEA_{}.xlsx'.format(key)))


def main_ora(save_folder, work_dir):
    path_pas = os.path.join(work_dir, '2023-09-13')
    # Get all .xlsx files in directory
    files = glob.glob(os.path.join(path_pas, 'ORA*.xlsx'))

    dict_pas = {}
    for file in files:
        names = file.split(os.sep)[-1].split('__')
        name_pos = names[1]
        name_comp_1 = names[2]
        name_comp_2 = names[4]

        df_tmp = pd.read_excel(file)

        # Apply cut-off parameters padj < 0.05 and abs(log2FC) > 1
        mask = df_tmp['p.adjust'] < 0.05
        pas_tmp = df_tmp.loc[mask, 'Description'].tolist()

        key = "{}__{}__ vs __{}".format(name_pos, name_comp_1, name_comp_2)
        dict_pas[key] = pas_tmp

    unique_dict_pas = get_unique_elements_for_each_key(input_dict=dict_pas)

    # create a dataframe for each endotype adding log2FoldChange and padj again
    # writer = pd.ExcelWriter(os.path.join(save_folder, 'Unique_DEGs_Dendrogram.xlsx'))
    for key in unique_dict_pas.keys():
        xlsx_filename = [s for s in files if key in s][0]
        df_init_tmp = pd.read_excel(xlsx_filename)

        df_new_tmp = pd.DataFrame(unique_dict_pas[key], columns=['Description'])
        df_result = pd.merge(df_new_tmp, df_init_tmp, on='Description', how='left')

        # save to .xlsx file
        df_result.to_excel(os.path.join(save_folder, 'Unique_ORA_{}.xlsx'.format(key)))


if __name__ == '__main__':
    input_dir = os.path.join('/Volumes', 'CH__data', 'Projects', 'Eyerich_AG_projects',
                              'BRAIN__Natalie_Garzorz-Stark_Peter_Seiringer', 'analysis', 'Molecular_subtypes',
                              'output', 'Pathways', 'Endotypes_Dendrogram')
    output_dir = os.path.join(input_dir, 'unique_Pas', str(date.today()))
    os.makedirs(output_dir, exist_ok=True)

    main(save_folder=output_dir, work_dir=input_dir)
    main_ora(save_folder=output_dir, work_dir=input_dir)
