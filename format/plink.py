#!/usr/bin/env python
# coding: utf-8
# https://zzz.bwh.harvard.edu/plink/data.shtml#ped

from finalreport import FinalReport
import pandas as pd


class Plink(object):

	def __init__(self) -> None:
		pass

	def handle(self) -> None:
		pass


"""
Создать map файл для Plink
https://www.cog-genomics.org/plink/1.9/formats#map
"""


def make_map(manifest):
    """
    Сделать map файл
    :param manifest: Dataframe манифест файла
    :return: Готовый Dataframe в формате map
    """
    manifest = manifest.sort_values(by='Name')
    # Переставить колонки и заменить названия половых и митохондриальной хромосомы
    permute_cols = manifest[['Chr', 'Name', 'MapInfo']].replace({'X': 30,
                                                                 'Y': 31,
                                                                 'MT': 33})
    # Вставить растояния в сантиморганидах
    permute_cols.insert(2, 'morgans', [0] * len(manifest))
    return permute_cols

#
# if __name__ == '__main__':
#     manifest_df = pd.read_csv('BovineSNP50_v3_A1_manifest.csv', usecols=['Chr', 'MapInfo', 'Name'])
#     make_map(manifest_df).to_csv('Bovine.map', sep='\t', header=False, index=False)



def make_ped(data, sex):
    replace_int = data.replace({0: 'A A',
                                1: 'A B',
                                2: 'B B',
                                5: '0 0'})
    replace_int.sort_index()
    ped_row = []
    for animal in replace_int.columns:
        markers = ' '.join(replace_int[animal].tolist())
        row = f'1 {animal} 0 0 {sex} 0 {markers}'
        ped_row.append(row)
    return '\n'.join(ped_row)


# if __name__ == '__main__':
#     df = pd.read_csv('../data_filter/filtered_data_Bovine50K.csv', index_col=0)
#     with open('Bovine.ped', 'w') as f:
#         f.write(make_ped(df, 1))
