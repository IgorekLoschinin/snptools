#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd

"""
Search for paternity according to ICAR recommendations
https://www.icar.org/Documents/GenoEx/ICAR%20Guidelines%20for%20Parentage%20Verification%20and%20Parentage%20Discovery%20based%20on%20SNP.pdf
"""


class Discovery(object):
    """ Search for paternity according to ICAR recommendations """

    def __init__(
            self, isag_markers: pd.Series | list | set | None = None
    ) -> None:
        """
        :param isag_markers: fixed sample of markers to confirm paternity
        """
        self.__isag_markers = isag_markers

        self.__info = []
        self.__perc_conflicts = None  # Number of conflicts

    @property
    def status(self) -> None | str:
        if self.__perc_conflicts is not None:
            if 0 < self.__perc_conflicts < 1:
                return 'Discovered'
            elif 1 < self.__perc_conflicts < 3:
                return 'Doubtful'
            elif self.__perc_conflicts >= 3:
                return 'Excluded'
            else:
                return None

    @property
    def info(self) -> None | str:
        return

    @property
    def nun_conflicts(self) -> None | pd.DataFrame:
        return self.__perc_conflicts

    def search_pater(
            self,
            *,
            data: pd.DataFrame,
            descendant: str,
            parents: list[str],
            snp_name_col: str
    ) -> bool:
        """ Search for paternity

        :param data: - SNP data for descendant and parent
        :param descendant: - Columns name of the descendant in the data
        :param parents: - Columns name or list name of the parents in the data
        :param snp_name_col: - SNP columns name is data
        :return: -
        """

        if self.__isag_markers is None:
            raise Exception('Error. No array of snp names to verify')

        sample_by_markers = data.loc[
            data[snp_name_col].isin(self.__isag_markers),
            [snp_name_col, descendant, *parents]
        ]

        # Фильтруем 5ки у потомка
        desc_marks = sample_by_markers.loc[
            sample_by_markers[descendant] != 5, [snp_name_col, descendant]
        ]

        # Согласно ICAR кол-во доступных маркеров должно быть выше 450
        if len(desc_marks) < 450:
            self.__info.append('Calf call rate is low')
            return False

        # Общие после фильтрации маркеры потенциальных предков
        sample_parents = sample_by_markers.loc[
            sample_by_markers[snp_name_col].isin(desc_marks[snp_name_col]),
            parents
        ]

        # Количество доступных маркеров у потенциальных предков
        prob_parents_same_n_markers = (sample_parents < 5).sum()

        # Количество конфликтов
        n_conflicts = (
            abs(sample_parents.sub(desc_marks[descendant], axis=0)) == 2
        ).sum()

        # Процент конфликтов
        self.__perc_conflicts = (
            (n_conflicts / prob_parents_same_n_markers) * 100
        )

        return True


if __name__ == '__main__':
    df = pd.read_csv('parentage_test.csv', sep=" ")
    parentage_markers = pd.read_pickle("isag_disc.pl")

    obj = Discovery(
        isag_markers=parentage_markers.markers
    )
    obj.search_pater(
        data=df,
        descendant="BY000041988163",
        parents=["EE10512586"],
        snp_name_col="SNP_Name"
    )
