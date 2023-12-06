#!/usr/bin/env python
# coding: utf-8
__author__ = "Igor Loschinin (igor.loschinin@gmail.com)"

import numpy as np
import pandas as pd


"""
https://www.icar.org/Documents/GenoEx/ICAR%20Guidelines%20for%20Parentage%20Verification%20and%20Parentage%20Discovery%20based%20on%20SNP.pdf    
"""


class Verification(object):
    """
    Verification of paternity according to ICAR recommendations.

    :argument isag_marks: Fixed sample of markers to confirm paternity.
    """

    def __init__(
            self, isag_marks: pd.Series | list | set | None = None
    ) -> None:
        self.__isag_marks = isag_marks

        # The minimum number of SNP available in the profile
        # of each animal and potential parent must be scaled (i.e.: 95%
        # truncated down)
        self.__min_num_snp = 0.95
        self.__num_conflicts = None  # Number of conflicts

    @property
    def status(self) -> None | str:
        if self.__num_conflicts is not None:
            if self.__num_conflicts <= 2:
                return 'Accept'
            elif 3 <= self.__num_conflicts <= 5:
                return 'Doubtful'
            elif self.__num_conflicts > 5:
                return 'Excluded'
            else:
                return None

    @property
    def num_conflicts(self) -> None | int:
        return self.__num_conflicts

    def check_on(
            self,
            data: pd.DataFrame,
            descendant: str,
            parent: str,
            snp_name_col: str
    ) -> None:
        """ Verification of paternity according to ICAR recommendations.

        :param data: SNP data for descendant and parent.
        :param descendant: Columns name of the descendant in the data.
        :param parent: Columns name of the parent in the data.
        :param snp_name_col: SNP column name in data.
        """

        if self.__isag_marks is None:
            raise ValueError('Error. No array of snp names to verify')

        num_isag_mark = len(self.__isag_marks)
        min_num_comm_snp = int(num_isag_mark - (2 * (num_isag_mark * 0.05)))

        sample_mark = data.loc[
            data[snp_name_col].isin(self.__isag_marks), [descendant, parent]
        ]

        # The number of markers is not 5ok
        desc_n_markers = (sample_mark[descendant] < 5).sum()
        parent_n_markers = (sample_mark[parent] < 5).sum()

        # According to ICAR, the number of markers not 5ok should be more
        # than 95%
        if (desc_n_markers < num_isag_mark * self.__min_num_snp) and \
                (parent_n_markers < num_isag_mark * self.__min_num_snp):
            raise Exception('Calf and parent have low call rate')

        comm_snp_no_missing = sample_mark.replace(5, np.nan).dropna()
        num_comm_markers = len(comm_snp_no_missing)

        if num_comm_markers < min_num_comm_snp:
            raise Exception('Pair call rate is low')

        self.__num_conflicts = (abs(
            comm_snp_no_missing[descendant] - comm_snp_no_missing[parent]
        ) == 2).sum()
