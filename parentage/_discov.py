#!/usr/bin/env python
# coding: utf-8
__author__ = "Igor Loschinin (igor.loschinin@gmail.com)"

import pandas as pd

"""
Search for paternity according to ICAR recommendations
https://www.icar.org/Documents/GenoEx/ICAR%20Guidelines%20for%20Parentage%20Verification%20and%20Parentage%20Discovery%20based%20on%20SNP.pdf
"""


class Discovery(object):
    """ Search for paternity according to ICAR recommendations

    :argument isag_markers: Fixed sample of markers to confirm paternity.
    """

    def __init__(
            self, isag_markers: pd.Series | list | set | None = None
    ) -> None:
        self.__isag_markers = isag_markers

        self.__num_conflicts = None  # Number of conflicts
        self.__perc_conflicts = None

    @property
    def status(self) -> None | str:
        """ The status of each parent discovered. """

        if self.__perc_conflicts is not None:
            if 0 <= self.__perc_conflicts < 1:
                return 'Discovered'
            elif 1 < self.__perc_conflicts < 3:
                return 'Doubtful'
            elif self.__perc_conflicts >= 3:
                return 'Excluded'
            else:
                return None

    @property
    def num_conflicts(self) -> None | int:
        return self.__num_conflicts

    @property
    def perc_conflicts(self) -> None | float:
        return self.__perc_conflicts

    def search_parent(
            self,
            data: pd.DataFrame,
            descendant: str,
            parents: str,
            snp_name_col: str
    ) -> None:
        """ Search for paternity.

        :param data: SNP data for descendant and parent.
        :param descendant: Columns name of the descendant in the data.
        :param parents: Columns name or list name of the parents in the data.
        :param snp_name_col: SNP columns name is data.
        """

        if self.__isag_markers is None:
            raise ValueError("Error. No array of snp names to verify")

        sample_by_markers = data.loc[
            data[snp_name_col].isin(self.__isag_markers),
            [snp_name_col, descendant, parents]
        ]

        # Filtering 5s from a descendent
        desc_marks = sample_by_markers.loc[
            sample_by_markers[descendant] != 5, [snp_name_col, descendant]
        ]

        # According to ICAR, the number of available markers must be
        # above 450
        if len(desc_marks) < 450:
            raise Exception("Calf call rate is low.")

        # Common after filtering markers of potential ancestors
        sample_parents = sample_by_markers.loc[
            sample_by_markers[snp_name_col].isin(desc_marks[snp_name_col]),
            parents
        ]

        # Number of available markers in potential ancestors
        prob_parents_same_n_markers = (sample_parents < 5).sum()

        # number of conflicts
        self.__num_conflicts = (
            abs(sample_parents.sub(desc_marks[descendant], axis=0)) == 2
        ).sum()

        # Percentage of conflicts
        self.__perc_conflicts = (
            (self.__num_conflicts / prob_parents_same_n_markers) * 100
        ).round(2)

    def __status_define(self) -> None:
        ...
