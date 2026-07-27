#!/usr/bin/env python
# coding: utf-8
__author__ = "Igor Loschinin (igor.loschinin@gmail.com)"

import numpy as np
import pandas as pd


"""
https://www.icar.org/Documents/GenoEx/ICAR%20Guidelines%20for%20Parentage%20Verification%20and%20Parentage%20Discovery%20based%20on%20SNP.pdf
https://www.icar.org/wp-content/uploads/documents/Guidelines-for-Parentage-Verification-and-Parentage-Discovery-Based-on-SNP-Genotypes-Feb.-2025.pdf    
"""


class Verification(object):
    """
    Verification of paternity according to ICAR recommendations.

    :argument isag_marks: Fixed sample of markers to confirm paternity.
    """

    def __init__(
            self,
            isag_marks: pd.Series | list | set | None = None
    ) -> None:
        self._isag_marks = isag_marks

        # ----- Attributes for Step 1 (Checking one parent) ----
        # The minimum number of SNP available in the profile
        # of each animal and potential parent must be scaled (i.e.: 95%
        # truncated down)
        self._min_num_snp = 0.95
        self._num_conflicts_sing = None  # Number of conflicts
        self._status_sing = None

        # ----- Attributes for Step 2 (Checking Trio/Mating) ----
        self._status_mat = None
        self._num_conflicts_mat = None
        self._num_common_mat = None

    @property
    def status_sing(self) -> None | str:
        return self._status_sing

    @property
    def num_conflicts_sing(self) -> None | int:
        return self._num_conflicts_sing

    @property
    def status_mat(self) -> None | str:
        return self._status_mat

    @property
    def num_conflicts_mat(self) -> None | int:
        return self._num_conflicts_mat

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

        if self._isag_marks is None:
            raise ValueError('Error. No array of snp names to verify')

        num_isag_mark = len(self._isag_marks)
        min_available = np.floor(num_isag_mark * self._min_num_snp)

        min_num_comm_snp = int(num_isag_mark - (2 * (num_isag_mark * 0.05)))

        sample_mark = data.loc[
            data[snp_name_col].isin(self._isag_marks), [descendant, parent]
        ]

        # The number of markers is not 5ok
        desc_n_markers = (sample_mark[descendant] < 5).sum()
        parent_n_markers = (sample_mark[parent] < 5).sum()

        # 1. Checking SNP availability (Call Rate)
        # According to ICAR, the number of markers not 5ok should be more
        # than 95%
        if (desc_n_markers < min_available) and \
                (parent_n_markers < min_available):
            self._status_sing = 'Not Verified'
            self._num_conflicts_sing = None
            return

        # 2. Search for common valid SNPs
        comm_snp_no_missing = sample_mark.replace(5, np.nan).dropna()
        num_comm_markers = len(comm_snp_no_missing)

        if num_comm_markers < min_num_comm_snp:
            self._status_sing = 'Not Checked'
            self._num_conflicts_sing = None
            return

        # 3. Counting conflicts (both homozygotes for different alleles:
        # 0 and 2)
        self._num_conflicts_sing = (abs(
            comm_snp_no_missing[descendant] - comm_snp_no_missing[parent]
        ) == 2).sum()

        # 4. Determination of status under new ICAR rules
        if self._num_conflicts_sing is not None:
            if 0 <= self._num_conflicts_sing <= 2:
                self._status_sing = 'Accept'
            elif 3 <= self._num_conflicts_sing <= 5:
                self._status_sing = 'Doubtful'
            elif self._num_conflicts_sing > 5:
                self._status_sing = 'Excluded'
            else:
                self._status_sing = None

    def check_mating(
            self,
            data: pd.DataFrame,
            descendant: str,
            sire: str,
            dam: str,
            snp_name_col: str
    ) -> None:
        """ Verify mating combination / Trio (Step 2). """
        if self._isag_marks is None:
            raise ValueError('Error. No array of snp names to verify')

        total_snps = len(self._isag_marks)
        min_available = np.floor(total_snps * 0.95)
        # Minimum common SNP for trio: Total - 3 * (Total - Min_Available)
        min_common_trio = total_snps - 3 * (total_snps - min_available)

        sample_data = data[
            data[snp_name_col].isin(self._isag_marks) & \
            ((data[descendant] != 5) & (data[sire] != 5) & (data[dam] != 5))
        ].set_index(snp_name_col)

        # 1. Finding common valid SNPs for the trio
        self._num_common_mat = sample_data.shape[0]

        if self._num_common_mat < min_common_trio:
            self._status_mat = 'Not Checked'
            self._num_conflicts_mat = None
            return

        # 2. Conflict counting (Parents are homozygous for ONE allele,
        # offspring is heterozygous)
        parents_homo_same = (sample_data[sire].isin([0, 2])) & \
                            (sample_data[dam].isin([0, 2])) & \
                            (sample_data[sire] == sample_data[dam])
        progeny_hetero = (sample_data[descendant] == 1)

        conflicts = parents_homo_same & progeny_hetero
        self._num_conflicts_mat = conflicts.sum()

        # 3. Определение статуса согласно ICAR Feb 2025
        if 0 <= self._num_conflicts_mat <= 3:
            self._status_mat = 'Mating Accepted'
        elif 4 <= self._num_conflicts_mat <= 7:
            self._status_mat = 'Mating Doubtful'
        elif self._num_conflicts_mat > 7:
            self._status_mat = 'Mating Excluded'
