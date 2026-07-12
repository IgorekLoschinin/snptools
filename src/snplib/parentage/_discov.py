#!/usr/bin/env python
# coding: utf-8
__author__ = "Igor Loschinin (igor.loschinin@gmail.com)"

import pandas as pd

"""
Search for paternity according to ICAR recommendations
https://www.icar.org/wp-content/uploads/documents/Guidelines-for-Parentage-Verification-and-Parentage-Discovery-Based-on-SNP-Genotypes-Feb.-2025.pdf
"""


class Discovery(object):
	""" Search for paternity according to ICAR recommendations (Revised Feb 2025)

	:argument isag_markers: Fixed sample of markers to confirm paternity.
	"""

	def __init__(
			self, isag_markers: pd.Series | list | set | None = None
	) -> None:
		self._isag_markers = isag_markers

		# Attributes for Step 1 (Find one parent)
		self.__num_common_sing = None
		self.__num_conflicts_sing = None
		self.__perc_conflicts_sing = None
		self.__status_sing = None

		# Attributes for Step 2 (Checking Trio/Mating)
		self.__num_common_mat = None
		self.__num_conflicts_mat = None
		self.__perc_conflicts_mat = None
		self.__status_mat = None

	@property
	def status_sing(self) -> pd.Series | str | None:
		""" The status of each parent discovered (Step 1). """
		return self.__status_sing

	@property
	def status_mating(self) -> str | None:
		""" The status of the mating combination (Step 2). """
		return self.__status_mat

	@property
	def num_conflicts_sing(self) -> pd.Series | int | None:
		return self.__num_conflicts_sing

	@property
	def perc_conflicts_sing(self) -> pd.Series | float | None:
		return self.__perc_conflicts_sing

	def search_parent(
			self,
			data: pd.DataFrame,
			descendant: str,
			parents: str,
			snp_name_col: str
	) -> None:
		""" Search for paternity (Step 1).

		:param data: SNP data for descendant and parent(s).
		:param descendant: Column name of the descendant in the data.
		:param parents: Column name of column names of the parents in
		the data.
		:param snp_name_col: SNP column name in data.
		"""
		if self._isag_markers is None:
			raise ValueError("Error. No array of snp names to verify")

		# 1. Sampling isag markers and clearing data from gaps (5)
		sample_by_markers = data.loc[
			data[snp_name_col].isin(self._isag_markers) & \
			((data[descendant] != 5) & (data[parents] != 5))
		].set_index(snp_name_col)

		self.__num_common_sing = sample_by_markers.shape[0]

		# 2. We count conflicts (both are homozygous for different alleles,
		# for example, 0 and 2)
		self.__num_conflicts_sing = (
			(sample_by_markers[parents] - sample_by_markers[descendant]).abs() == 2
		).sum()

		# 3. We calculate the percentage of conflicts
		if self.__num_common_sing != 0:
			self.__perc_conflicts_sing = (
				(self.__num_conflicts_sing / self.__num_common_sing) * 100
			).round(2)

		# 4. We assign statuses according to the rules of ICAR
		if self.__num_common_sing < 400:
			self.__status_sing = 'Not Discovered'
		else:
			if self.__perc_conflicts_sing is not None:
				if 0 <= self.__perc_conflicts_sing <= 0.5:
					self.__status_sing = 'Discovered'
				elif 0.5 < self.__perc_conflicts_sing <= 2.0:
					self.__status_sing = 'Possible'
				elif self.__perc_conflicts_sing > 2:
					self.__status_sing = 'Excluded'
				else:
					self.__status_sing = None

	def search_mating(
			self,
			data: pd.DataFrame,
			descendant: str,
			sire: str,
			dam: str,
			snp_name_col: str
	) -> None:
		""" Verify mating combination / Trio (Step 2).

		:param data: SNP data for descendant, sire, and dam.
		:param descendant: Column name of the descendant.
		:param sire: Column name of the sire.
		:param dam: Column name of the dam.
		:param snp_name_col: SNP column name in data.
		"""
		if self._isag_markers is None:
			raise ValueError("Error. No array of snp names to verify")

		sample_by_markers = data.loc[
			data[snp_name_col].isin(self._isag_markers)
		].set_index(snp_name_col)

		desc_vals = sample_by_markers[descendant]
		sire_vals = sample_by_markers[sire]
		dam_vals = sample_by_markers[dam]

		# Общие валидные SNP для трио
		valid_desc = desc_vals != 5
		valid_sire = sire_vals != 5
		valid_dam = dam_vals != 5
		common_valid = valid_desc & valid_sire & valid_dam

		self.__num_common_mat = common_valid.sum()

		if self.__num_common_mat < 400:
			self.__num_conflicts_mat = 0
			self.__perc_conflicts_mat = None
			self.__status_mat = 'Not Checked'
			return

		# Конфликты для трио: оба родителя гомозиготы по ОДНОМУ аллелю,
		# а потомок гетерозиготен
		parents_homozygous_same = (sire_vals.isin([0, 2])) & (
			dam_vals.isin([0, 2])) & (sire_vals == dam_vals)
		progeny_heterozygous = (desc_vals == 1)

		conflicts = parents_homozygous_same & progeny_heterozygous & common_valid
		self.__num_conflicts_mat = conflicts.sum()

		perc_conflicts = (self.__num_conflicts_mat / self.__num_common_mat) * 100
		self.__perc_conflicts_mat = round(perc_conflicts, 2)

		# Статусы для Шага 2 согласно ICAR Feb 2025
		if 0 < perc_conflicts <= 1.0:
			self.__status_mat = 'Confirmed'
		elif 1 < perc_conflicts <= 4.0:
			self.__status_mat = 'Possible'
		elif perc_conflicts > 4.0:
			self.__status_mat = 'Excluded'
		else:
			self.__status_mat = None
