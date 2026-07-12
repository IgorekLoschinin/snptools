#!/usr/bin/env python
# coding: utf-8
__author__ = "Igor Loschinin (igor.loschinin@gmail.com)"

from . import DIR_DATA
from snptools.src.snplib.parentage import (
	Verification,
	isag_verif
)

import pytest
import pandas as pd


@pytest.fixture
def data() -> pd.DataFrame:
	return pd.read_csv(DIR_DATA / "parentage_test_verf.csv", sep=" ")


@pytest.fixture
def obj_verification() -> Verification:
	return Verification(isag_marks=isag_verif().markers)


class TestVerification(object):

	def test_check_on_successfully(
		self, data: pd.DataFrame, obj_verification: Verification
	) -> None:

		assert obj_verification.check_on(
			data=data,
			descendant="BY000041988163",
			parent="EE10512586",
			snp_name_col="SNP_Name"
		) is None
		assert obj_verification.num_conflicts_sing == 31
		assert obj_verification.status_sing == "Excluded"

	def test_check_on_1(self, data: pd.DataFrame) -> None:
		"""
		The test checks the exception for missing token data for verification.
		"""
		obj_verification = Verification()

		with pytest.raises(
			ValueError, match="Error. No array of snp names to verify"
		):
			obj_verification.check_on(
				data=data,
				descendant="BY000041988163",
				parent="EE10512586",
				snp_name_col="SNP_Name"
			)
		assert obj_verification.status_sing is None
		assert obj_verification.num_conflicts_sing is None

	def test_check_on_2(
		self, data: pd.DataFrame, obj_verification: Verification
	) -> None:
		"""
		Exception for low call rate in both animals.
		"""

		assert obj_verification.check_on(
			data=data[:-100],
			descendant="BY000041988163",
			parent="EE10512586",
			snp_name_col="SNP_Name"
		) is None
		assert obj_verification.status_sing == 'Not Verified'
		assert obj_verification.num_conflicts_sing is None

	def test_check_on_3(
		self, data: pd.DataFrame, obj_verification: Verification
	) -> None:
		"""
		Exception when paired call rate is below threshold.
		"""

		data.loc[228:, 'BY000041988163'] = 5
		data.loc[239:, 'EE10512586'] = 5

		assert obj_verification.check_on(
			data=data,
			descendant="BY000041988163",
			parent="EE10512586",
			snp_name_col="SNP_Name"
		) is None
		assert obj_verification.status_sing == 'Not Checked'
		assert obj_verification.num_conflicts_sing is None

	def test_search_parent_4(
		self, data: pd.DataFrame, obj_verification: Verification
	) -> None:
		"""
		Test if the transmitted animal names are not in the dataframe.
		"""

		# For descendant
		with pytest.raises(KeyError):
			obj_verification.check_on(
				data=data,
				descendant="BY00004198816",
				parent="EE10512586",
				snp_name_col="SNP_Name"
			)
		assert obj_verification.status_sing is None
		assert obj_verification.num_conflicts_sing is None

		# For parents
		with pytest.raises(KeyError):
			obj_verification.check_on(
				data=data,
				descendant="BY000041988163",
				parent="EE105125864",
				snp_name_col="SNP_Name"
			)
		assert obj_verification.status_sing is None
		assert obj_verification.num_conflicts_sing is None

	def test_search_parent_5(
		self, data: pd.DataFrame, obj_verification: Verification
	) -> None:
		"""
		Test when all snp data is not read - equal to 5
		"""
		data[["BY000041988163", "EE10512586"]] = 5

		obj_verification.check_on(
			data=data,
			descendant="BY000041988163",
			parent="EE10512586",
			snp_name_col="SNP_Name"
		)
		assert obj_verification.status_sing == 'Not Verified'
		assert obj_verification.num_conflicts_sing is None

	def test_search_parent_6(
			self, data: pd.DataFrame, obj_verification: Verification
	) -> None:
		"""
		Test when there is a complete match
		"""
		data[["BY000041988163", "EE10512586"]] = 2

		obj_verification.check_on(
			data=data,
			descendant="BY000041988163",
			parent="EE10512586",
			snp_name_col="SNP_Name"
		)
		assert obj_verification.status_sing == "Accept"
		assert obj_verification.num_conflicts_sing == 0

	def test_search_parent_7(
			self, data: pd.DataFrame, obj_verification: Verification
	) -> None:
		"""
		Test when there is a Doubtful
		"""
		data.loc[70:, "EE10512586"] = 1

		obj_verification.check_on(
			data=data,
			descendant="BY000041988163",
			parent="EE10512586",
			snp_name_col="SNP_Name"
		)
		assert obj_verification.status_sing == "Doubtful"
		assert obj_verification.num_conflicts_sing == 5

	# # ==================== STEP 2 TESTS ====================
	#
	# def test_check_mating_accepted(self, data: pd.DataFrame,
	# 							   obj_verification: Verification) -> None:
	# 	"""Test 0-3 conflicts -> Mating Accepted."""
	# 	data_copy = data.copy()
	# 	snp_col = "SNP_Name"
	# 	desc = "BY000041988163"
	# 	sire = "EE10512586"
	#
	# 	# Создаем мать с 0 конфликтами
	# 	sample = data_copy.set_index(snp_col)
	# 	desc_vals = sample[desc]
	# 	sire_vals = sample[sire]
	# 	dam_vals = sire_vals.copy()
	#
	# 	# Исправляем потенциальные конфликты (оба родителя гомозиготны, потомок гетерозиготен)
	# 	mask = (desc_vals == 1) & (sire_vals.isin([0, 2])) & (
	# 				desc_vals != 5) & (sire_vals != 5)
	# 	dam_vals[mask] = 2 - sire_vals[mask]
	#
	# 	data_copy["Dam_OK"] = dam_vals.values
	# 	data_copy = data_copy.reset_index()
	#
	# 	obj_verification.check_mating(data_copy, desc, sire, "Dam_OK", snp_col)
	# 	assert obj_verification.status_step2 == 'Mating Accepted'
	# 	assert obj_verification.num_conflicts_step2 <= 3
	#
	# def test_check_mating_doubtful(self, data: pd.DataFrame,
	# 							   obj_verification: Verification) -> None:
	# 	"""Test 4-7 conflicts -> Mating Doubtful."""
	# 	data_copy = data.copy()
	# 	snp_col = "SNP_Name"
	# 	desc = "BY000041988163"
	# 	sire = "EE10512586"
	#
	# 	sample = data_copy.set_index(snp_col)
	# 	desc_vals = sample[desc]
	# 	sire_vals = sample[sire]
	# 	dam_vals = sire_vals.copy()
	#
	# 	# Сначала убираем все конфликты
	# 	mask = (desc_vals == 1) & (sire_vals.isin([0, 2])) & (
	# 				desc_vals != 5) & (sire_vals != 5)
	# 	dam_vals[mask] = 2 - sire_vals[mask]
	#
	# 	# Теперь искусственно создаем ровно 5 конфликтов
	# 	conflict_indices = data_copy.index[
	# 						   (desc_vals == 1) & (sire_vals.isin([0, 2])) & (
	# 									   desc_vals != 5) & (sire_vals != 5)][
	# 					   :5]
	# 	dam_vals.iloc[conflict_indices] = sire_vals.iloc[
	# 		conflict_indices]  # Делаем родителей одинаковыми гомозиготами
	#
	# 	data_copy["Dam_Doubt"] = dam_vals.values
	# 	data_copy = data_copy.reset_index()
	#
	# 	obj_verification.check_mating(data_copy, desc, sire, "Dam_Doubt",
	# 								  snp_col)
	# 	assert obj_verification.status_step2 == 'Mating Doubtful'
	# 	assert obj_verification.num_conflicts_step2 == 5
	#
	# def test_check_mating_excluded(self, data: pd.DataFrame,
	# 							   obj_verification: Verification) -> None:
	# 	"""Test >7 conflicts -> Mating Excluded."""
	# 	data_copy = data.copy()
	# 	snp_col = "SNP_Name"
	# 	desc = "BY000041988163"
	# 	sire = "EE10512586"
	#
	# 	# Мать идентична отцу, что гарантированно создаст много конфликтов там, где потомок гетерозиготен
	# 	data_copy["Dam_Excl"] = data_copy[sire]
	#
	# 	obj_verification.check_mating(data_copy, desc, sire, "Dam_Excl",
	# 								  snp_col)
	#
	# 	markers = obj_verification._Verification__isag_marks
	# 	num_conflicts, _, status = _calculate_step2(data_copy, desc, sire,
	# 												"Dam_Excl", snp_col,
	# 												markers)
	#
	# 	assert obj_verification.status_step2 == status
	# 	assert obj_verification.num_conflicts_step2 == num_conflicts
	# 	assert obj_verification.status_step2 == 'Mating Excluded'
	#
	# def test_check_mating_not_checked(self, data: pd.DataFrame,
	# 								  obj_verification: Verification) -> None:
	# 	"""Status 'Not Checked' when trio common SNPs < 165."""
	# 	data_copy = data.copy()
	# 	data_copy["Dam_Fake"] = data_copy["EE10512586"]
	#
	# 	# Обрезаем данные, чтобы общих SNP для трио точно не хватило
	# 	obj_verification.check_mating(data_copy[:-150], "BY000041988163",
	# 								  "EE10512586", "Dam_Fake", "SNP_Name")
	#
	# 	assert obj_verification.status_step2 == 'Not Checked'
	# 	assert obj_verification.num_conflicts_step2 is None
	#