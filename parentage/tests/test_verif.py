#!/usr/bin/env python
# coding: utf-8

from . import DIR_DATA
from .. import Verification, isag_verif

import pytest
import pandas as pd


@pytest.fixture
def data() -> pd.DataFrame:
	return pd.read_csv(DIR_DATA / "parentage_test_verf.csv", sep=" ")


@pytest.fixture
def obj_verification() -> Verification:
	return Verification(isag_marks=isag_verif().markers)


class TestVerification(object):

	def test_check_on_0(
		self, data: pd.DataFrame, obj_verification: Verification
	) -> None:
		"""
		successfully!
		"""

		assert obj_verification.check_on(
			data=data,
			descendant="BY000041988163",
			parent="EE10512586",
			snp_name_col="SNP_Name"
		) is None
		assert obj_verification.num_conflicts == 31
		assert obj_verification.status == "Excluded"

	def test_check_on_1(self, data: pd.DataFrame) -> None:
		"""
		тест проверяет исключение на отсутствие данных о маркерах для
		верификации
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
		assert obj_verification.status is None
		assert obj_verification.num_conflicts is None

	def test_check_on_2(
		self, data: pd.DataFrame, obj_verification: Verification
	) -> None:
		"""
		исключение на низкий call rate у обоих животных
		"""

		with pytest.raises(
			Exception, match="Calf and parent have low call rate"
		):
			obj_verification.check_on(
				data=data[:-100],
				descendant="BY000041988163",
				parent="EE10512586",
				snp_name_col="SNP_Name"
			)
		assert obj_verification.status is None
		assert obj_verification.num_conflicts is None

	def test_check_on_3(
		self, data: pd.DataFrame, obj_verification: Verification
	) -> None:
		"""
		исключение когда парный call rate ниже порогового значения
		"""

		data.loc[228:, 'BY000041988163'] = 5
		data.loc[239:, 'EE10512586'] = 5

		with pytest.raises(
			Exception, match="Pair call rate is low"
		):
			obj_verification.check_on(
				data=data,
				descendant="BY000041988163",
				parent="EE10512586",
				snp_name_col="SNP_Name"
			)
		assert obj_verification.status is None
		assert obj_verification.num_conflicts is None

	def test_search_parent_4(
		self, data: pd.DataFrame, obj_verification: Verification
	) -> None:
		"""
		тест если передаваемые названия жифотных отсутствуют в датафрейме
		"""

		# For descendant
		with pytest.raises(KeyError):
			obj_verification.check_on(
				data=data,
				descendant="BY00004198816",
				parent="EE10512586",
				snp_name_col="SNP_Name"
			)
		assert obj_verification.status is None
		assert obj_verification.num_conflicts is None

		# For parents
		with pytest.raises(KeyError):
			obj_verification.check_on(
				data=data,
				descendant="BY000041988163",
				parent="EE105125864",
				snp_name_col="SNP_Name"
			)
		assert obj_verification.status is None
		assert obj_verification.num_conflicts is None

	def test_search_parent_5(
		self, data: pd.DataFrame, obj_verification: Verification
	) -> None:
		"""
		тест когда все данные снп не прочитанны - равны 5
		"""
		data[["BY000041988163", "EE10512586"]] = 5

		with pytest.raises(
			Exception, match="Calf and parent have low call rate"
		):
			obj_verification.check_on(
				data=data,
				descendant="BY000041988163",
				parent="EE10512586",
				snp_name_col="SNP_Name"
			)
		assert obj_verification.status is None
		assert obj_verification.num_conflicts is None

	def test_search_parent_6(
			self, data: pd.DataFrame, obj_verification: Verification
	) -> None:
		"""
		тест когда полное совпадение
		"""
		data[["BY000041988163", "EE10512586"]] = 2

		obj_verification.check_on(
			data=data,
			descendant="BY000041988163",
			parent="EE10512586",
			snp_name_col="SNP_Name"
		)
		assert obj_verification.status == "Accept"
		assert obj_verification.num_conflicts == 0
