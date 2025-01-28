#!/usr/bin/env python
# coding: utf-8
__author__ = "Igor Loschinin (igor.loschinin@gmail.com)"

from . import DIR_DATA
from snptools.src.snplib.parentage import Verification, isag_verif

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
		assert obj_verification.num_conflicts == 31
		assert obj_verification.status == "Excluded"

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
		assert obj_verification.status is None
		assert obj_verification.num_conflicts is None

	def test_check_on_2(
		self, data: pd.DataFrame, obj_verification: Verification
	) -> None:
		"""
		Exception for low call rate in both animals.
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
		Exception when paired call rate is below threshold.
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
		Test when all snp data is not read - equal to 5
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
		Test when there is a complete match
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
