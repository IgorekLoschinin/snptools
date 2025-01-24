#!/usr/bin/env python
# coding: utf-8
__author__ = "Igor Loschinin (igor.loschinin@gmail.com)"

from . import DIR_DATA
from snptools.src.snplib.parentage import Discovery, isag_disc

import pytest
import pandas as pd


@pytest.fixture
def data() -> pd.DataFrame:
	return pd.read_csv(DIR_DATA / "parentage_test_disc.csv", sep=" ")


@pytest.fixture
def obj_discovery() -> Discovery:
	return Discovery(isag_markers=isag_disc().markers)


class TestDiscovery(object):

	def test_search_parent_successfully(
		self, data: pd.DataFrame, obj_discovery: Discovery
	) -> None:

		assert obj_discovery.search_parent(
			data=data,
			descendant="BY000041988163",
			parents="EE10512586",
			snp_name_col="SNP_Name"
		) is None
		assert obj_discovery.num_conflicts == 77
		assert obj_discovery.status == "Excluded"
		assert obj_discovery.perc_conflicts == 14.86

	def test_search_parent_1(self, data: pd.DataFrame) -> None:
		"""
		An exception is thrown for the absence of data with isag markers
		"""
		obj_discovery = Discovery()

		with pytest.raises(
			ValueError, match="Error. No array of snp names to verify"
		):
			obj_discovery.search_parent(
				data=data,
				descendant="BY000041988163",
				parents="EE10512586",
				snp_name_col="SNP_Name"
			)
		assert obj_discovery.status is None
		assert obj_discovery.num_conflicts is None
		assert obj_discovery.perc_conflicts is None

	def test_search_parent_2(
		self, data: pd.DataFrame, obj_discovery: Discovery
	) -> None:
		"""
		Exception when the number of markers required to confirm paternity is
		less than the established value.
		"""

		with pytest.raises(
			Exception, match="Calf call rate is low."
		):
			obj_discovery.search_parent(
				data=data[:-100],
				descendant="BY000041988163",
				parents="EE10512586",
				snp_name_col="SNP_Name"
			)
		assert obj_discovery.status is None
		assert obj_discovery.num_conflicts is None
		assert obj_discovery.perc_conflicts is None

	def test_search_parent_3(
		self, data: pd.DataFrame, obj_discovery: Discovery
	) -> None:
		"""
		Test if the transmitted animal names are not in the dataframe.
		"""

		# For descendant
		with pytest.raises(KeyError):
			obj_discovery.search_parent(
				data=data,
				descendant="BY00004198816",
				parents="EE10512586",
				snp_name_col="SNP_Name"
			)
		assert obj_discovery.status is None
		assert obj_discovery.num_conflicts is None
		assert obj_discovery.perc_conflicts is None

		# For parents
		with pytest.raises(KeyError):
			obj_discovery.search_parent(
				data=data,
				descendant="BY000041988163",
				parents="EE105125864",
				snp_name_col="SNP_Name"
			)
		assert obj_discovery.status is None
		assert obj_discovery.num_conflicts is None
		assert obj_discovery.perc_conflicts is None

	def test_search_parent_4(
		self, data: pd.DataFrame, obj_discovery: Discovery
	) -> None:
		"""
		Test when all snp data is not read - equal to 5.
		"""
		data[["BY000041988163", "EE10512586"]] = 5

		with pytest.raises(
			Exception, match="Calf call rate is low."
		):
			obj_discovery.search_parent(
				data=data,
				descendant="BY000041988163",
				parents="EE10512586",
				snp_name_col="SNP_Name"
			)
		assert obj_discovery.status is None
		assert obj_discovery.num_conflicts is None
		assert obj_discovery.perc_conflicts is None

	def test_search_parent_5(
			self, data: pd.DataFrame, obj_discovery: Discovery
	) -> None:
		"""
		Test when there is a complete match.
		"""
		data[["BY000041988163", "EE10512586"]] = 2

		obj_discovery.search_parent(
			data=data,
			descendant="BY000041988163",
			parents="EE10512586",
			snp_name_col="SNP_Name"
		)
		assert obj_discovery.status == "Discovered"
		assert obj_discovery.num_conflicts == 0
		assert obj_discovery.perc_conflicts == 0.0

	def test_search_parent_6(
			self, data: pd.DataFrame, obj_discovery: Discovery
	) -> None:
		"""
		Partial match test.
		"""
		data.loc[202:, "EE10512586"] = 1

		obj_discovery.search_parent(
			data=data,
			descendant="BY000041988163",
			parents="EE10512586",
			snp_name_col="SNP_Name"
		)
		assert obj_discovery.status == "Doubtful"
		assert obj_discovery.num_conflicts == 14
		assert obj_discovery.perc_conflicts == 2.70
