#!/usr/bin/env python
# coding: utf-8

from . import DIR_DATA
from .. import Discovery, isag_disc

import pytest
import pandas as pd


@pytest.fixture
def data() -> pd.DataFrame:
	return pd.read_csv(DIR_DATA / "parentage_test_disc.csv", sep=" ")


@pytest.fixture
def obj_discovery() -> Discovery:
	return Discovery(isag_markers=isag_disc().markers)


class TestDiscovery(object):

	def test_search_parent_0(
		self, data: pd.DataFrame, obj_discovery: Discovery
	) -> None:
		"""
		successfully!
		"""

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
		когда выбрасывается исключение на отсутствие данных с маркерами исаг
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
		исключение когда количество необходимых маркеров для подтверждения
		отцовства меньше установленного значения
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
		тест если передаваемые названия жифотных отсутствуют в датафрейме
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
		тест когда все данные снп не прочитанны - равны 5
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
		тест когда полное совпадение
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
		тест на неполное совпание
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
