#!/usr/bin/env python
# coding: utf-8
__author__ = "Igor Loschinin (igor.loschinin@gmail.com)"

from . import DIR_DATA
from snptools.src.snplib.parentage import (
	Discovery,
	isag_disc
)

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
		assert obj_discovery.status_sing == "Excluded"
		assert obj_discovery.num_conflicts_sing == 77
		assert obj_discovery.perc_conflicts_sing == 14.89

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
		assert obj_discovery.status_sing is None
		assert obj_discovery.num_conflicts_sing is None
		assert obj_discovery.perc_conflicts_sing is None

	def test_search_parent_2(
		self, data: pd.DataFrame, obj_discovery: Discovery
	) -> None:
		"""
		Exception when the number of markers required to confirm paternity is
		less than the established value.
		"""

		assert obj_discovery.search_parent(
			data=data[:-150],
			descendant="BY000041988163",
			parents="EE10512586",
			snp_name_col="SNP_Name"
		) is None
		assert obj_discovery.status_sing == 'Not Discovered'
		assert obj_discovery.num_conflicts_sing == 58
		assert obj_discovery.perc_conflicts_sing == 15.72

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
		assert obj_discovery.status_sing is None
		assert obj_discovery.num_conflicts_sing is None
		assert obj_discovery.perc_conflicts_sing is None

		# For parents
		with pytest.raises(KeyError):
			obj_discovery.search_parent(
				data=data,
				descendant="BY000041988163",
				parents="EE105125864",
				snp_name_col="SNP_Name"
			)
		assert obj_discovery.status_sing is None
		assert obj_discovery.num_conflicts_sing is None
		assert obj_discovery.perc_conflicts_sing is None

	def test_search_parent_4(
		self, data: pd.DataFrame, obj_discovery: Discovery
	) -> None:
		"""
		Test when all snp data is not read - equal to 5.
		"""
		data[["BY000041988163", "EE10512586"]] = 5

		assert obj_discovery.search_parent(
			data=data,
			descendant="BY000041988163",
			parents="EE10512586",
			snp_name_col="SNP_Name"
		) is None
		assert obj_discovery.status_sing == 'Not Discovered'
		assert obj_discovery.num_conflicts_sing == 0
		assert obj_discovery.perc_conflicts_sing is None

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
		assert obj_discovery.status_sing == "Discovered"
		assert obj_discovery.num_conflicts_sing == 0
		assert obj_discovery.perc_conflicts_sing == 0.0

	def test_search_parent_6(
			self, data: pd.DataFrame, obj_discovery: Discovery
	) -> None:
		"""
		Partial match test.
		"""
		data.loc[150:, "EE10512586"] = 1

		obj_discovery.search_parent(
			data=data,
			descendant="BY000041988163",
			parents="EE10512586",
			snp_name_col="SNP_Name"
		)
		assert obj_discovery.status_sing == "Possible"
		assert obj_discovery.num_conflicts_sing == 5
		assert obj_discovery.perc_conflicts_sing == 0.97

	# def test_search_mating_confirmed(
	# 		self, data: pd.DataFrame, obj_discovery: Discovery
	# ) -> None:
	# 	"""
	# 	Test mating with 0 conflicts (Confirmed).
	# 	"""
	# 	data_copy = data.copy()
	# 	snp_name_col = "SNP_Name"
	# 	descendant = "BY000041988163"
	# 	sire = "EE10512586"
	#
	# 	data_copy = data_copy.set_index(snp_name_col)
	# 	desc_vals = data_copy[descendant]
	# 	sire_vals = data_copy[sire]
	#
	# 	# Create a dam with 0 conflicts
	# 	dam_vals = sire_vals.copy()
	# 	mask = (desc_vals == 1) & (sire_vals.isin([0, 2]))
	# 	dam_vals[mask] = 2 - sire_vals[mask]
	#
	# 	data_copy["Dam_0"] = dam_vals
	# 	data_copy = data_copy.reset_index()
	#
	# 	obj_discovery.search_mating(
	# 		data=data_copy,
	# 		descendant=descendant,
	# 		sire=sire,
	# 		dam="Dam_0",
	# 		snp_name_col=snp_name_col
	# 	)
	#
	# 	markers = obj_discovery._Discovery__isag_markers
	# 	num_common, num_conflicts, perc_conflicts = _calculate_step2(
	# 		data_copy, descendant, sire, "Dam_0", snp_name_col, markers
	# 	)
	#
	# 	assert obj_discovery.status_step2 == "Confirmed"
	# 	assert obj_discovery.num_conflicts_step2 == num_conflicts
	# 	assert obj_discovery.perc_conflicts_step2 == perc_conflicts
	#
	# def test_search_mating_possible(
	# 		self, data: pd.DataFrame, obj_discovery: Discovery
	# ) -> None:
	# 	"""
	# 	Test mating with ~1.63% conflicts (Possible).
	# 	"""
	# 	data_copy = data.copy()
	# 	snp_name_col = "SNP_Name"
	# 	descendant = "BY000041988163"
	# 	sire = "EE10512586"
	#
	# 	data_copy = data_copy.set_index(snp_name_col)
	# 	desc_vals = data_copy[descendant]
	# 	sire_vals = data_copy[sire]
	#
	# 	# Start with 0 conflicts dam
	# 	dam_vals = sire_vals.copy()
	# 	mask = (desc_vals == 1) & (sire_vals.isin([0, 2]))
	# 	dam_vals[mask] = 2 - sire_vals[mask]
	#
	# 	# Force 10 conflicts
	# 	conflict_snp_names = data_copy.index[
	# 		(desc_vals == 1) & (sire_vals.isin([0, 2])) & (desc_vals != 5) & (
	# 					sire_vals != 5)]
	# 	dam_vals.loc[conflict_snp_names[:10]] = sire_vals.loc[
	# 		conflict_snp_names[:10]]
	#
	# 	data_copy["Dam_Possible"] = dam_vals
	# 	data_copy = data_copy.reset_index()
	#
	# 	obj_discovery.search_mating(
	# 		data=data_copy,
	# 		descendant=descendant,
	# 		sire=sire,
	# 		dam="Dam_Possible",
	# 		snp_name_col=snp_name_col
	# 	)
	#
	# 	markers = obj_discovery._Discovery__isag_markers
	# 	num_common, num_conflicts, perc_conflicts = _calculate_step2(
	# 		data_copy, descendant, sire, "Dam_Possible", snp_name_col, markers
	# 	)
	#
	# 	assert obj_discovery.status_step2 == "Possible"
	# 	assert obj_discovery.num_conflicts_step2 == num_conflicts
	# 	assert obj_discovery.perc_conflicts_step2 == perc_conflicts
	#
	# def test_search_mating_excluded(
	# 		self, data: pd.DataFrame, obj_discovery: Discovery
	# ) -> None:
	# 	"""
	# 	Test mating with high conflicts (Excluded).
	# 	"""
	# 	data_copy = data.copy()
	# 	snp_name_col = "SNP_Name"
	# 	descendant = "BY000041988163"
	# 	sire = "EE10512586"
	#
	# 	# Dam is identical to sire, which will cause many conflicts
	# 	data_copy["Dam_Excluded"] = data_copy[sire]
	#
	# 	obj_discovery.search_mating(
	# 		data=data_copy,
	# 		descendant=descendant,
	# 		sire=sire,
	# 		dam="Dam_Excluded",
	# 		snp_name_col=snp_name_col
	# 	)
	#
	# 	markers = obj_discovery._Discovery__isag_markers
	# 	num_common, num_conflicts, perc_conflicts = _calculate_step2(
	# 		data_copy, descendant, sire, "Dam_Excluded", snp_name_col, markers
	# 	)
	#
	# 	assert obj_discovery.status_step2 == "Excluded"
	# 	assert obj_discovery.num_conflicts_step2 == num_conflicts
	# 	assert obj_discovery.perc_conflicts_step2 == perc_conflicts
	#
	# def test_search_mating_not_checked(
	# 		self, data: pd.DataFrame, obj_discovery: Discovery
	# ) -> None:
	# 	"""
	# 	Status 'Not Checked' when the number of common markers for trio is less than 400.
	# 	"""
	# 	data_copy = data.copy()
	# 	data_copy["Dam_Fake"] = data_copy["EE10512586"]
	#
	# 	obj_discovery.search_mating(
	# 		data=data_copy[:-250],
	# 		# Slicing to ensure < 400 common SNPs for trio
	# 		descendant="BY000041988163",
	# 		sire="EE10512586",
	# 		dam="Dam_Fake",
	# 		snp_name_col="SNP_Name"
	# 	)
	#
	# 	assert obj_discovery.status_step2 == "Not Checked"
	# 	assert obj_discovery.num_conflicts_step2 == 0
	# 	assert obj_discovery.perc_conflicts_step2 is None
