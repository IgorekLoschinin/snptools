#!/usr/bin/env python
# coding: utf-8
__author__ = "Igor Loschinin (igor.loschinin@gmail.com)"

import numpy as np

from . import DIR_DATA
from snptools.src.snplib.parentage import (
	Discovery,
	isag_disc
)

import pytest
import pandas as pd


@pytest.fixture
def obj_discovery() -> Discovery:
	return Discovery(isag_markers=isag_disc().markers)


@pytest.fixture
def data() -> pd.DataFrame:
	return pd.read_csv(DIR_DATA / "parentage_test_disc.csv", sep=" ")


def generate_dam_vector(
		base_df: pd.DataFrame,
		n_conflicts: int,
		descendant_col: str,
		sire_col: str,
		n_missing: int = 0,
		seed: int = 42,
) -> pd.Series:
	""" Generates vector results for ISAG markers.

	:param descendant_col: Descendant column name
	:param sire_col: Sibling column name
	:param base_df: DataFrame with the columns SNP_Name, descendant, and sire
	:param n_conflicts: Number of constraints to create
	:param n_missing: Number of scores equivalent to 5 (missing)
	:param seed: Seeds for reproducibility
	:return: pd.Series with the index SNP_Name
	"""
	rng = np.random.default_rng(seed)
	df = base_df.copy()

	# Candidates for conflicts: descendant==1 and sire ∈ {0,2}
	conflict_candidates = df.index[
		(df[descendant_col] == 1) & (df[sire_col].isin([0, 2]))
		].tolist()

	if n_conflicts > len(conflict_candidates):
		raise ValueError(
			f"Cannot create {n_conflicts} conflicts: "
			f"only {len(conflict_candidates)} candidates available"
		)

	conflict_idx = rng.choice(
		conflict_candidates, size=n_conflicts, replace=False
	)
	non_conflict_idx = np.setdiff1d(df.index, conflict_idx)

	dam = np.empty(len(df), dtype=int)

	# 1) In conflict positions: dam = sire (both homozygotes, descendant 1)
	dam[conflict_idx] = df.loc[conflict_idx, sire_col].values

	# 2) In other positions: valid values without conflict
	for i in non_conflict_idx:
		desc = df.loc[i, descendant_col]
		sire = df.loc[i, sire_col]

		# All valid combinations (sire, dam, desc) without conflict
		valid_dams = [
			d
			for d in (0, 1, 2)
			if not (sire == d and sire in (0, 2) and desc == 1)
		]
		dam[i] = rng.choice(valid_dams)

	# 3) Missing (5)
	if n_missing > 0:
		missing_idx = rng.choice(df.index, size=n_missing, replace=False)
		dam[missing_idx] = 5

	return pd.Series(dam, index=df["SNP_Name"], name="dam")


@pytest.fixture
def data_mat(request) -> pd.DataFrame:

	desc_col = "BY000041988163"
	sire_col = "EE10512586"

	n_conflicts, n_missing, seed = request.param

	base = pd.read_csv(DIR_DATA / "parentage_test_disc.csv", sep=" ")

	dam_vec = generate_dam_vector(
		base,
		descendant_col=desc_col,
		sire_col=sire_col,
		n_conflicts=n_conflicts,
		n_missing=n_missing,
		seed=seed
	).reset_index()

	trio = base[[desc_col, sire_col]].copy()
	trio[["SNP_Name", "dam"]] = dam_vec.values

	return trio[["SNP_Name", "dam", desc_col, sire_col]]


class TestDiscovery(object):

	snp_name_col = "SNP_Name"
	descendant = "BY000041988163"
	sire = "EE10512586"
	dam = "dam"

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

	@pytest.mark.parametrize("data_mat", [(0, 200, 1)], indirect=True)
	def test_search_mating_not_checked(
			self, data_mat: pd.DataFrame, obj_discovery: Discovery
	) -> None:

		obj_discovery.search_mating(
			data=data_mat,
			descendant=self.descendant,
			sire=self.sire,
			dam=self.dam,
			snp_name_col=self.snp_name_col
		)

		assert obj_discovery.status_mating == "Not Checked"
		assert obj_discovery.num_conflicts_mating == 0
		assert obj_discovery.perc_conflicts_mating is None

	@pytest.mark.parametrize("data_mat", [(0, 0, 2)], indirect=True)
	def test_search_mating_none_status(
			self, data_mat: pd.DataFrame, obj_discovery: Discovery
	) -> None:

		obj_discovery.search_mating(
			data=data_mat,
			descendant=self.descendant,
			sire=self.sire,
			dam=self.dam,
			snp_name_col=self.snp_name_col
		)

		assert obj_discovery.status_mating is None
		assert obj_discovery.num_conflicts_mating == 0
		assert obj_discovery.perc_conflicts_mating == 0

	@pytest.mark.parametrize("data_mat", [(3, 0, 3)], indirect=True)
	def test_search_mating_confirmed(
			self, data_mat: pd.DataFrame, obj_discovery: Discovery
	) -> None:

		obj_discovery.search_mating(
			data=data_mat,
			descendant=self.descendant,
			sire=self.sire,
			dam=self.dam,
			snp_name_col=self.snp_name_col
		)

		assert obj_discovery.status_mating == 'Confirmed'
		assert obj_discovery.num_conflicts_mating == 2
		assert obj_discovery.perc_conflicts_mating == 0.39

	@pytest.mark.parametrize("data_mat", [(10, 0, 4)], indirect=True)
	def test_search_mating_possible(
			self, data_mat: pd.DataFrame, obj_discovery: Discovery
	) -> None:

		obj_discovery.search_mating(
			data=data_mat,
			descendant=self.descendant,
			sire=self.sire,
			dam=self.dam,
			snp_name_col=self.snp_name_col
		)

		assert obj_discovery.status_mating == 'Possible'
		assert obj_discovery.num_conflicts_mating == 9
		assert obj_discovery.perc_conflicts_mating == 1.74

	@pytest.mark.parametrize("data_mat", [(30, 0, 5)], indirect=True)
	def test_search_mating_excluded(
			self, data_mat: pd.DataFrame, obj_discovery: Discovery
	) -> None:

		obj_discovery.search_mating(
			data=data_mat,
			descendant=self.descendant,
			sire=self.sire,
			dam=self.dam,
			snp_name_col=self.snp_name_col
		)

		assert obj_discovery.status_mating == 'Excluded'
		assert obj_discovery.num_conflicts_mating == 24
		assert obj_discovery.perc_conflicts_mating == 4.64
