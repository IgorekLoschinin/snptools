#!/usr/bin/env python
# coding: utf-8
__author__ = "Igor Loschinin (igor.loschinin@gmail.com)"

import numpy as np
import pandas as pd
import pytest

from snptools.src.snplib.parentage import (
	Verification,
	isag_verif
)
from . import DIR_DATA


@pytest.fixture
def data() -> pd.DataFrame:
	return pd.read_csv(DIR_DATA / "parentage_test_verf.csv", sep=" ")


@pytest.fixture
def obj_verification() -> Verification:
	return Verification(isag_marks=isag_verif().markers)


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

	base = pd.read_csv(DIR_DATA / "parentage_test_verf.csv", sep=" ")

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


class TestVerification(object):

	snp_name_col = "SNP_Name"
	descendant = "BY000041988163"
	sire = "EE10512586"
	dam = "dam"

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

	# ==================== STEP 2 TESTS ====================

	@pytest.mark.parametrize("data_mat", [(0, 50, 1)], indirect=True)
	def test_search_mating_not_checked(
			self, data_mat: pd.DataFrame, obj_verification: Verification
	) -> None:
		obj_verification.check_mating(
			data=data_mat,
			descendant=self.descendant,
			sire=self.sire,
			dam=self.dam,
			snp_name_col=self.snp_name_col
		)

		assert obj_verification.status_mat == 'Not Checked'
		assert obj_verification.num_conflicts_mat is None

	@pytest.mark.parametrize("data_mat", [(3, 0, 3)], indirect=True)
	def test_search_mating_confirmed(
			self, data_mat: pd.DataFrame, obj_verification: Verification
	) -> None:
		obj_verification.check_mating(
			data=data_mat,
			descendant=self.descendant,
			sire=self.sire,
			dam=self.dam,
			snp_name_col=self.snp_name_col
		)

		assert obj_verification.status_mat == 'Mating Accepted'
		assert obj_verification.num_conflicts_mat == 1

	@pytest.mark.parametrize("data_mat", [(6, 0, 4)], indirect=True)
	def test_search_mating_doubtful(
			self, data_mat: pd.DataFrame, obj_verification: Verification
	) -> None:
		obj_verification.check_mating(
			data=data_mat,
			descendant=self.descendant,
			sire=self.sire,
			dam=self.dam,
			snp_name_col=self.snp_name_col
		)

		assert obj_verification.status_mat == 'Mating Doubtful'
		assert obj_verification.num_conflicts_mat == 6

	@pytest.mark.parametrize("data_mat", [(30, 0, 5)], indirect=True)
	def test_search_mating_excluded(
			self, data_mat: pd.DataFrame, obj_verification: Verification
	) -> None:
		obj_verification.check_mating(
			data=data_mat,
			descendant=self.descendant,
			sire=self.sire,
			dam=self.dam,
			snp_name_col=self.snp_name_col
		)

		assert obj_verification.status_mat == 'Mating Excluded'
		assert obj_verification.num_conflicts_mat == 21
