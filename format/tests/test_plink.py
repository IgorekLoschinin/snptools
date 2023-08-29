#!/usr/bin/env python
# coding: utf-8

from . import DIR_FILES
from .. import (
	make_ped,
	make_map,
	make_fam,
	make_lgen
)

import pytest
import pandas as pd


@pytest.fixture
def data_ped() -> pd.DataFrame:
	"""
	    SAMPLE_ID  SNP
	0        1100  025
	1        1101  022
	2        1102  052
	3        1103  022

	:return: - Data in formate snp
	"""
	return pd.read_pickle(DIR_FILES / "fplink/ped/file.pl")


class TestPlinkFormatMap(object):

	def test_map(self) -> None:
		data = pd.read_csv(
			DIR_FILES / "fplink/map/file_bovinesnp50.csv",

		)

		g = make_map(data)


class TestPlinkFormatPed(object):

	def test_ped(self, data_ped: pd.DataFrame) -> None:

		f = make_ped(
			data=data_ped,
			sid_col="SAMPLE_ID",
			fid_col="SAMPLE_ID",
			snp_col="SNP"
		)

		print()


class TestPlinkFormatFam(object):

	def test_fam(self) -> None:
		...


class TestPlinkFormatLgen(object):

	def test_lgen(self) -> None:
		...
