#!/usr/bin/env python
# coding: utf-8
__author__ = "Igor Loschinin (igor.loschinin@gmail.com)"

from . import DIR_FILES
from snptools.src.snplib.format import make_ped

import pytest
import pandas as pd


@pytest.fixture
def data_ped(request) -> pd.DataFrame | None:
	return pd.read_pickle(DIR_FILES / f"fplink/ped/{request.param}")


class TestPlinkFormatPed(object):

	@pytest.mark.parametrize("data_ped", ["file.pl"], indirect=True)
	def test_ped_true(self, data_ped: pd.DataFrame) -> None:
		assert not make_ped(
			data_ped,
			"SAMPLE_ID",
			"SNP",
			fid_col="SAMPLE_ID"
		).empty

		assert not make_ped(
			data_ped,
			"SAMPLE_ID",
			"SNP"
		).empty

	def test_ped_empty(self) -> None:
		assert make_ped(
			pd.DataFrame(columns=["SAMPLE_ID", "SNP"]),
			"SAMPLE_ID",
			"SNP"
		).empty

		assert make_ped(
			pd.DataFrame(columns=["SAMPLE_ID", "SNP"]),
			"SAMPLE_ID",
			"SNP",
			fid_col="SAMPLE_ID"
		).empty

	@pytest.mark.parametrize("data_ped", ["file.pl"], indirect=True)
	def test_ped_raise_columns(self, data_ped: pd.DataFrame) -> None:
		# SID_COL
		with pytest.raises(
				KeyError, match="Data has not in name columns!"
		):
			make_ped(
				data=data_ped,
				sid_col="SAMPLE_ID1",
				fid_col="SAMPLE_ID",
				snp_col="SNP"
			)

		# SNP_COL
		with pytest.raises(
				KeyError, match="Data has not in name columns!"
		):
			make_ped(
				data_ped,
				"SAMPLE_ID",
				"SNP1",
				fid_col="SAMPLE_ID"
			)

		# FID_COL
		with pytest.raises(
				KeyError, match="Data has not in name columns SAMPLE_ID1!"
		):
			make_ped(
				data_ped,
				"SAMPLE_ID",
				"SNP",
				fid_col="SAMPLE_ID1"
			)

	@pytest.mark.parametrize("data_ped", ["file2.pl"], indirect=True)
	def test_ped_raises_underscope_sid(self, data_ped: pd.DataFrame) -> None:

		# SID_COL
		with pytest.raises(
				Exception,
				match="Replace in 'Sample ID' columns '_' on another a simbols"
		):
			res = make_ped(
				data_ped,
				"SAMPLE_ID",
				"SNP"
			)

	@pytest.mark.parametrize("data_ped", ["file3.pl"], indirect=True)
	def test_ped_raises_underscope_fid(self, data_ped: pd.DataFrame) -> None:

		# FID_COL
		with pytest.raises(
				Exception,
				match="Replace in 'Family ID' columns '_' on another a simbols"
		):
			res = make_ped(
				data_ped,
				"SAMPLE_ID",
				"SNP",
				fid_col="FAMILY_ID"
			)

	@pytest.mark.parametrize("data_ped", ["file4.pl"], indirect=True)
	def test_ped_check_data(self, data_ped: pd.DataFrame) -> None:
		res = make_ped(
			data_ped,
			"SAMPLE_ID",
			"SNP",
			fid_col="FAMILY_ID",
			father_col="father",
			mother_col="mother",
			sex_col="sex"
		)

		res2 = make_ped(
			data_ped,
			"SAMPLE_ID",
			"SNP",
			fid_col="FAMILY_ID",
		)

		assert all(res.father.values == list('1234'))
		assert all(res.mother.values == list('5678'))
		assert all(res.sex.values == list('1210'))
		assert all(res2.father.values == list('0000'))
		assert all(res2.mother.values == list('0000'))
		assert all(res2.sex.values == list('0000'))
