#!/usr/bin/env python
# coding: utf-8
__author__ = "Igor Loschinin (igor.loschinin@gmail.com)"

from . import DIR_FILES
from snptools.src.snplib.format import make_fam

import pytest
import pandas as pd


@pytest.fixture
def data_fam(request) -> pd.DataFrame | None:
	return pd.read_pickle(DIR_FILES / f"fplink/fam/{request.param}")


class TestPlinkFormatPed(object):

	@pytest.mark.parametrize("data_fam", ["file.pl"], indirect=True)
	def test_fam_true(self, data_fam: pd.DataFrame) -> None:
		assert not make_fam(
			data_fam,
			"SAMPLE_ID",
			"SAMPLE_ID"
		).empty

		assert not make_fam(
			data_fam,
			"SAMPLE_ID",
			"SAMPLE_ID"
		).empty

	def test_fam_empty(self) -> None:
		assert make_fam(
			pd.DataFrame(columns=["SAMPLE_ID", "SNP"]),
			"SAMPLE_ID",
		).empty

		assert make_fam(
			pd.DataFrame(columns=["SAMPLE_ID", "SNP"]),
			"SAMPLE_ID",
			"SAMPLE_ID",
		).empty

	@pytest.mark.parametrize("data_fam", ["file.pl"], indirect=True)
	def test_fam_raise_columns(self, data_fam: pd.DataFrame) -> None:
		# SID_COL
		with pytest.raises(
				KeyError, match="Data has not in name columns SAMPLE_ID1!"
		):
			make_fam(
				data_fam,
				"SAMPLE_ID1",
				"SAMPLE_ID",
			)

		# FID_COL
		with pytest.raises(
				KeyError, match="Data has not in name columns SAMPLE_ID1!"
		):
			make_fam(
				data_fam,
				"SAMPLE_ID",
				"SAMPLE_ID1"
			)

	@pytest.mark.parametrize("data_fam", ["file2.pl"], indirect=True)
	def test_fam_raises_underscope_sid(self, data_fam: pd.DataFrame) -> None:

		# SID_COL
		with pytest.raises(
				Exception,
				match="Replace in 'Sample ID' columns '_' on another a simbols"
		):
			make_fam(
				data_fam,
				"SAMPLE_ID",
				"SAMPLE_ID"
			)

	@pytest.mark.parametrize("data_fam", ["file3.pl"], indirect=True)
	def test_fam_raises_underscope_fid(self, data_fam: pd.DataFrame) -> None:

		# FID_COL
		with pytest.raises(
				Exception,
				match="Replace in 'Family ID' columns '_' on another a simbols"
		):
			make_fam(
				data_fam,
				"SAMPLE_ID",
				"FAMILY_ID"
			)

	@pytest.mark.parametrize("data_fam", ["file4.pl"], indirect=True)
	def test_fam_check_data(self, data_fam: pd.DataFrame) -> None:
		res = make_fam(
			data_fam,
			"SAMPLE_ID",
			"FAMILY_ID",
			father_col="father",
			mother_col="mother",
			sex_col="sex",
			pheno_col="pheno"
		)

		res2 = make_fam(
			data_fam,
			"SAMPLE_ID",
			"FAMILY_ID",
		)

		assert all(res.father.values == list('1234'))
		assert all(res.mother.values == list('5678'))
		assert all(res.sex.values == list('1210'))
		assert all(res.pheno.values == ['12', '13', '14', '15'])

		assert all(res2.father.values == list('0000'))
		assert all(res2.mother.values == list('0000'))
		assert all(res2.sex.values == list('0000'))
		assert all(res2.pheno.values == ['-9', '-9', '-9', '-9'])
