#!/usr/bin/env python
# coding: utf-8
__author__ = "Igor Loschinin (igor.loschinin@gmail.com)"

from . import DIR_FILES
from snptools.src.snplib.format import make_lgen

import pytest
import pandas as pd


@pytest.fixture
def data_lgen(request) -> pd.DataFrame:
	return pd.read_pickle(DIR_FILES / f"fplink/lgen/{request.param}")


class TestPlinkFormatLgen(object):

	@pytest.mark.parametrize("data_lgen", ["file.pl"], indirect=True)
	def test_lgen_true(self, data_lgen: pd.DataFrame) -> None:
		assert not make_lgen(
			data_lgen,
			"Sample ID",
			"SNP Name",
			["Allele1 - AB", "Allele2 - AB"]
		).empty

	def test_lgen_empty(self) -> None:
		assert make_lgen(
			pd.DataFrame(columns=[
				"Sample ID", "SNP Name", "Allele1 - AB", "Allele2 - AB"
			]),
			"Sample ID",
			"SNP Name",
			["Allele1 - AB", "Allele2 - AB"]
		).empty

	@pytest.mark.parametrize("data_lgen", ["file.pl"], indirect=True)
	def test_lgen_raise_columns(self, data_lgen: pd.DataFrame) -> None:

		with pytest.raises(
			Exception,
			match="Replace in 'Sample ID' columns '_' on another a simbols"
		):
			res1 = data_lgen.copy(deep=True)
			res1["Sample ID"] = res1["Sample ID"] + "_"

			make_lgen(
				res1,
				"Sample ID",
				"SNP Name",
				["Allele1 - AB", "Allele2 - AB"]
			)

		with pytest.raises(
			Exception,
			match="Replace in 'Family ID' columns '_' on another a simbols"
		):
			res1 = data_lgen.copy(deep=True)
			res1["Family ID"] = res1["Sample ID"] + "_"

			make_lgen(
				res1,
				"Sample ID",
				"SNP Name",
				["Allele1 - AB", "Allele2 - AB"],
				fid_col="Family ID"
			)

		# SID
		with pytest.raises(KeyError):
			make_lgen(
				data_lgen,
				"Sample ID1",
				"SNP Name",
				["Allele1 - AB", "Allele2 - AB"],
				fid_col="Family ID"
			)

		# FID_COL
		with pytest.raises(KeyError):
			make_lgen(
				data_lgen,
				"Sample ID",
				"SNP Name",
				["Allele1 - AB", "Allele2 - AB"],
				fid_col="Family ID"
			)

		# SNP name
		with pytest.raises(KeyError):
			make_lgen(
				data_lgen,
				"Sample ID",
				"SNP Name1",
				["Allele1 - AB", "Allele2 - AB"]
			)

		# Alleles
		with pytest.raises(KeyError):
			make_lgen(
				data_lgen,
				"Sample ID",
				"SNP Name",
				["Allele1 - AB1", "Allele2 - AB1"]
			)
