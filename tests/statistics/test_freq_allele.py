#!/usr/bin/env python
# coding: utf-8
__author__ = "Igor Loschinin (igor.loschinin@gmail.com)"

from . import DIR_DATA
from snptools.src.snplib.statistics import allele_freq

import pytest
import numpy as np
import pandas as pd


@pytest.fixture
def data_df() -> pd.DataFrame:
	return pd.read_pickle(DIR_DATA / "freq/file.pl")
	# [0.   , 0.9  , 0.889]


def data_str() -> list[tuple]:
	return [
		('2212120', 0.714),
		('02011015010000500', 0.2),
		('01110152120222512', 0.6)
	]


class TestAlleleFreq(object):

	def test_allele_freq_df_dtype_obj(self, data_df: pd.DataFrame) -> None:
		data_df.SNP = data_df.SNP.astype(str)
		result = allele_freq(data=data_df, id_col="SNP_NAME", seq_col="SNP")

		assert isinstance(result, pd.DataFrame) and not result.empty
		assert result.SNP.round(6).isin([0.000, 0.900, 0.889]).all()

	def test_allele_freq_df_dtype_int(self, data_df: pd.DataFrame) -> None:
		data_df.SNP = data_df.SNP.astype("int8")
		result = allele_freq(data=data_df, id_col="SNP_NAME", seq_col="SNP")

		assert isinstance(result, pd.DataFrame) and not result.empty
		assert result.SNP.round(6).isin([0.000, 0.900, 0.889]).all()

	def test_allele_freq_df_dtype_float(self, data_df: pd.DataFrame) -> None:
		data_df.SNP = data_df.SNP.astype("float32")
		result = allele_freq(data=data_df, id_col="SNP_NAME", seq_col="SNP")

		assert isinstance(result, pd.DataFrame) and not result.empty
		assert result.SNP.round(6).isin([0.000, 0.900, 0.889]).all()

	def test_allele_freq_df_data_rand_simbols(
			self, data_df: pd.DataFrame
	) -> None:
		data_df.SNP = [
			np.random.choice(["A", "C", "G", "T"])
			for _ in range(data_df.SNP.shape[0])
		]
		assert allele_freq(
			data=data_df, id_col="SNP_NAME", seq_col="SNP"
		) is None

	def test_allele_freq_df_empty(self) -> None:
		with pytest.raises(KeyError):
			allele_freq(
				data=pd.DataFrame(), id_col="SNP_NAME", seq_col="SNP"
			)

	def test_allele_freq_df_empty_only_columns(self) -> None:
		result = allele_freq(
			data=pd.DataFrame(columns=["SNP_NAME", "SNP"]),
			id_col="SNP_NAME",
			seq_col="SNP"
		)

		assert isinstance(result, pd.DataFrame) and result.empty

	def test_allele_freq_df_raises(self, data_df: pd.DataFrame) -> None:
		with pytest.raises(KeyError):
			allele_freq(data=data_df, id_col="SNP_NAME")
			allele_freq(data=data_df, seq_col="SNP")
			allele_freq(data=data_df)

	@pytest.mark.parametrize("data, obs_value", data_str())
	def test_allele_freq_str(self, data: str, obs_value: float) -> None:
		assert allele_freq(data=data) == obs_value

	def test_allele_freq_non_type(self) -> None:
		assert allele_freq(data=1423) is None
