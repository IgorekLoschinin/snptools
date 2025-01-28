#!/usr/bin/env python
# coding: utf-8
__author__ = "Igor Loschinin (igor.loschinin@gmail.com)"

from . import DIR_DATA
from snptools.src.snplib.statistics import call_rate

import pytest
import numpy as np
import pandas as pd


@pytest.fixture
def data_df(request) -> pd.DataFrame:
	match request.param:
		case "cra":
			return pd.read_pickle(DIR_DATA / "cr/file_cra.pl")

		case "crm":
			return pd.read_pickle(DIR_DATA / "cr/file_crm.pl")


@pytest.fixture
def data_str() -> list[str]:
	return ['02011015010000500', '01110152120222512']


class TestCallRateAnimal(object):

	@pytest.mark.parametrize("data_df", ["cra"], indirect=True)
	def test_cra_datafame_dtype_obj(self, data_df: pd.DataFrame) -> None:
		data_df.SNP = data_df.SNP.astype(str)
		result = call_rate(data=data_df, id_col="SAMPLE_ID", snp_col="SNP")

		assert isinstance(result, pd.DataFrame) and not result.empty
		assert result.SNP.round(6).isin([0.882353, 0.882353]).all()

	@pytest.mark.parametrize("data_df", ["cra"], indirect=True)
	def test_cra_datafame_dtype_int(self, data_df: pd.DataFrame) -> None:
		data_df.SNP = data_df.SNP.astype("int8")
		result = call_rate(data=data_df, id_col="SAMPLE_ID", snp_col="SNP")

		assert isinstance(result, pd.DataFrame) and not result.empty
		assert result.SNP.round(6).isin([0.882353, 0.882353]).all()

	@pytest.mark.parametrize("data_df", ["cra"], indirect=True)
	def test_cra_datafame_dtype_float(self, data_df: pd.DataFrame) -> None:
		data_df.SNP = data_df.SNP.astype("float32")
		result = call_rate(data=data_df, id_col="SAMPLE_ID", snp_col="SNP")

		assert isinstance(result, pd.DataFrame) and not result.empty
		assert result.SNP.round(6).isin([0.882353, 0.882353]).all()

	@pytest.mark.parametrize("data_df", ["cra"], indirect=True)
	def test_cra_datafame_dtype_random_simbols(
		self, data_df: pd.DataFrame
	) -> None:
		data_df.SNP = [
			np.random.choice(["A", "C", "G", "T"])
			for _ in range(data_df.SNP.shape[0])
		]
		result = call_rate(data=data_df, id_col="SAMPLE_ID", snp_col="SNP")

		assert result is None

	def test_cra_datafame_empty1(self) -> None:
		with pytest.raises(KeyError):
			call_rate(data=pd.DataFrame(), id_col="SAMPLE_ID", snp_col="SNP")

	def test_cra_datafame_empty2(self) -> None:
		result = call_rate(
			data=pd.DataFrame(columns=["SAMPLE_ID", "SNP"]),
			id_col="SAMPLE_ID",
			snp_col="SNP"
		)

		assert isinstance(result, pd.DataFrame) and result.empty

	@pytest.mark.parametrize("data_df", ["cra"], indirect=True)
	def test_cra_datafame_fail(self, data_df: pd.DataFrame) -> None:
		with pytest.raises(KeyError):
			call_rate(data=data_df, id_col="SAMPLE_ID")
			call_rate(data=data_df, snp_col="SNP")
			call_rate(data=data_df)

	def test_cra_str_int(self, data_str: list[str]) -> None:
		for sequence in data_str:
			assert call_rate(data=sequence) == 0.882353

	def test_cra_str_simbols(self) -> None:
		data_str = ['GCATGAGGTATACTCTA', 'CGCCATGCTGTATATCC']

		for sequence in data_str:
			assert call_rate(data=sequence) is None

	def test_cra_str_empty(self) -> None:
		assert call_rate(data="") is None

	def test_cra_str_mixid(self) -> None:
		assert call_rate(data="GCATGAG3G4T6A67TACTCTA") is None


class TestCallRateMarker(object):

	@pytest.mark.parametrize("data_df", ["crm"], indirect=True)
	def test_crm_datafame_dtype_obj(self, data_df: pd.DataFrame) -> None:
		data_df.SNP = data_df.SNP.astype(str)
		result = call_rate(data=data_df, id_col="SNP_NAME", snp_col="SNP")

		assert isinstance(result, pd.DataFrame) and not result.empty
		assert result.SNP.round(6).isin([0.727273, 0.909091, 0.818182]).all()

	@pytest.mark.parametrize("data_df", ["crm"], indirect=True)
	def test_crm_datafame_dtype_int(self, data_df: pd.DataFrame) -> None:
		data_df.SNP = data_df.SNP.astype("int8")
		result = call_rate(data=data_df, id_col="SNP_NAME", snp_col="SNP")

		assert isinstance(result, pd.DataFrame) and not result.empty
		assert result.SNP.round(6).isin([0.727273, 0.909091, 0.818182]).all()

	@pytest.mark.parametrize("data_df", ["crm"], indirect=True)
	def test_crm_datafame_dtype_float(self, data_df: pd.DataFrame) -> None:
		data_df.SNP = data_df.SNP.astype("float32")
		result = call_rate(data=data_df, id_col="SNP_NAME", snp_col="SNP")

		assert isinstance(result, pd.DataFrame) and not result.empty
		assert result.SNP.round(6).isin([0.727273, 0.909091, 0.818182]).all()

	@pytest.mark.parametrize("data_df", ["crm"], indirect=True)
	def test_crm_datafame_dtype_random_simbols(
		self, data_df: pd.DataFrame
	) -> None:
		data_df.SNP = [
			np.random.choice(["A", "C", "G", "T"])
			for _ in range(data_df.SNP.shape[0])
		]
		result = call_rate(data=data_df, id_col="SNP_NAME", snp_col="SNP")

		assert result is None

	def test_crm_datafame_empty1(self) -> None:
		with pytest.raises(KeyError):
			call_rate(data=pd.DataFrame(), id_col="SNP_NAME", snp_col="SNP")

	def test_crm_datafame_empty2(self) -> None:
		result = call_rate(
			data=pd.DataFrame(columns=["SNP_NAME", "SNP"]),
			id_col="SNP_NAME",
			snp_col="SNP"
		)

		assert isinstance(result, pd.DataFrame) and result.empty

	@pytest.mark.parametrize("data_df", ["crm"], indirect=True)
	def test_crm_datafame_fail(self, data_df: pd.DataFrame) -> None:
		with pytest.raises(KeyError):
			call_rate(data=data_df, id_col="SNP_NAME")
			call_rate(data=data_df, snp_col="SNP")
			call_rate(data=data_df)

	def test_crm_str_simbols(self) -> None:
		data_str = ['GCATGAGGTATACTCTA', 'CGCCATGCTGTATATCC']

		for sequence in data_str:
			assert call_rate(data=sequence) is None

	def test_crm_str_empty(self) -> None:
		assert call_rate(data="") is None

	def test_crm_str_mixid(self) -> None:
		assert call_rate(data="GCATGAG3G4T6A67TACTCTA") is None
