#!/usr/bin/env python
# coding: utf-8
__author__ = "Igor Loschinin (igor.loschinin@gmail.com)"

from . import DIR_FILES
from snptools.src.snplib.format import Snp

import pytest
import pandas as pd


@pytest.fixture
def data_fr(request) -> pd.DataFrame:
	return pd.read_csv(DIR_FILES / f"fsnp/{request.param}", sep="\t")


@pytest.fixture
def obj_snp(request) -> Snp:
	return Snp(fmt=request.param)


class TestSNP(object):

	@pytest.mark.parametrize(
		"obj_snp, data_fr", [("uga", 'file1.txt')], indirect=True
	)
	def test_snp_process_uga_true(
			self, obj_snp: Snp, data_fr: pd.DataFrame
	) -> None:

		obj_snp.process(data_fr)
		assert obj_snp.data is not None and not obj_snp.data.empty
		assert obj_snp.data.SNP.isin([
			'02011015010000500', '01110152120222512'
		]).all()

	@pytest.mark.parametrize("obj_snp", ["uga"], indirect=True)
	def test_snp_process_uga_empty(self, obj_snp: Snp) -> None:

		obj_snp.process(pd.DataFrame(columns=[
			'SNP Name', 'Sample ID', 'Allele1 - AB', 'Allele2 - AB',
			'GC Score', 'GT Score'
		]))
		assert obj_snp.data is not None and obj_snp.data.empty

	@pytest.mark.parametrize(
		"obj_snp, data_fr",
		[("uga", 'file1.txt'), (None, 'file1.txt')],
		indirect=True
	)
	def test_snp_process_raises(
			self, obj_snp: Snp, data_fr: pd.DataFrame
	) -> None:

		with pytest.raises(KeyError):
			obj_snp.process(pd.DataFrame(columns=[
				'SNP Name1', 'Sample ID1', 'Allele1 - AB', 'Allele2 - AB',
				'GC Score', 'GT Score'
			]))

		assert obj_snp.data is None

	@pytest.mark.parametrize(
		"obj_snp, data_fr", [(None, 'file1.txt')], indirect=True
	)
	def test_snp_process_df(
			self, obj_snp: Snp, data_fr: pd.DataFrame
	) -> None:

		obj_snp.process(data_fr)
		assert obj_snp.data is not None and not obj_snp.data.empty

	@pytest.mark.parametrize("obj_snp", [None], indirect=True)
	def test_snp_process_df_empty(self, obj_snp: Snp) -> None:

		obj_snp.process(pd.DataFrame(columns=[
			'SNP Name', 'Sample ID', 'Allele1 - AB', 'Allele2 - AB',
			'GC Score', 'GT Score'
		]))
		assert obj_snp.data is not None and obj_snp.data.empty

	@pytest.mark.parametrize(
		"obj_snp, data_fr", [("uga", 'file1.txt')], indirect=True
	)
	def test_snp_to_file_uga1(
			self, obj_snp: Snp, data_fr: pd.DataFrame, tmp_path
	) -> None:
		"""
		The name sample_id is one length
		"""

		_dir_sub = tmp_path / "sub"
		_dir_sub.mkdir()
		_file_save = _dir_sub / "data_snp.csv"

		obj_snp.process(data_fr)
		assert obj_snp.data is not None and not obj_snp.data.empty

		obj_snp.to_file(_file_save)
		assert _file_save.is_file() and _file_save.exists()
		assert (
			_file_save.read_text() ==
			"14814 02011015010000500\n14815 01110152120222512\n"
		)

	@pytest.mark.parametrize(
		"obj_snp, data_fr", [("uga", 'file2.txt')], indirect=True
	)
	def test_snp_to_file_uga2(
			self, obj_snp: Snp, data_fr: pd.DataFrame, tmp_path
	) -> None:
		"""
		The name sample_id of different length
		"""

		_dir_sub = tmp_path / "sub"
		_dir_sub.mkdir()
		_file_save = _dir_sub / "data_snp.csv"

		obj_snp.process(data_fr)
		assert obj_snp.data is not None and not obj_snp.data.empty

		obj_snp.to_file(_file_save)
		assert _file_save.is_file() and _file_save.exists()
		assert (
			_file_save.read_text() ==
			"14814qwert 02011015010000500\n14815      01110152120222512\n"
		)
