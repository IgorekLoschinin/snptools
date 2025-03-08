#!/usr/bin/env python
# coding: utf-8
__author__ = "Igor Loschinin (igor.loschinin@gmail.com)"


from . import DIR_FILES
from snptools.src.snplib.finalreport import FinalReport

import pytest


@pytest.fixture
def report(request) -> FinalReport:
	return FinalReport(**request.param)


class TestFinalReport(object):

	@pytest.mark.parametrize("report", [{'allele': None}], indirect=True)
	def test_handle_1(self, report: FinalReport) -> None:
		""" If both files do not exist """

		assert not report.handle(
			DIR_FILES / "fr/f.txt", DIR_FILES / "fr/f.xlsx",
		)

	@pytest.mark.parametrize("report", [{'allele': None}], indirect=True)
	def test_handle_8(self, report: FinalReport) -> None:
		""" If files exist """

		assert report.handle(
			DIR_FILES / "fr/file1.txt", DIR_FILES / "fr/file1.xlsx",
		)

		assert len(report.header) != 0 and not report.snp_data.empty

	@pytest.mark.parametrize("report", [{'allele': None}], indirect=True)
	def test_handle_2(self, report: FinalReport) -> None:
		""" If the file to convert does not exist """

		assert not report.handle(
			DIR_FILES / "fr/file1.txt", DIR_FILES / "fr/f.xlsx",
		)

	@pytest.mark.parametrize("report", [{'allele': None}], indirect=True)
	def test_handle_3(self, report: FinalReport) -> None:
		""" If the data does not contain header data """

		report.handle(
			DIR_FILES / "fr/file2.txt", DIR_FILES / "fr/file2.xlsx",
		)

		assert len(report.header) == 0 and not report.snp_data.empty

	@pytest.mark.parametrize("report", [{'allele': None}], indirect=True)
	def test_handle_4(self, report: FinalReport) -> None:
		""" If the file contains only header and field names """

		report.handle(
			DIR_FILES / "fr/file3.txt", DIR_FILES / "fr/file3.xlsx",
		)

		assert report.snp_data is not None and report.snp_data.empty

	@pytest.mark.parametrize("report", [{'allele': None}], indirect=True)
	def test_handle_5(self, report: FinalReport) -> None:
		""" If the data file is empty """

		report.handle(
			DIR_FILES / "fr/file5.txt", DIR_FILES / "fr/file5.xlsx",
		)

		assert len(report.header) == 0
		assert report.snp_data is not None and report.snp_data.empty

	@pytest.mark.parametrize("report", [{'allele': None}], indirect=True)
	def test_handle_6(self, report: FinalReport) -> None:
		""" If the conversion file is empty """

		assert report.handle(
			DIR_FILES / "fr/file6.txt", DIR_FILES / "fr/file6.xlsx",
		)

		assert not report.snp_data.empty
		assert len(report.header) != 0

	@pytest.mark.parametrize("report", [{'allele': None}], indirect=True)
	def test_handle_7(self, report: FinalReport) -> None:
		""" If the data file is not needed to convert ID name """

		report.handle(DIR_FILES / "fr/file4.txt", None)

		assert not report.snp_data.empty
		assert len(report.header) != 0

	@pytest.mark.parametrize("report", [{'allele': None}], indirect=True)
	def test_allele_none(self, report: FinalReport) -> None:
		report.handle(DIR_FILES / "fr/file4.txt", None)

		_fields = [
			'SNP Name', 'Sample ID', 'Allele1 - Forward', 'Allele2 - Forward',
			'Allele1 - Top', 'Allele2 - Top', 'Allele1 - AB', 'Allele2 - AB',
			'GC Score', 'X', 'Y'
		]

		assert report.snp_data.columns.difference(_fields).empty

	@pytest.mark.parametrize("report", [{'allele': "AB"}], indirect=True)
	def test_sample_allele_ab(self, report: FinalReport) -> None:
		report.handle(DIR_FILES / "fr/file4.txt", None)

		_fields = [
			'SNP Name', 'Sample ID', 'Allele1 - AB', 'Allele2 - AB',
			'GC Score', 'X', 'Y'
		]

		assert report.snp_data.columns.difference(_fields).empty

	@pytest.mark.parametrize("report", [{'allele': "Forward"}], indirect=True)
	def test_sample_allele_forward(self, report: FinalReport) -> None:
		report.handle(DIR_FILES / "fr/file4.txt", None)

		_fields = [
			'SNP Name', 'Sample ID', 'Allele1 - Forward', 'Allele2 - Forward',
			'GC Score', 'X', 'Y'
		]

		assert report.snp_data.columns.difference(_fields).empty

	@pytest.mark.parametrize("report", [{'allele': "Top"}], indirect=True)
	def test_sample_allele_top(self, report: FinalReport) -> None:
		report.handle(DIR_FILES / "fr/file4.txt", None)

		_fields = [
			'SNP Name', 'Sample ID', 'Allele1 - Top', 'Allele2 - Top',
			'GC Score', 'X', 'Y'
		]

		assert report.snp_data.columns.difference(_fields).empty

	@pytest.mark.parametrize(
		"report", [{'allele': ["AB", "Top"]}], indirect=True
	)
	def test_sample_allele_list1(self, report: FinalReport) -> None:
		report.handle(DIR_FILES / "fr/file4.txt", None)

		_fields = [
			'SNP Name', 'Sample ID', 'Allele1 - Top', 'Allele2 - Top',
			'Allele1 - AB', 'Allele2 - AB', 'GC Score', 'X', 'Y'
		]

		assert report.snp_data.columns.difference(_fields).empty

	@pytest.mark.parametrize("report", [{'allele': ["AB"]}], indirect=True)
	def test_sample_allele_list2(self, report: FinalReport) -> None:
		report.handle(DIR_FILES / "fr/file4.txt", None)

		_fields = [
			'SNP Name', 'Sample ID', 'Allele1 - AB', 'Allele2 - AB',
			'GC Score', 'X', 'Y'
		]

		assert report.snp_data.columns.difference(_fields).empty

	@pytest.mark.parametrize(
		"report", [{'allele': ("AB", "Top")}], indirect=True
	)
	def test_sample_allele_tuple(self, report: FinalReport) -> None:
		report.handle(DIR_FILES / "fr/file4.txt", None)

		_fields = [
			'SNP Name', 'Sample ID', 'Allele1 - Top', 'Allele2 - Top',
			'Allele1 - AB', 'Allele2 - AB', 'GC Score', 'X', 'Y'
		]

		assert report.snp_data.columns.difference(_fields).empty

	@pytest.mark.parametrize(
		"report", [{'allele': {"AB", "Top"}}], indirect=True
	)
	def test_sample_allele_set(self, report: FinalReport) -> None:
		report.handle(DIR_FILES / "fr/file4.txt", None)

		_fields = [
			'SNP Name', 'Sample ID', 'Allele1 - Top', 'Allele2 - Top',
			'Allele1 - AB', 'Allele2 - AB', 'GC Score', 'X', 'Y'
		]

		assert report.snp_data.columns.difference(_fields).empty

	@pytest.mark.parametrize("report", [{'allele': "GG"}], indirect=True)
	def test_sample_allele_not_exist(self, report: FinalReport) -> None:

		with pytest.raises(
				Exception, match="Error. Allele GG not in data."
		):
			report.handle(DIR_FILES / "fr/file4.txt", None)

	@pytest.mark.parametrize("report", [{'allele': ["AB"]}], indirect=True)
	def test_on_ru_symbol(self, report: FinalReport) -> None:
		""" An error is checked if the name of the number is Kirilitsa """

		with pytest.raises(
			Exception, match="Error. Unique keys contain Cyrillic alphabet."
		):
			report.handle(
				DIR_FILES / "fr/file7.txt", DIR_FILES / "fr/file7.xlsx"
			)

		assert report.snp_data is None and len(report.header) == 0

	@pytest.mark.parametrize("report", [{'allele': ["AB"]}], indirect=True)
	def test_on_nan_in_sid_file_xlsx(self, report: FinalReport) -> None:
		""" Checking for processing empty values in SID """

		report.handle(
			DIR_FILES / "fr/file8.txt",
			DIR_FILES / "fr/file8.xlsx"
		)

		assert report.snp_data is not None and not report.snp_data.empty
		assert report.snp_data['Sample ID'].isna().any()

	@pytest.mark.parametrize("report", [{'allele': "AB"}], indirect=True)
	def test_on_nan_in_fr(self, report: FinalReport) -> None:
		""" Checking for missing values in SID """
		report.handle(
			DIR_FILES / "fr/file9.txt",
			DIR_FILES / "fr/file9.xlsx"
		)

		assert report.snp_data is not None and not report.snp_data.empty
		assert report.snp_data['Sample ID'].isna().any()

	@pytest.mark.parametrize(
		"report",
		[{'allele': "AB", 'usecols': ['SNP Name', 'Sample ID']}],
		indirect=True
	)
	def test_raise_use_allele_and_usecols(self, report: FinalReport) -> None:
		""" Checking for allele and usecols sharing error """

		with pytest.raises(
			Exception, match="Error. Usecols is used for allele is none."
		):
			report.handle(
				DIR_FILES / "fr/file1.txt", DIR_FILES / "fr/file1.xlsx"
			)

		assert report.snp_data is None and len(report.header) == 0

	@pytest.mark.parametrize(
		"report",
		[{'usecols': ['SNP Name', 'Sample ID']}],
		indirect=True
	)
	def test_using_usecols(self, report: FinalReport) -> None:
		""" Checking for use setting usecols """

		_fields = ['SNP Name', 'Sample ID']

		report.handle(
			DIR_FILES / "fr/file1.txt", DIR_FILES / "fr/file1.xlsx"
		)

		assert not report.snp_data.empty and len(report.header) != 0
		assert report.snp_data.columns.difference(_fields).empty

	@pytest.mark.parametrize(
		"report",
		[{
			'usecols': ['SNP Name', 'Sample ID'],
			'dtype': {'SNP Name': 'category'}
		}],
		indirect=True
	)
	def test_using_dtype_and_combo_usecols(self, report: FinalReport) -> None:
		""" Checking the use of dtype """

		_fields = ['SNP Name', 'Sample ID']

		report.handle(
			DIR_FILES / "fr/file1.txt", DIR_FILES / "fr/file1.xlsx"
		)

		assert not report.snp_data.empty and len(report.header) != 0
		assert report.snp_data.columns.difference(_fields).empty
		assert report.snp_data['SNP Name'].dtypes == 'category'

	@pytest.mark.parametrize(
		"report",
		[{
			'allele': "AB",
			'dtype': {'SNP Name': 'category'}
		}],
		indirect=True
	)
	def test_using_dtype_and_combo_usecols(self, report: FinalReport) -> None:
		""" Checking the use combo dtype and allele """

		report.handle(
			DIR_FILES / "fr/file1.txt", DIR_FILES / "fr/file1.xlsx"
		)

		assert not report.snp_data.empty and len(report.header) != 0
		assert report.snp_data['SNP Name'].dtypes == 'category'
