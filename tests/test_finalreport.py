#!/usr/bin/env python
# coding: utf-8

from . import PATH_DIR_FILES

from pathlib import Path
from finalreport import FinalReport


class TestFinalReport(object):

	def test_handle_1(self) -> None:
		""" if both files do not exist """
		obj = FinalReport()

		print(Path().cwd())

		assert not obj.handle(
			Path(f"{PATH_DIR_FILES}/fr/f.txt"),
			Path(f"{PATH_DIR_FILES}/fr/f.xlsx")
		)

	def test_handle_2(self) -> None:
		""" if the file to convert does not exist """
		obj = FinalReport()

		assert not obj.handle(
			Path(f"{PATH_DIR_FILES}/fr/file1.txt"),
			Path(f"{PATH_DIR_FILES}/fr/f.xlsx")
		)

	def test_handle_3(self) -> None:
		""" if the data does not contain header data """
		obj = FinalReport()

		obj.handle(
			Path(f"{PATH_DIR_FILES}/fr/file2.txt"),
			Path(f"{PATH_DIR_FILES}/fr/file2.xlsx")
		)

		assert len(obj.header) == 0 and \
		       not obj.snp_data.empty

	def test_handle_4(self) -> None:
		""" if the file contains only header and field
		names """
		obj = FinalReport()

		obj.handle(
			Path(f"{PATH_DIR_FILES}/fr/file3.txt"),
			Path(f"{PATH_DIR_FILES}/fr/file3.xlsx")
		)

		assert obj.snp_data is not None and obj.snp_data.empty

	def test_handle_5(self) -> None:
		""" if the data file is empty """
		obj = FinalReport()

		assert not obj.handle(
			Path(f"{PATH_DIR_FILES}/fr/file5.txt"),
			Path(f"{PATH_DIR_FILES}/fr/file5.xlsx")
		)

	def test_handle_6(self) -> None:
		""" if the conversion file is empty """
		obj = FinalReport()

		assert obj.handle(
			Path(f"{PATH_DIR_FILES}/fr/file6.txt"),
			Path(f"{PATH_DIR_FILES}/fr/file6.xlsx")
		)

	def test_handle_7(self) -> None:
		""" if the data file is not needed convert ID name """
		obj = FinalReport()

		obj.handle(
			Path(f"{PATH_DIR_FILES}/fr/file4.txt"),
			None
		)

		assert not obj.snp_data.empty
		assert len(obj.header) != 0

	def test_handle_8(self) -> None:
		""" if files exist """
		obj = FinalReport()

		assert obj.handle(
			Path(f"{PATH_DIR_FILES}/fr/file1.txt"),
			Path(f"{PATH_DIR_FILES}/fr/file1.xlsx")
		)

	def test_sample_allele_ab(self) -> None:
		obj = FinalReport(allele="AB")

		assert obj.handle(
			Path(f"{PATH_DIR_FILES}/fr/file4.txt"),
			None
		)

	def test_sample_allele_forward(self) -> None:
		obj = FinalReport(allele="Forward")

		assert obj.handle(
			Path(f"{PATH_DIR_FILES}/fr/file4.txt"),
			None
		)

	def test_sample_allele_top(self) -> None:
		obj = FinalReport(allele="Top")

		assert obj.handle(
			Path(f"{PATH_DIR_FILES}/fr/file4.txt"),
			None
		)

	def test_sample_allele_list1(self) -> None:
		obj = FinalReport(allele=["AB", "Top"])

		assert obj.handle(
			Path(f"{PATH_DIR_FILES}/fr/file4.txt"),
			None
		)

	def test_sample_allele_list2(self) -> None:
		obj = FinalReport(allele=["AB"])

		assert obj.handle(
			Path(f"{PATH_DIR_FILES}/fr/file4.txt"),
			None
		)

	def test_sample_allele_not_exist(self) -> None:
		obj = FinalReport(allele="GG")

		assert obj.handle(
			Path(f"{PATH_DIR_FILES}/fr/file4.txt"),
			None
		)
