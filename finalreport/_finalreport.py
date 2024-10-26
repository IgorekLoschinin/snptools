#!/usr/bin/env python
# coding: utf-8
__author__ = "Igor Loschinin (igor.loschinin@gmail.com)"
__all__ = ("FinalReport",)

from pathlib import Path
from functools import reduce

import re
import pandas as pd


class FinalReport(object):
	""" File that contains SNP information.

	:argument allele: A variant form of a single nucleotide polymorphism
		(SNP), a specific polymorphic site or a whole gene detectable at
		a locus.  Type: 'AB', 'Forward', 'Top', 'Plus', 'Design'
	:argument sep: Delimiter to use. Default value: "\\t"

	Example:
		[Header]
		GSGT Version	2.0.4
		Processing Date	10/14/2021 4:02 PM
		Content		BovineSNP50_v3_A1.bpm
		Num SNPs	53218
		Total SNPs	53218
		Num Samples	3
		Total Samples	3
		[Data]
		SNP Name  Sample ID  Allele1 - AB  Allele2 - AB  GC Score  GT Score
		ABCA12	1	A	A	0.4048	0.8164
		APAF1	1	B	B	0.9067	0.9155
		...
	"""

	__PATTERN_HEADER = re.compile(r'(^\[Header\])')
	__PATTERN_DATA = re.compile(r'(^\[Data\])')

	def __init__(
			self,
			allele: str | list | None = None,
			sep: str = "\t"
	) -> None:
		self._delimiter = sep
		self._full_data = None

		self.__header = {}
		self.__snp_data = None
		self.__allele = allele
		self._map_rn = None

	@property
	def header(self) -> dict:
		return self.__header

	@property
	def snp_data(self) -> pd.DataFrame | None:
		return self.__snp_data

	def handle(
			self, file_rep: Path | str, conv_file: Path | str = None
	) -> bool:
		""" Processes the FinalReport.txt file. Highlights meta information
		and data.

		:param file_rep: The file FinalReport.txt or another name.
		:param conv_file: The file that contains IDs of registration numbers
			of animals.
		:return: Returns true if file processing was successful, false if
			there were errors.
		"""

		try:

			if isinstance(file_rep, str):
				file_rep = Path(file_rep)

			if not file_rep.is_file() and not file_rep.exists():
				return False

			# Processing conversion file
			if conv_file is not None:
				if isinstance(conv_file, str):
					conv_file = Path(conv_file)

				if not conv_file.is_file() and not conv_file.exists():
					return False

				self.__convert_s_id(conv_file)

			# Processing report file
			if not self.read(file_rep):
				return False

			if self._full_data is None:
				raise Exception("Not data in file FinalReport.txt")

			self.__handler_header()
			self.__handler_data()

			if self._map_rn is not None:
				self.__snp_data['Sample ID'] = \
					self.__snp_data['Sample ID'].map(
						dict(zip(self._map_rn.SID, self._map_rn.UNIQ_KEY))
					)

		except Exception as e:
			raise e

		return True

	def read(self, file_rep: Path) -> bool:
		""" Reading data from the final_report file

		:param file_rep: path, pointer to the file to be read.
		:return: Returns true if the read was successful, false if it failed.
		"""
		try:
			if len(data := file_rep.read_text()) != 0:
				self._full_data = data.strip().split("\n")
				return True

			self._full_data = None

		except Exception as e:
			return False

		return True

	def __handler_header(self) -> None:
		""" Processes data from a file, selects meta-information. """

		for line in self._full_data:
			if self.__class__.__PATTERN_DATA.findall(line):
				return

			if self.__class__.__PATTERN_HEADER.findall(line):
				continue

			key = line.strip().split("\t")[0]
			value = line.strip().split("\t")[1]

			self.__header[key] = value

	def __handler_data(self) -> None:
		""" Processes data and forms an array for further processing. """

		temp = 1
		for line in self._full_data:
			if self.__class__.__PATTERN_DATA.findall(line):
				break
			temp += 1

		names_col = self.__sample_by_allele(
			self._full_data[temp].split(f"{self._delimiter}")
		)

		if names_col is None:
			raise Exception(f"Error. Allele {self.__allele} not in data.")

		self.__snp_data = pd.DataFrame(
			[
				item_data.split(f"{self._delimiter}")
				for item_data in self._full_data[temp + 1:]
			],
			columns=self._full_data[temp].split(f"{self._delimiter}")
		)[names_col]

	def __sample_by_allele(self, names: list[str]) -> list[str] | None:
		""" Method that generates a list of field names choosing which alleles
		to keep

		:param names: List of field names in the report file.
		:return: Returns a filtered list of fields by alleles.
		"""

		allele_templ = r'(^Allele\d\s[:-]\s{}\b)'

		match self.__allele:
			case None:
				return names

			case str():
				allele_pattern = re.compile(
					allele_templ.format(self.__allele)
				)

			case list() | tuple() | set():
				allele_pattern = re.compile(
					allele_templ.format("|".join(self.__allele))
				)
			case _:
				return None

		lst_allele = reduce(
			lambda i, j: i + j,
			[allele_pattern.findall(item) for item in names]
		)

		if len(lst_allele) == 0:
			return None

		exclude_alleles = [
			item for item in names
			if item.startswith("Allele") and item not in lst_allele
		]

		return list(filter(
			lambda x: True if x not in exclude_alleles else False, names
		))

	def __convert_s_id(self, path_file: Path) -> None:
		"""Converts sample id which is in FinalReport to animal registration
		number.

		:param path_file: xlsx file with animal numbers label
		"""

		self._map_rn = pd.read_excel(
			path_file,
			header=None,
			names=['SID', 'UNIQ_KEY', 'SEX'],
			dtype={'SID': str},
			index_col=False
		)

		if self._map_rn.empty:
			self._map_rn = None
			return

		self._map_rn.SID = self._map_rn.SID.str.strip()
		self._map_rn.UNIQ_KEY = self._map_rn.UNIQ_KEY.str.strip()

		if self._check_on_ru_symbols(self._map_rn.UNIQ_KEY):
			raise Exception("Error. Unique keys contain Cyrillic alphabet.")

		if self._map_rn.UNIQ_KEY.isna().any():
			self._map_rn.fillna('unknown', inplace=True)

	@staticmethod
	def _check_on_ru_symbols(seq: pd.Series) -> bool | None:
		"""

		:param seq:
		:return:
		"""

		return any(seq.apply(lambda x: bool(re.search('[а-яА-Я]', x))))


