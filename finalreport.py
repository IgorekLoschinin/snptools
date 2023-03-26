#!/usr/bin/env python
# coding: utf-8
import pandas as pd

from pathlib import Path


class FinalReport(object):
	""" File that contains SNP information

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

	__PATTERN_HEADER = "[Header]"
	__PATTERN_DATA = "[Data]"
	__LIST_ALLELES = ("AB", "Forward", "Top", "Plus", "Design")

	def __init__(
			self,
			allele: str | list | None = "AB",
			sep: str = "\t"
	) -> None:
		"""
		:param allele: - One of the possible bases/sequences that can occur at
			the variant - a locus in the genome where there are differences
			between individuals. Type: 'AB', 'Forward', 'Top'
		:param sep: - Delimiter to use. Default value: "\t"
		"""

		self._delimiter = sep
		self._full_data = None

		self.__header = dict()
		self.__snp_data = None
		self.__allele = allele
		self.__map_rn = None

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
		and data

		:param file_rep: - The file FinalReport.txt or another name.
		:param conv_file: - The file that contains IDs of registration numbers
			of animals

		:return: - Returns true if file processing was successful, false if
			there were errors
		"""

		try:

			if isinstance(file_rep, str):
				file_rep = Path(file_rep)

			if not file_rep.is_file() and not file_rep.exists():
				return False

			if conv_file is not None:
				if isinstance(conv_file, str):
					conv_file = Path(conv_file)

				if not conv_file.is_file() and not conv_file.exists():
					return False

				self.__convert_s_id(conv_file)

			if not self.read(file_rep):
				return False

			self.__handler_header()
			self.__handler_data()

			if self.__map_rn is not None:
				self.__snp_data['Sample ID'] = \
					self.__snp_data['Sample ID'].map(
						dict(zip(self.__map_rn.SID, self.__map_rn.UNIQ_KEY))
					)

			return True

		except Exception as e:
			return False

	def read(self, file_rep: Path) -> bool:
		""" Reading data from the final_report file

		:param file_rep: - path, pointer to the file to be read
		:return: - Returns true if the read was successful, false if it failed
		"""
		try:
			self._full_data = file_rep.read_text().strip().split("\n")
			return True

		except Exception as e:
			return False

	def __handler_header(self) -> None:
		""" Processes data from a file, selects meta-information """

		for line in self._full_data:
			if line == self.__class__.__PATTERN_DATA:
				return

			if line != self.__class__.__PATTERN_HEADER:
				key = line.strip().split("\t")[0]
				value = line.strip().split("\t")[1]

				self.__header[key] = value

	def __handler_data(self) -> None:
		""" Processes data and forms an array for further processing """

		temp = 1
		for line in self._full_data:
			if line == self.__class__.__PATTERN_DATA:
				break
			temp += 1

		self.__snp_data = pd.DataFrame(
			[
				item_data.split(f"{self._delimiter}")
				for item_data in self._full_data[temp + 1:]
			],
			columns=self._full_data[temp].split(f"{self._delimiter}")
		)

	def __convert_s_id(self, path_file: Path) -> None:
		"""Converts sample id which is in FinalReport to animal registration
		number

		:param path_file: - xlsx file with animal numbers label
		"""

		self.__map_rn = pd.read_excel(
			path_file,
			header=None,
			names=['SID', 'UNIQ_KEY'],
			dtype={'SID': str}
		)

		if self.__map_rn.empty:
			self.__map_rn = None

		else:
			if self.__map_rn.UNIQ_KEY.isna().any():
				self.__map_rn.fillna('unknown', inplace=True)
