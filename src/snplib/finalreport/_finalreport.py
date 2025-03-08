#!/usr/bin/env python
# coding: utf-8
__author__ = "Igor Loschinin (igor.loschinin@gmail.com)"
__all__ = ("FinalReport",)

import re
from functools import reduce
from pathlib import Path

import pandas as pd
from numpy import nan


class FinalReport(object):
	""" File that contains SNP information. File processing is triggered by the
	handle method. If values in 'SID' or 'UNIQ_KEY' were missing in the xlsx
	conversion file, the processed data will contain NAN values.

	:param allele: A variant form of a single nucleotide polymorphism (SNP), a
		specific polymorphic site or a whole gene detectable at a locus. Type:
		'AB', 'Forward', 'Top', 'Plus', 'Design'.
	:param sep: Delimiter to use. Default value: "\\t".
	:param usecols: Selection of fields for reading. Accelerates processing
		and reduces memory.
	:param dtype: Data type(s) to apply to either the whole dataset or
		individual columns. E.g., {'a': np.float64, 'b': np.int32, 'c': 'Int64'}.

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

	__PATTERN_HEADER = re.compile(r'(^\[Header])')
	__PATTERN_DATA = re.compile(r'(^\[Data])')

	__slots__ = (
		"_delimiter",
		"__allele",
		"__usecols",
		"__dtype",
		"__snp_data",
		"__header",
		"_map_rn",
	)

	def __init__(
			self,
			allele: str | list | None = None,
			usecols: list[str] | None = None,
			dtype: dict | None = None,
			sep: str = "\t"
	) -> None:
		self._delimiter = sep
		self.__allele = allele
		self.__usecols = usecols
		self.__dtype = dtype

		# self._full_data = None
		self.__snp_data: pd.DataFrame | None = None
		self.__header = {}
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

			if self.__allele is not None and self.__usecols is not None:
				raise Exception("Error. Usecols is used for allele is none.")

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

			# # Processing report file
			self.__handler_header(file_rep)
			self.__handler_data(file_rep)

			if not self.__snp_data.empty and self._map_rn is not None:
				self.__snp_data['Sample ID'] = \
					self.__snp_data['Sample ID'].map(
						dict(zip(self._map_rn.SID, self._map_rn.UNIQ_KEY))
					)

		except Exception as e:
			raise e

		return True

	def __handler_header(self, file_rep: Path) -> None:
		""" Processes data from a file, selects meta-information.

		:param file_rep: path, pointer to the file to be read.
		"""

		with open(file_rep, 'r') as file:

			for line in file:
				if self.__class__.__PATTERN_DATA.findall(line.strip()):
					return

				if self.__class__.__PATTERN_HEADER.findall(line.strip()) or\
					len(line.strip()) == 0:
					continue

				key = line.strip().split("\t")[0]
				value = line.strip().split("\t")[1]

				self.__header[key] = value

	def __handler_data(self, file_rep: Path) -> None:
		""" Processes data and forms an array for further processing.

		:param file_rep: path, pointer to the file to be read.
		"""

		with open(file_rep, 'r') as file:

			# Search for the data start index and skip
			for line in file:
				if self.__class__.__PATTERN_DATA.findall(line.strip()):
					break

			# line column
			orig_name_col = file.readline().strip().split(self._delimiter)

			if self.__allele is None and self.__usecols is None:
				self.__snp_data = pd.read_csv(
					file,
					sep=self._delimiter,
					header=None,
					names=orig_name_col,
					dtype=self.__dtype,
					low_memory=True,
					na_filter=True
				)

				return

			sub_n_col = self.__processing_columns(orig_name_col)
			self.__snp_data = pd.read_csv(
				file,
				sep=self._delimiter,
				header=None,
				names=orig_name_col,
				usecols=sub_n_col,
				dtype=self.__dtype,
				low_memory=True,
				na_filter=True
			)

			return

	def __processing_columns(self, lst_col: list[str]) -> list[str] | None:
		""" Processing the line with all the names of the fields and the
		sample of them.

		:param lst_col: List of all fields.
		:return: Returns a tuple with a list of names of selected fields.
		"""

		if self.__usecols is not None:
			check_n_col = [
				item for item in self.__usecols if item in lst_col
			]

			# Check on empty list
			if check_n_col:
				return self.__usecols

			raise Exception(
				f"Error. The USECOLS list contains not true fields."
			)

		# processing alleles
		sample_n_col = self.__sample_by_allele(lst_col)
		if sample_n_col is None:
			raise Exception(
				f"Error. Allele {self.__allele} not in data."
			)

		return sample_n_col

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

	@staticmethod
	def _check_on_ru_symbols(seq: pd.Series) -> bool | None:
		""" Checial verification of the Cyrillic

		:param seq: Squeezed for verification.
		:return: Truth if there are no symbols of Cyril and there is a lie if
			there is.
		"""

		return seq.apply(
			lambda x: bool(re.search('[а-яА-Я]', x)) if x is not nan else x
		).any()
