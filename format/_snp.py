# !/usr/bin/env python
# coding: utf-8

__author__ = "Igor Loschinin (igor.loschinin@gmail.com)"

from pathlib import Path
from snpTools import _FIELDS_ILLUMIN, _MAP_FIELDS, _ALLELE_CODE

import pandas as pd


class Snp(object):
	""" The process of converting genomic map data - FinalReport.txt obtained
	from Illumin. Recoding allele data into quantitative data, saving in the
	format necessary for calculating gblup on blupf90. """

	_FIELDS = ['SNP_NAME', 'SAMPLE_ID', 'SNP']
	_F_DTYPE = dict(zip(_FIELDS, (str for _ in range(len(_FIELDS)))))

	def __init__(self, fmt: str | None = "uga") -> None:
		"""
		:param fmt: - Data format to use snp in plink and blupf90. Default
			value "uga".
		"""
		self._format_data = fmt
		self.__data_snp = None

	@property
	def data(self) -> pd.DataFrame:
		return self.__data_snp

	def process(self, data: pd.DataFrame) -> bool:
		""" Data processing and formatting. Calculation of statistical
		information

		:param data: - Data from FinalReport file. Example:
			SNP Name  Sample ID  Allele1 - AB  Allele2 - AB  GC Score  GT Score
			ABCA12	14814	A	A	0.4048	0.8164
			ARS-BFGL-BAC-13031	14814	B	B	0.9083	0.8712
			ARS-BFGL-BAC-13039	14814	A	A	0.9005	0.9096
			ARS-BFGL-BAC-13049	14814	A	B	0.9295	0.8926

		:return: - Returns true if the data was formatted successfully and
		statistical information was calculated, false if an error.
		"""

		if not all(list(map(lambda x: x in data.columns, _FIELDS_ILLUMIN))):
			return False

		self.__data_snp = data.rename(columns=_MAP_FIELDS)
		self.__data_snp['SNP'] = \
			self.__data_snp[['ALLELE1', 'ALLELE2']].\
			sum(axis=1).\
			map(_ALLELE_CODE)

		self.__data_snp = self.__data_snp[Snp._FIELDS].astype(Snp._F_DTYPE)

		if self._format_data is not None and self._format_data == "uga":
			self.__data_snp = self._format_uga(
				self.__data_snp[['SAMPLE_ID', 'SNP']]
			)

	@staticmethod
	def _format_uga(data: pd.DataFrame) -> pd.DataFrame:
		""" Data format to use snp in plink and blupf90 """

		return data.groupby(by='SAMPLE_ID').sum().reset_index()

	def to_file(self, file_path: str | Path) -> None:
		""" Saving data to a file

		:param file_path: - Path to file
		"""

		if isinstance(file_path, str):
			file_path = Path(file_path)

		if self._format_data is not None and self._format_data == "uga":

			max_len = self.__data_snp["SAMPLE_ID"].str.len().max()
			self.__data_snp["SAMPLE_ID"] = self.__data_snp["SAMPLE_ID"].\
				apply(lambda x: "".join([x, " " * (max_len - len(x))]))

			self.__data_snp.\
				apply(" ".join, axis=1).\
				to_csv(file_path, sep=" ", index=False, header=False)

			self.__data_snp["SAMPLE_ID"] = \
				self.__data_snp["SAMPLE_ID"].str.strip()

			return None

		self.__data_snp.\
			to_csv(file_path, sep=" ", index=False, header=False)
