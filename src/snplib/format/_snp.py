# !/usr/bin/env python
# coding: utf-8

__author__ = "Igor Loschinin (igor.loschinin@gmail.com)"

from pathlib import Path
from .__settings import FIELDS_ILLUMIN, MAP_FIELDS

import pandas as pd


class Snp(object):
	""" The process of converting genomic map data - FinalReport.txt obtained
	from Illumin. Recoding allele data into quantitative data, saving in the
	format necessary for calculating gblup on blupf90.

	:argument fmt: Data format to use snp in plink and blupf90. Default
		value "uga". """

	_ALLELE_CODE = {
		'AA': 0, 'AB': 1, 'BA': 1, 'BB': 2, '--': 5
	}

	_FIELDS = ['SNP_NAME', 'SAMPLE_ID', 'SNP']
	_F_DTYPE = dict(zip(_FIELDS, (str for _ in range(len(_FIELDS)))))

	def __init__(self, fmt: str | None = "uga") -> None:
		self._format_data = fmt
		self.__data_snp = None

	@property
	def data(self) -> pd.DataFrame | None:
		return self.__data_snp

	def process(self, data: pd.DataFrame) -> None:
		""" Data processing and formatting. Calculation of statistical
		information

		:param data: Data from FinalReport file. Example:
			SNP Name  Sample ID  Allele1 - AB  Allele2 - AB  GC Score  GT Score
			ABCA12	14814	A	A	0.4048	0.8164
			ARS-BFGL-BAC-13031	14814	B	B	0.9083	0.8712
			ARS-BFGL-BAC-13039	14814	A	A	0.9005	0.9096
			ARS-BFGL-BAC-13049	14814	A	B	0.9295	0.8926

		:return: Returns true if the data was formatted successfully and
			statistical information was calculated, false if an error.
		"""

		if not all(list(map(lambda x: x in data.columns, FIELDS_ILLUMIN))):
			raise KeyError(
				'The name of the fields does not match the finalreport.txt '
				'file from Illumina'
			)

		self.__data_snp = data.rename(columns=MAP_FIELDS)
		self.__data_snp['SNP'] = \
			self.__data_snp[['ALLELE1', 'ALLELE2']].\
			sum(axis=1).\
			map(Snp._ALLELE_CODE)

		self.__data_snp = self.__data_snp[Snp._FIELDS].astype(Snp._F_DTYPE)

		if self._format_data is not None and self._format_data == "uga":
			self.__data_snp = self._format_uga(
				self.__data_snp[['SAMPLE_ID', 'SNP']]
			)

	@staticmethod
	def _format_uga(data: pd.DataFrame) -> pd.DataFrame:
		""" Data format to use snp in plink and blupf90. """

		return data.groupby(by='SAMPLE_ID').sum().reset_index()

	def to_file(self, file_path: str | Path) -> None:
		""" Saving data to a file.

		:param file_path: Path to file
		"""

		if isinstance(file_path, str):
			file_path = Path(file_path)

		if self._format_data is not None and self._format_data == "uga":

			max_len = self.__data_snp["SAMPLE_ID"].str.len().max()

			self.__data_snp.\
				apply(
					lambda x: " ".join([
						self._add_space(x.iloc[0], max_len), x.iloc[1]
					]),
					axis=1
				).\
				to_csv(file_path, index=False, header=False)

			self.__data_snp["SAMPLE_ID"] = \
				self.__data_snp["SAMPLE_ID"].str.strip()

			return None

		self.__data_snp.to_csv(file_path, sep=" ", index=False)

	@staticmethod
	def _add_space(value: str, max_len: int) -> str:
		""" Adding spaces up to the maximum length of the value in the
		sample_id data.

		:param value: Sample_id value
		:param max_len: Max len sample_id value
		:return: Return replacing value
		"""
		return "".join([value, " " * (max_len - len(value))])
