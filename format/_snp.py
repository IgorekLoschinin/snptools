#!/usr/bin/env python
# coding: utf-8

__author__ = "Igor Loschinin (igor.loschinin@gmail.com)"

from pathlib import Path

import pandas as pd


_FIELDS_ILLUMIN = ['SNP Name', 'Sample ID', 'Allele1 - AB', 'Allele2 - AB']
_RENAME_FIELDS = ['SNP_NAME', 'SAMPLE_ID', 'ALLELE1', 'ALLELE2']
_MAP_FIELDS = dict(zip(_FIELDS_ILLUMIN, _RENAME_FIELDS))


class Snp(object):
	""" The process of converting genomic map data - FinalReport.txt obtained
	from Illumin. Recoding allele data into quantitative data, saving in the
	format necessary for calculating gblup on blupf90. """

	_ALLELE_CODE = {
		'AA': 0, 'AB': 1, 'BA': 1, 'BB': 2, '--': 5
	}

	def __init__(self) -> None:
		self.__data_snp = None

		self.__data_call_rate = None

	@property
	def data(self) -> pd.DataFrame:
		return self.__data_snp

	@property
	def c_r(self) -> pd.DataFrame:
		return self.__data_call_rate

	def process(self, data: pd.DataFrame) -> bool:
		""" Data processing and formatting. Calculation of statistical
		information.

		:param data: - Data from FinalReport file. Example:
			SNP Name  Sample ID  Allele1 - AB  Allele2 - AB  GC Score  GT Score
			ABCA12	14814	A	A	0.4048	0.8164
			ARS-BFGL-BAC-13031	14814	B	B	0.9083	0.8712
			ARS-BFGL-BAC-13039	14814	A	A	0.9005	0.9096
			ARS-BFGL-BAC-13049	14814	A	B	0.9295	0.8926

		:return: - Returns true if the data was formatted successfully and
		statistical information was calculated, false if an error.
		"""
		try:
			data_snp = data.rename(columns=_MAP_FIELDS)

			self.__data_snp = data_snp[['SNP_NAME', 'SAMPLE_ID']]
			self.__data_snp['SNP'] = \
				data_snp[['ALLELE1', 'ALLELE2']].\
				sum(axis=1).\
				map(self.__class__._ALLELE_CODE).\
				astype(str)

		except Exception as e:
			return False

		return True

	def to_uga(self, file: Path = None) -> bool:
		...

		# self.__data_snp = \
		# 	self.__data_snp.groupby(by='SAMPLE_ID').sum().reset_index()
		# .apply(" ".join, axis=1).to_csv(
		# 	"./safasdfassdfasd.txt", index=False, header=False)