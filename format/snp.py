#!/usr/bin/env python
# coding: utf-8

import pandas as pd

from format import _MAP_FIELDS


class Snp(object):
	""" The process of converting genomic map data - FinalReport.txt obtained
	from Illumin. Recoding allele data into quantitative data, saving in the
	format necessary for calculating gblup on blupf90. """

	_MAP_ALLELE = {
		'AA': 0, 'AB': 1, 'BA': 1, 'BB': 2, '--': 5
	}

	def __init__(self) -> None:
		self.__data_snp = pd.DataFrame()

		self.__data_call_rate = pd.DataFrame()
		self.__data_maf = pd.DataFrame()
		self.__data_hwe = pd.DataFrame()

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

			self.__data_snp['SAMPLE_ID'] = data_snp.SAMPLE_ID
			self.__data_snp['SNP'] = \
				data_snp[['ALLELE1', 'ALLELE2']].\
				sum(axis=1).\
				map(self.__class__._MAP_ALLELE).\
				astype(str)

			self.__data_snp = \
				self.__data_snp.groupby(by='SAMPLE_ID').sum().reset_index()

		except Exception as e:
			return False

		if not self.__data_snp.empty:

			self.__call_rate()
			if self.__data_call_rate.empty:
				return False

			# self.__minor_alllele_freq()
			# if self.__data_maf.empty:
			# 	return False
			#
			# self. __hardy_weinberg_equilibrium()
			# if self.__data_hwe.empty:
			# 	return False

		return True

	@property
	def data(self) -> pd.DataFrame:
		return self.__data_snp

	@property
	def c_r(self) -> pd.DataFrame:
		return self.__data_call_rate

	@property
	def maf(self) -> None:
		return None

	@property
	def hwe(self) -> None:
		return None

	def __call_rate(self) -> None:
		""" The call rate for a given SNP is defined as the proportion of
		individuals in the study for which the corresponding SNP information
		is not missing. In the following example, we filter using a call rate
		of 95%, meaning we retain SNPs for which there is less than 5% missing
		data. """

		self.__data_call_rate['UNIQ_KEY'] = self.__data_snp.SAMPLE_ID
		self.__data_call_rate['CALL_RATE'] = self.__data_snp.SNP.apply(
			lambda x: 1 - (x.count('5') / len(x))
		)

	def __minor_alllele_freq(self) -> None:
		""" The minor allele frequency is therefore the frequency at which the
		minor allele occurs within a population. """

	def __hardy_weinberg_equilibrium(self) -> None:
		"""
		pip install snphwe https://github.com/jeremymcrae/snphwe
		"""

	def __allele_freq(self) -> None:
		pass

	def __verification_parent(self) -> None:
		pass

	# def linkage_disequilibrium(self):
	# https://mr-dictionary.mrcieu.ac.uk/term/ld/
	# 	df = self.replace(5, np.nan).T.astype(np.float16)
	# 	return df.corr().applymap(lambda x: round(x * x, 3))
