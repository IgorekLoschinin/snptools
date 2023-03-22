#!/usr/bin/env python
# coding: utf-8

from finalreport import FinalReport


class Snp(object):

	def __init__(self) -> None:
		pass

	def handle(self) -> None:
		pass

	def call_rate(self) -> None:
		pass

	def marker_call_rate(self) -> None:
		pass

	def maf(self) -> None:
		pass

	def hwe(self) -> None:
		"""
		pip install snphwe https://github.com/jeremymcrae/snphwe
		"""

	def allele_freq(self) -> None:
		pass

	def verification_parent(self) -> None:
		pass


	# def linkage_disequilibrium(self):
	# 	df = self.replace(5, np.nan).T.astype(np.float16)
	# 	return df.corr().applymap(lambda x: round(x * x, 3))