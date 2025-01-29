#!/usr/bin/env python
# coding: utf-8
__author__ = "Igor Loschinin (igor.loschinin@gmail.com)"

import pandas as pd


def call_rate(
		data: pd.DataFrame | str,
		id_col: str = None,
		snp_col: str = None
) -> pd.DataFrame | float | None:
	""" The call rate for a given SNP is defined as the proportion of
	individuals in the study for which the corresponding SNP information is
	not missing. In the following example, we filter using a call rate of 95%,
	meaning we retain SNPs for which there is less than 5% missing data.

	Of the say, 54K markers in the chip, 50K have been genotyped for a
	particular animal, the “call rate animal” is 50K/54K=93%.

	Of the say, 900 animals genotyped for marker CL635944_160.1, how many
	have actually been successfully read? Assume that 600 have been read, then
	the “call rate marker” is 600/900 = 67%.

	:param data: Pre-processed data on which the call rate is calculated.
	:param id_col: The name of the column with the id of the animals or
		markers.
	:param snp_col: The name of the column with the snp sequence.
	:return: Return dataframe with call rates for each animal if a dataframe
		is transmitted. The number if the snp sequence is passed as a string.
		None if there were errors.
	"""

	if isinstance(data, pd.DataFrame):
		try:
			if data[snp_col].dtype.hasobject:
				if not data[snp_col].str.isdigit().all():
					return None

				return data[[id_col, snp_col]].\
					groupby(by=id_col)[snp_col].\
					apply(lambda x: 1 - ((x == "5").sum() / len(x))).\
					reset_index()

			return data[[id_col, snp_col]]. \
				groupby(by=id_col)[snp_col]. \
				apply(lambda x: 1 - ((x == 5).sum() / len(x))). \
				reset_index()

		except Exception as e:
			raise e

	elif isinstance(data, str):
		if not data.isdigit():
			return None

		return round(1 - (data.count('5') / len(data)), 6)

	else:
		return None
