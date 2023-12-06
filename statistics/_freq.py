#!/usr/bin/env python
# coding: utf-8
__author__ = "Igor Loschinin (igor.loschinin@gmail.com)"

import pandas as pd


def allele_freq(
		data: pd.DataFrame | str, id_col: str = None, seq_col: str = None
) -> pd.DataFrame | float | None:
	""" The allele frequency represents the incidence of a gene variant in a
	population.

	:param data: Data array.
	:param id_col: Columns with snp names.
	:param seq_col: Columns with value snp in format ucg - 0, 1, 2, 5.
	:return: Return the alleles frequency.
	"""

	if isinstance(data, pd.DataFrame):
		try:
			if data[seq_col].dtype.hasobject:
				if not data[seq_col].str.isdigit().all():
					return None

				return data.\
					loc[data[seq_col] != "5", [id_col, seq_col]]. \
					groupby(by=id_col)[seq_col]. \
					apply(lambda x: x.astype("int8").sum() / (2 * x.count())).\
					reset_index().\
					round(3)

			return data.\
				loc[data[seq_col] != 5, [id_col, seq_col]].\
				groupby(by=id_col)[seq_col].\
				apply(lambda x: x.sum() / (2 * x.count())).\
				reset_index().\
				round(3)

		except Exception as e:
			raise e

	elif isinstance(data, str):
		if not data.isdigit():
			return None

		sam_seq = tuple(
			map(int, filter(lambda x: x if x != "5" else None, data))
		)
		return round(sum(sam_seq) / (2 * len(sam_seq)), 3)

	else:
		return None


def minor_allele_freq(value: float) -> float:
	""" The minor allele frequency is therefore the frequency at which the
	minor allele occurs within a population.

	:param value: Allele frequency
	:return: Return the minor alleles frequency
	"""

	if value > 0.5:
		return round(1 - value, 3)

	return round(value, 3)
