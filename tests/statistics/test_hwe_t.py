#!/usr/bin/env python
# coding: utf-8
__author__ = "Igor Loschinin (igor.loschinin@gmail.com)"

from snptools.src.snplib.statistics import hwe_test

import pytest
import pandas as pd


class TestHWE(object):

	@pytest.mark.parametrize(
		"seq, freq",
		[
			('2212120', 0.714),
			('02011015010000500', 0.2),
			('01110152120222512', 0.6),
			('00005005005', 0.0),
			('22521212222', 0.9),
			('52221521222', 0.889)
		]
	)
	def test_hweT_true(self, seq: str, freq: float) -> None:
		"""
		 check snphwe gives expected p-values
		"""

		_seq_snp = pd.Series(list(map(int, seq)))

		assert hwe_test(_seq_snp, freq)

	@pytest.mark.parametrize("seq, freq", [('000000000102', 0.125)])
	def test_hweT_false(self, seq: str, freq: float) -> None:
		"""
		 check snphwe gives expected p-values
		"""

		_seq_snp = pd.Series(list(map(int, seq)))

		assert not hwe_test(_seq_snp, freq)
