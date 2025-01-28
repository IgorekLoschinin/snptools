#!/usr/bin/env python
# coding: utf-8
__author__ = "Igor Loschinin (igor.loschinin@gmail.com)"

from snptools.src.snplib.statistics import minor_allele_freq as maf

import pytest


class TestMinorAlleleFreq(object):

	@pytest.mark.parametrize("value, res", [
		(0.0, 0.0), (0.9, 0.1), (0.889, 0.111),
		(0.714, 0.286), (0.22, 0.22), (0.45, 0.45), (0.6, 0.4)
	])
	def test_minor_allele_freq(self, value: float, res: float) -> None:
		assert maf(value) == res
