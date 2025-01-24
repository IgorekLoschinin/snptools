#!/usr/bin/env python
# coding: utf-8
__author__ = "Igor Loschinin (igor.loschinin@gmail.com)"

from snptools.src.snplib.statistics import hwe

import pytest


class TestHWE(object):

	def test_snphwe(self) -> None:
		"""
		 check snphwe gives expected p-values
		"""
		assert hwe(500, 10, 5000) == 0.6515718999145375
		assert hwe(1000, 20, 5000) == 1.2659849194317374e-05

	def test_snphwe_odd_inputs(self) -> None:
		"""
		check snphwe with odd inputs
		"""
		# should raise errors with odd inputs

		with pytest.raises(ValueError, match="snphwe: zero genotypes"):
			hwe(0, 0, 0)

		with pytest.raises(ValueError, match="snphwe: negative allele count"):
			hwe(-5, 10, 1000)

	def test_snphwe_large_input(self) -> None:
		"""
		check snphwe doesn't give errors with large sample sizes
		"""
		assert hwe(200000, 200000, 200000) == 0.0

	def test_snphwe_uncertain_genotypes(self) -> None:
		"""
		check uncertain genotypes give correct p-values
		"""
		assert hwe(4989.99999, 494999.999, 9.9999999) == 0.5702231983054381
