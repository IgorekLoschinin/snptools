#!/usr/bin/env python
# coding: utf-8
__author__ = "Igor Loschinin (igor.loschinin@gmail.com)"

from ._snphwe import hwe, hwe_test
from ._callrate import call_rate
from ._freq import allele_freq, minor_allele_freq


__all__ = [
	"call_rate",
	"allele_freq",
	"minor_allele_freq",
	"hwe",
	"hwe_test"
]
