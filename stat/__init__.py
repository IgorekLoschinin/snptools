#!/usr/bin/env python
# coding: utf-8

from ._snphwe import hwe, hwe_test
from ._callrate import call_rate as cr
from ._freq import allele_freq, minor_allele_freq


__all__ = [
	"cr",
	"allele_freq",
	"minor_allele_freq",
	"hwe",
	"hwe_test"
]
