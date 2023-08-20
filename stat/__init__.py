#!/usr/bin/env python
# coding: utf-8

from ._snphwe import hwe, hwe_test
from ._callrate import call_rate as cr
from ._freq import allele_freq, minor_allele_freq as maf


__all__ = [
	"cr",
	"allele_freq",
	"maf",
	"hwe",
	"hwe_test"
]
