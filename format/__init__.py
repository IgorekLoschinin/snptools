#!/usr/bin/env python
# coding: utf-8

from ._snp import Snp
from ._plink import (
	make_map,
	make_ped,
	make_lgen,
	make_fam
)

__all__ = [
	"Snp",
	"make_map",
	"make_ped",
	"make_fam",
	"make_lgen"
]
