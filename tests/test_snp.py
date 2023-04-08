#!/usr/bin/env python
# coding: utf-8

from pathlib import Path

from . import PATH_DIR_FILES
from ..format.snp import Snp
from ..finalreport import FinalReport


class TestSNP(object):
	DATA_FILE = Path(f"{PATH_DIR_FILES}/fsnp/file1.txt")
	CONV_FILE = Path(f"{PATH_DIR_FILES}/fsnp/file1.xlsx")

	def test_snp_process(self) -> None:

		self._fl = FinalReport()
		self._snp = Snp()

		self._fl.handle(
			Path(self.__class__.DATA_FILE), Path(self.__class__.CONV_FILE)
		)
		assert self._snp.process(self._fl.snp_data)
		assert not self._snp.data.empty and self._snp.data.SNP.isin([
			'02011015010000500', '01110152120222512'
		]).all()

		assert not self._snp.c_r.empty and self._snp.c_r.CALL_RATE.round(6).isin(
			[0.882353, 0.882353]
		).all()

		print()
