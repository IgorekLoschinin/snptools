#!/usr/bin/env python
# coding: utf-8

from . import PATH_DIR_FILES
from .._snp import Snp
from snpTools._finalreport import FinalReport


class TestSNP(object):

	def test_snp_process(self) -> None:

		self._fl = FinalReport()
		self._snp = Snp()

		self._fl.handle(
			PATH_DIR_FILES / "fsnp/file1.txt",
			PATH_DIR_FILES / "fsnp/file1.xlsx"
		)
		assert self._snp.process(self._fl.snp_data)
		assert not self._snp.data.empty and self._snp.data.SNP.isin([
			'02011015010000500', '01110152120222512'
		]).all()

		assert not self._snp.c_r.empty and self._snp.c_r.CALL_RATE.round(6).isin(
			[0.882353, 0.882353]
		).all()

		print()
