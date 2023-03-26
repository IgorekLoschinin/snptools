#!/usr/bin/env python
# coding: utf-8

from pathlib import Path
from finalreport import FinalReport
from format.snp import Snp

DATA_FILE = Path("./tests/files/fsnp/FinalReport.txt")
CONV_FILE = Path("./tests/files/fsnp/FinalReport.xlsx")


# def test_final_report():
# 	obj_fr = FinalReport()
# 	obj_fr.handle(DATA_FILE, CONV_FILE)
#
# 	obj_snp = Snp()
# 	assert obj_snp.process(obj_fr.snp_data)
