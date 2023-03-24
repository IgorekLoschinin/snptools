#!/usr/bin/env python
# coding: utf-8

from pathlib import Path
from finalreport import FinalReport

DATA_FILE = (
	(Path("./tests/files/fr/data_50k.txt"), Path("./tests/files/fr/data_50k.xlsx")),
	(Path("./tests/files/fr/data_LD7k.txt"), Path("./tests/files/fr/data_LD7k.xlsx")),
)


def test_final_report():
	obj_fr = FinalReport()

	for data_file, conv_file in DATA_FILE:
		assert obj_fr.handle(data_file, conv_file)
		assert obj_fr.snp_data is not None
