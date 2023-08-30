#!/usr/bin/env python
# coding: utf-8

from . import DIR_FILES
from .. import make_lgen

import pytest
import pandas as pd


@pytest.fixture
def data_lgen() -> pd.DataFrame:
	return pd.read_csv(DIR_FILES / "fplink/lgen/file.csv")


class TestPlinkFormatLgen(object):

	def test_lgen(self, data_lgen) -> None:
		f = make_lgen(data_lgen)

		print()
