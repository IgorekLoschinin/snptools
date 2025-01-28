#!/usr/bin/env python
# coding: utf-8
__author__ = "Igor Loschinin (igor.loschinin@gmail.com)"

from . import DIR_FILES
from snptools.src.snplib.format import make_map

import pytest
import pandas as pd


@pytest.fixture
def data_map() -> pd.DataFrame:
	return pd.read_csv(DIR_FILES / "fplink/map/file_bovinesnp50.csv")


class TestPlinkFormatMap(object):

	def test_map_true(self, data_map) -> None:

		res = make_map(data_map)
		assert not res.empty

	def test_map_raise(self, data_map) -> None:
		with pytest.raises(
				KeyError, match="Manifest has no data to build map format!"
		):
			make_map(data_map)
			make_map(pd.DataFrame())
			make_map(
				pd.DataFrame(columns=['Chr', 'Name', 'MapInfo', 'morgans'])
			)

		with pytest.raises(
				KeyError, match="Manifest has no data to build map format!"
		):
			make_map(pd.DataFrame())

	def test_map_empty(self) -> None:
		assert make_map(
				pd.DataFrame(columns=['Chr', 'Name', 'MapInfo', 'morgans'])
			).empty
