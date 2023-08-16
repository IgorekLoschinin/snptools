#!/usr/bin/env python
# coding: utf-8

from pathlib import Path

import pandas as pd


def isag_disc() -> pd.DataFrame:
	return pd.read_pickle(Path(__file__).parent.joinpath("isag_disc.pl"))


def isag_verif() -> pd.DataFrame:
	return pd.read_pickle(Path(__file__).parent.joinpath("isag_verif.pl"))
