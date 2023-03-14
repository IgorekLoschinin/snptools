#!/usr/bin/python

from pathlib import Path
from typing import Self


class FinalReport(object):

	def __init__(self, allele: str = "AB") -> None:
		self.__header = dict()
		self.__data = list()
		self.__allele = allele  # {'Forward', 'Top', 'AB'}

	def read(self, file_rep: Path, sep: str = "\t") -> bool | Self:
		pass

	def write(self, file_rep: Path, ext: str = "txt") -> bool | Self:
		pass

	def handler_header(self) -> dict | None:
		pass

	def handler_data(self) -> list | None:
		pass
