#!/usr/bin/env python
# coding: utf-8

from pathlib import Path

import pandas as pd


class FinalReport(object):
	"""  """

	__PATTERN_HEADER = "[Header]"
	__PATTERN_DATA = "[Data]"

	def __init__(self, allele: str = "AB", sep: str = "\t") -> None:
		self._delimiter = sep
		self.__full_data = None

		self.__header = dict()
		self.__snp_data = None
		self.__allele = allele  # {'Forward', 'Top', 'AB'}

	@property
	def header(self) -> dict:
		return self.__header

	@property
	def snp_data(self) -> pd.DataFrame | None:
		return self.__snp_data

	def handle(
			self, file_rep: Path | str, conv_file: Path | str
	) -> bool:
		"""  """

		if isinstance(file_rep, str):
			file_rep = Path(file_rep)

		if isinstance(conv_file, str):
			conv_file = Path(conv_file)

		if not self.read(file_rep):
			# logger
			print('Error')
			return False

		try:
			self.__handler_header()
			self.__handler_data()
			self.__convert_num(conv_file)

			return True

		except Exception as e:
			print(f'{e}')
			return False

	def read(self, file_rep: Path) -> bool:
		""" Reading data from the final_report file

		:param file_rep: - path, pointer to the file to be read
		:return: - Returns true if the read was successful, false if it failed
		"""
		try:
			self.__full_data = file_rep.read_text().split("\n")
			return True

		except Exception as e:
			# logger
			return False

	def __handler_header(self) -> None:
		"""  """

		for item in self.__full_data:
			if item == self.__class__.__PATTERN_DATA:
				return

			if item != self.__class__.__PATTERN_HEADER:
				key = item.strip().split("\t")[0]
				value = item.strip().split("\t")[1]

				self.__header[key] = value

	def __handler_data(self) -> None:
		"""  """

		temp = 1
		for item in self.__full_data:
			if item == self.__class__.__PATTERN_DATA:
				break
			temp += 1

		self.__snp_data = pd.DataFrame(
			[
				item_data.split(f"{self._delimiter}")
				for item_data in self.__full_data[temp + 1:]
			],
			columns=self.__full_data[temp].split(f"{self._delimiter}")
		).dropna()

	def __convert_num(self, path_file: Path) -> None:
		"""

		:param path_file:
		:return:
		"""

		map_rn = pd.read_excel(
			path_file, header=None, names=['SID', 'UNIQ_KEY'], dtype={'SID': str}
		).dropna()

		self.__snp_data = self.__snp_data['Sample ID'].map(
			dict(zip(map_rn.SID, map_rn.UNIQ_KEY))
		)

	def to_txt(self, file_rep: Path, ext: str = "txt", sep=" ") -> bool:
		"""  """



if __name__ == "__main__":
	f = FinalReport()
	f.handle("./tests/data/data_50k.txt", "./tests/data/data_50k.xlsx")

