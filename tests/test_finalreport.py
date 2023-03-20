#!/usr/bin/env python
# coding: utf-8
import time
from functools import wraps

from . import DATA_FILE, NUM_CONV_FILE
from finalreport import FinalReport


def timeit(func):
	@wraps(func)
	def timeit_wrapper(*args, **kwargs):
		start_time = time.perf_counter()
		result = func(*args, **kwargs)
		end_time = time.perf_counter()
		total_time = end_time - start_time
		print(f'Function {func.__name__} Took {total_time:.10} seconds')
		return result

	return timeit_wrapper


@timeit
def test_final_report():
	obj_fr = FinalReport()
	assert obj_fr.handle(DATA_FILE, NUM_CONV_FILE)
