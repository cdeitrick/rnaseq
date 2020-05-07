from pathlib import Path
from typing import *
import pytest
try:
	import analysis
except ModuleNotFoundError:
	from . import analysis

def test_get_strain_count_column():
	columns = ["WT_1_counts", "WT_1_length", "A422T_3_counts", "A422T_2_length"]

	expected = ['WT_1_counts']

	result = analysis.get_strain_count_columns(columns, 'WT')
	assert result == expected


if __name__ == "__main__":
	test_get_strain_count_column()