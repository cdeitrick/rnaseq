from pathlib import Path
from typing import *
import pandas
import math
import itertools
import argparse
import matplotlib.pyplot as plt
import seaborn


def squareform(pairwise_values: Dict[Tuple[str, str], float]) -> pandas.DataFrame:
	""" Converts a dictionary with all pairwise values for a set of points into a square matrix representation.
	"""
	keys = sorted(set(itertools.chain.from_iterable(pairwise_values.keys())))
	_square_map = dict()
	for left in keys:
		series = dict()
		for right in keys:
			value = pairwise_values.get((left, right), 0)
			series[right] = value
		_square_map[left] = series
	return pandas.DataFrame(_square_map)


def calculate_distance(left: pandas.Series, right: pandas.Series):
	elements = [(i - j) ** 2 for i, j in zip(left.values, right.values)]
	elements = [i for i in elements if not math.isnan(i)]
	distance = sum(elements)
	distance = math.sqrt(distance)
	return distance


def generate_expresson_matrix(table_expression: pandas.DataFrame) -> pandas.DataFrame:
	df = table_expression.pivot(index = 'locusTag', columns = 'strain', values = 'log2FoldChange')

	return df


def generate_distance_matrix_from_expression(matrix: pandas.DataFrame):
	"""
		The matrix should have strains/samples as columns and locus tags as the index.
	"""
	items = matrix.columns
	distances = dict()
	for left in items:
		left_series = matrix[left]
		for right in items:
			right_series = matrix[right]
			if left == right: continue
			distance = calculate_distance(left_series, right_series)
			distances[left, right] = distance
	return squareform(distances)


def create_parser(args: List[str] = None) -> argparse.Namespace:
	parser = argparse.ArgumentParser()

	parser.add_argument(
		"--filename",
		help = "The concatenated differential expression table.",
		type = Path
	)

	parser.add_argument(
		"--output",
		help = "Path to the generated figure.",
		type = Path
	)

	if args:
		args = [str(i) for i in args]
		args = parser.parse_args(args)
	else:
		args = parser.parse_args()

	return args


def generate_heatmap(distances: pandas.DataFrame, kind: str):
	fig, ax = plt.subplots(figsize = (12, 12))

	label_size_x = 24
	label_size_y = 24
	if kind == 'strain':
		label_size_ticks = 24
	else:
		label_size_ticks = 16
	add_annotations = kind == 'strain'

	ax = seaborn.heatmap(
		distances, square = True, ax = ax,
		annot = add_annotations,
		fmt = '.0f',
		annot_kws = {'fontSize': 24},
		cmap = "Greens",
		linecolor = "#FFFFFF",
		linewidth = 1
	)
	ax.set_title("Euclidean Distance between samples", fontsize = 28)
	ax.set_xlabel(kind, fontsize = label_size_x)
	ax.set_ylabel(kind, fontsize = label_size_y)

	for tick in ax.xaxis.get_major_ticks():
		tick.label.set_fontsize(label_size_ticks)
	for tick in ax.yaxis.get_major_ticks():
		tick.label.set_fontsize(label_size_ticks)

	plt.show()


def main(args):
	table = pandas.read_csv(args.filename, sep = "\t")

	kind = 'strain' if 'strain' in table.columns else 'sample'

	if kind == 'strain':
		matrix = generate_expresson_matrix(table)
	else:
		matrix = table.set_index('locusTag')
	distances = generate_distance_matrix_from_expression(matrix)

	generate_heatmap(distances, kind)


if __name__ == "__main__":
	project_folder = Path.home() / "storage" / "projects" / "tils" / "RNAseq" / "differential_expression_trimmed"
	data_folder = project_folder / "data"
	filename_table = data_folder / "differential_expression.tsv"
	filename_table = data_folder / "abundance.all.matrix.tsv"

	args = [
		"--filename", filename_table,
		"--output", "heatmap.png"
	]
	args = create_parser(args)
	main(args)
