from pathlib import Path
from typing import *
import pandas
import math
from pprint import pprint

column_label_pvalue = 'padj'

def get_filenames(expression_folder: Path) -> Dict[Tuple[str,str], Path]:
	folder = expression_folder / "results_deseq"
	files = list(i for i in folder.iterdir() if i.suffix == '.csv')

	data = dict()
	for filename in files:
		left, right = filename.stem.split('_')
		data[left, right] = filename

	return data


def normalize_series(series: pandas.Series) -> pandas.Series:
	mean = series.mean()
	sig = series.std()
	result = series.apply(lambda s: (s - mean) / sig)
	return result


def read_table(filename: Path) -> pandas.DataFrame:
	table = pandas.read_csv(filename)
	table.columns = ['locusTag'] + list(table.columns)[1:]

	#table = table.set_index('locusTag')
	table['isSignificant'] = table['pvalue'].apply(lambda s: s < 0.05)
	return table


def volcanoplotplt(table: pandas.DataFrame):
	column_label_fold_change = "log2FoldChange"
	column_label_pvalue = "padj"

	import matplotlib.pyplot as plt
	fig, ax = plt.subplots(figsize = (10, 10))

	ax.scatter(table[column_label_fold_change].tolist(), table[column_label_pvalue].tolist())
	plt.show()





def main():
	project_folder = Path.home() / "storage" / "projects" / "tils"
	expression_folder = project_folder / "RNAseq" / "differential_expression2"

	filenames = get_filenames(expression_folder)

	pprint(filenames)
	table = read_table(filenames['WT', 'A244T'])

	volcanoplot(table)


if __name__ == "__main__":
	main()
