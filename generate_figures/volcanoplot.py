from pathlib import Path
from typing import *
import pandas
import math
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
import statsmodels.api as sm
from loguru import logger
from statsmodels.regression import linear_model
from statsmodels.sandbox.stats.multicomp import TukeyHSDResults  # Used to add a typing annotation to tukeyhsd()
from statsmodels.stats.multicomp import MultiComparison


def assign_colors():
	pass


def volcanoplot(table: pandas.DataFrame):
	colors = ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628", "#f781bf"]
	strains = ["WT", "A244T", "N274Y", "N445K", "P421L", "S3C", "T41T", "tRNA"]
	palette = dict(zip(strains, colors))
	palette['not significant'] = "#777777"
	color = [(strain if (pvalue < 0.05 and abs(fold) > 2) else "not significant") for strain, pvalue, fold in
		zip(table['strain'].tolist(), table['padj'].tolist(), table['log2FoldChange'].tolist())]
	table['color'] = color

	column_label_fold_change = "log2FoldChange"
	import plotly.express as px
	fig = px.scatter(
		table,
		x = column_label_fold_change, y = "-log10(pvalue)",
		hover_data = ['locusTag'],
		color = 'color',
		color_discrete_map = palette,
		title = "Gene Expression Relative to WT",
		width = 1920,
		height = 1080
	)
	# fig.write_image("volcanoplot.png")
	fig.show()


def boxplot(table: pandas.DataFrame):
	# plt.style.use('fivethirtyeight')
	fig, ax = plt.subplots(figsize = (12, 10))
	fig: plt.Figure

	ax: plt.Axes = table.boxplot('log2FoldChange', by = 'chromosome', ax = ax)
	ax.set_title("")
	ax.set_ylabel("log2FoldChange")
	ax.set_facecolor("#FFFFFF")
	ax.grid(False)
	plt.tight_layout()


def main():
	project_folder = Path.home() / "storage" / "projects" / "eisha_rna"
	project_folder = Path.home() / "storage" / "projects" / "tils" / "RNAseq" / "differential_expression_trimmed"
	results_folder = project_folder / "results"
	data_folder = project_folder / "data"
	filename_annotations = project_folder / "data"  / "annotations.tsv"
	filename_table = data_folder / "differential_expression.tsv"

	table = pandas.read_csv(filename_table, sep = "\t")
	annotation_table = pandas.read_csv(filename_annotations, sep = "\t")

	print(table.info())
	print()
	print(annotation_table.info())
	table['-log10(pvalue)'] = table['padj'].apply(lambda s: -math.log10(s))
	table = table.merge(annotation_table, how = "left", left_on = 'locusTag', right_on = 'locusTag')
	# table = table[table['seq_type'] == 'plasmid']
	volcanoplot(table)

	table = table[table['padj'] < 0.05]
	table = table[table['log2FoldChange'].abs() > 2]

	#plt.savefig("chromosome_expression.png")
	plt.show()



if __name__ == "__main__":
	main()
