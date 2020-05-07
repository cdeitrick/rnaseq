from pathlib import Path
from typing import *
import pandas
import matplotlib.pyplot as plt
import seaborn

def main():
	#folder = Path.home() / "storage" / "projects" / "eisha_rna"
	project_folder = Path.home() / "storage" / "projects" / "tils" / "RNAseq" / "differential_expression"
	filename_expression = project_folder / "data" / "differential_expression.tsv"
	filename_counts = project_folder / "data" / "abundance.all.tsv"
	filename_annotations = project_folder / "data" / "annotations.tsv"

	table_expression = pandas.read_csv(filename_expression, sep = "\t")
	table_expression = table_expression[(table_expression['padj'] < 0.05) & (table_expression['log2FoldChange'].abs() > 2)]
	significant_locus_tags = table_expression['locusTag'].tolist()
	table_expression = pandas.read_csv(filename_expression, sep = "\t")
	table_expression = table_expression[table_expression['locusTag'].isin(significant_locus_tags)]
	table_expression.to_csv(
		project_folder / "data" / "differential_expression.significant.tsv",
		sep = "\t", index = False
	)

	table_counts = pandas.read_csv(filename_counts, sep = "\t")

	table_annotations = pandas.read_csv(filename_annotations, sep = "\t")

	#table_expression = table_expression.merge(table_annotations, left_on = 'locusTag', right_on = 'locusTag')

	unique_chromosomes = table_annotations['chromosome'].unique()
	unique_strains = table_expression['strain'].unique()

	#table_counts = table_counts.merge(table_annotations, left_on = "locusTag", right_on = "locusTag")
	#table_counts = table_counts[table_counts['locusTag'].isin(table_expression['locusTag'].tolist())]
	#print(table_counts.info())
	#fig, ax = plt.subplots(figsize = (12,10))
	"""
	seaborn.boxplot(
		x = 'strain',
		y = 'counts',
		hue = 'chromosome',
		data = table_counts,
		ax = ax
	)

	seaborn.boxplot(
		x = 'strain',
		y = 'log2FoldChange',
		hue = 'chromosome',
		data = table_expression,
		ax = ax
	)
	"""

	#plt.show()



if __name__ == "__main__":
	main()