from pathlib import Path
from typing import *
import pandas
import matplotlib.pyplot as plt
import seaborn
from loguru import logger

TableType = pandas.DataFrame  # So that the typing annotations are shorter


def anova():
	pass


def make_boxplot(table_expression: pandas.DataFrame):
	fig, ax = plt.subplots(figsize = (12, 10))
	ax = seaborn.boxplot(x = 'strain', y = 'log2FoldChange', hue = 'chromosome', data = table_expression, ax = ax)
	ax.set_title("Gene Expression by Chromosome/Plasmid")
	plt.tight_layout()
	plt.savefig("gene_expression_by_strain_and_chromosome.png")
	plt.show()


def get_strain_count_columns(labels: List[str], strain) -> List[str]:
	result = [i for i in labels if (strain in i and 'counts' in i)]
	return result


def read_table_counts(filename: Path) -> pandas.DataFrame:
	""" Reads in the counts table and renames the columns to be more consistant with the other tables."""
	table_counts = pandas.read_csv(filename, sep = "\t")

	# Sort by name so the columns are in a consistent order
	table_counts = table_counts[sorted(table_counts.columns)]
	# table_counts = table_counts.set_index(['strain', 'replicate', 'locusTag'])
	return table_counts


def read_table_expression(filename: Path) -> pandas.DataFrame:
	""" Reads the expression table and applies a few filters. """
	table_expression = pandas.read_csv(filename, sep = "\t")  # .set_index('locusTag')
	# Filter by significance and the amount expressed.
	table_expression_filtered = table_expression[table_expression['padj'] < 0.05]
	table_expression_filtered = table_expression_filtered[table_expression_filtered['log2FoldChange'] > 2]
	logger.debug(f"Length of the filtered expression table: {len(table_expression_filtered)}")

	# Need to include all locus tags that were significant in at least one sample
	# So the `table_expression_filtered` table has all of the significant locus tags, but only
	# for the strains they were significant in.
	significant_locus_tags = table_expression_filtered['locusTag'].tolist()
	table_expression = table_expression[table_expression['locusTag'].isin(significant_locus_tags)]
	logger.debug(f"Significant Expression table: {len(table_expression)}")
	# Add annotations to the epression table so we can select rows based on which chromosome a gene is on.
	# table_expression = table_expression.merge(table_annotations, how = "left", left_index = True,
	#	right_index = True)

	# Remove columns that aren't needed

	# 'locusTag' should already be the index
	# expression_columns = ['log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj', 'strain', 'chromosome']
	keep_columns = ['log2FoldChange', 'strain', 'padj', 'locusTag']
	table_expression = table_expression[keep_columns]

	# May as well sort the columns
	table_expression = table_expression[sorted(table_expression.columns)]

	return table_expression


def add_id_to_label(strain: str, replicate: int, columns: List[str]) -> List[str]:
	""" Renames the labels so they are unique to a specific strain and replicate. """
	cols = [f"{strain}_{replicate}_{i}" for i in columns]
	return cols

def validate_columns(table:pandas.DataFrame, expected:List[str]):
	try:
		assert sorted(table.columns) == expected
	except AssertionError:
		message = "The columns in the table do not math the expected labels."
		logger.error(message)
		logger.error(f"\tExpected: {expected}")
		logger.error(f"\tActual: {sorted(table.columns)}")
		raise AssertionError(message)
class AnnotateExpression:
	""" Generates a table combining the gene expression table with the read count table. """

	def __init__(self):
		pass

	def merge_tables(self, table_expression: TableType, table_counts: TableType, label_wt: str,
			label_other: str) -> TableType:
		""" Adds read count information for the locus tags in the expression table.
			Notes
			-----
			The expression table should have these columns:
			- `log2FoldChange`
			- `strain`
			- `padj`
			The index should be named `locusTag`.
			Both counts tables should have these columns:
			- `strain`
			- `replicate`
			- `locusTag`
			- `length`
			- `counts`
			- `sampleId`
			The index for this table should not be set.
		"""

		# Make sure that the tables are properly formatted.
		columns_expression = ['log2FoldChange', 'padj', 'strain']
		columns_counts = ['counts', 'length', 'locusTag', 'replicate', 'sampleId', 'strain', 'targetId']
		validate_columns(table_expression, columns_expression)
		validate_columns(table_counts, columns_counts)
		assert table_expression.index.name == 'locusTag'

		# Filter out locus tags in the counts table that are not present in the expression table.
		locus_tags = list(table_expression.index)
		table_counts = table_counts[table_counts['locusTag'].isin(locus_tags)]

		# Split the caounts table by strain
		table_counts_wt, table_counts_other = self.split_table_counts(table_counts, label_wt, label_other)

		# The count tables have three replicates per sample, so add those and the mean value for that sample.
		replicates = table_counts['replicate'].unique()

		for replicate in replicates:
			counts_replicate_wt = table_counts_wt[table_counts_wt['replicate'] == replicate]
			counts_replicate_other = table_counts_other[table_counts_other['replicate'] == replicate]

			# Make sure the index for both counts tables is set as `locusTag`
			counts_replicate_wt = counts_replicate_wt.set_index('locusTag')
			counts_replicate_other = counts_replicate_other.set_index('locusTag')

			# Remove columns that aren't needed anymore
			counts_replicate_wt = counts_replicate_wt[['length', 'counts']]
			counts_replicate_other = counts_replicate_other[['length', 'counts']]

			# Rename the columns so they are uniqur to the sample and replicate.
			counts_replicate_wt.columns = add_id_to_label(label_wt, replicate, counts_replicate_wt.columns)
			counts_replicate_other.columns = add_id_to_label(label_other, replicate, counts_replicate_other.columns)

			# Merge the tables
		#	print('Before merging: ', len(table_expression))
			table_expression = table_expression.merge(counts_replicate_wt, left_index = True, right_index = True,
				how = 'left')
			#print('During merging: ', len(table_expression))
			table_expression = table_expression.merge(counts_replicate_other, left_index = True, right_index = True,
				how = 'left')
			#print('After merging: ', len(table_expression))

		# Add the mean of the counts across replicates to the table

		replicate_columns_wt = get_strain_count_columns(table_expression.columns, label_wt)
		replicate_columns_other = get_strain_count_columns(table_expression.columns, label_other)

		counts_wt = table_expression[replicate_columns_wt]
		count_other = table_expression[replicate_columns_other]

		table_expression[f"{label_wt}_{label_other}_counts"] = counts_wt.mean(axis = 1)
		table_expression[f"{label_other}_counts"] = count_other.mean(axis = 1)


		table_expression = table_expression
		return table_expression

	@staticmethod
	def split_table_counts(table_counts: TableType, label_wt: str, label_other: str) -> \
			Tuple[TableType, TableType]:
		""" Splits the counts table by strain and filters out locus tags not in the expression table.
			The counts table should be indexed yb strain, replicate, locus_tag
		"""

		table_counts_wt = table_counts[table_counts['strain'] == label_wt]
		table_counts_other = table_counts[table_counts['strain'] == label_other]

		return table_counts_wt, table_counts_other

	def run(self, filename_expression: Path, filename_counts: Path, filename_annotations: Path):
		label_wt = "WT"
		# label_strain = "A244T"

		table_expression = read_table_expression(filename_expression)
		table_counts = read_table_counts(filename_counts)
		strains = table_expression['strain'].unique()
		# table_annotations = pandas.read_csv(filename_annotations, sep = "\t").set_index('locusTag')
		logger.debug(strains)
		strain_tables = list()
		for label_strain in strains:
			table_expression_strain = table_expression[table_expression['strain'] == label_strain]
			table_expression_strain = table_expression_strain.set_index('locusTag')

			# We only need the counts for WT and the current strain.
			# We can also omit any locus tags not in the expression table.
			table_expression_strain = self.merge_tables(table_expression_strain, table_counts, label_wt, label_strain)
			strain_tables.append(table_expression_strain.reset_index())

		first_table = strain_tables[0]

		for t in strain_tables:
			first_table = first_table.merge(t, left_on = "locusTag", right_on = "locusTag")

		return first_table


def main():
	project_folder = Path.home() / "storage" / "projects" / "tils" / "RNAseq" / "differential_expression"
	project_folder = Path.home() / "storage" / "projects" / "eisha_rna"
	filename_expression = project_folder / "data" / "differential_expression.tsv"
	filename_abundance = project_folder / "data" / "abundance.all.tsv"
	filename_annotations = project_folder / "data" / "annotations.tsv"
	filename_result = project_folder / "data" / "rna_expression.tsv"

	app = AnnotateExpression()
	_table_expression = app.run(filename_expression, filename_abundance, filename_annotations)
	print(_table_expression.info())
	_table_expression.to_csv(filename_result, sep = "\t")


if __name__ == "__main__":
	main()
