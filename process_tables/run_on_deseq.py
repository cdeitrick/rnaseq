"""
	This script combines the deseq tables with annotations from the reference genome.
"""

from pathlib import Path
import argparse
from typing import *
import pandas


def combine_deseq_tables(filenames: List[Path]) -> pandas.DataFrame:
	tables = list()
	for filename in filenames:
		strain = filename.stem.split('.')[0]
		t = pandas.read_csv(filename, sep = "\t")
		t['strain'] = strain
		tables.append(t)

	df = pandas.concat(tables)
	df.index.name = "locusTag"

	return df


def get_significant_locus_tags(filename_expression: Path) -> List[str]:
	table = pandas.read_csv(filename_expression, sep = "\t")
	table = table[(table['log2FoldChange'].abs() > 2) & (table['padj'] < 0.05) & (table['chromosome'] != 'Plasmid1')]
	significant_locus_tags = table['locusTag'].tolist()

	return significant_locus_tags


def transform_table(table: pandas.DataFrame) -> pandas.DataFrame:
	strains = table['strain'].unique()
	tables = list()
	for strain in strains:
		strain_table = table[table['strain'] == strain]
		strain_table = strain_table.sort_values(by = 'locusTag').set_index('locusTag')
		strain_table = strain_table[['log2FoldChange', 'padj']]
		strain_table.columns = [f"{strain}_expression", f"{strain}_pvalue"]
		tables.append(strain_table)
	# Should make sure the tables have the same index
	for table in tables:
		assert list(table.index) == list(tables[0].index)
	df = pandas.concat(tables, axis = 1)
	return df


def create_parser(args: List[str] = None) -> argparse.Namespace:
	parser = argparse.ArgumentParser()
	parser.add_argument(
		"--folder",
		help = "The folder with the pairwise deseq results.",
		type = Path
	)

	parser.add_argument(
		"--output",
		help = "The folder to generate the new tables in.",
		type = Path
	)

	parser.add_argument(
		"--annotations",
		help = "The annotations table.",
		type = Path
	)

	if args:
		args = [str(i) for i in args]
		args = parser.parse_args(args)
	else:
		args = parser.parse_args()
	return args


def main():

	#project_folder = Path.home() / "storage" / "projects" / "tils" / "RNAseq" / "differential_expression_trimmed"
	project_folder = Path.home() / "storage" / "projects" / "eisha_rna"
	args = [
		"--folder", project_folder / "results" / "deseq_results",
		"--output", project_folder / "data",
		"--annotations", project_folder / "data" / "annotations.tsv"
	]
	#args = None
	args = create_parser(args)

	filename_expression = args.output / "differential_expression.tsv"
	filename_annotations = args.annotations
	filename_full = filename_expression.with_suffix('.bystrain.all.tsv')
	filename_significant = filename_expression.with_suffix('.bystrain.significant.tsv')
	annotation_columns = ['chromosome', 'name', 'locusTag', 'annotation']

	table_annotations = pandas.read_csv(filename_annotations, sep = "\t")
	table_annotations = table_annotations[annotation_columns]

	filenames = list(args.folder.iterdir())
	filenames = [i for i in filenames if i.stem.split('.')[1] == 'WT']

	table_expression = combine_deseq_tables(filenames)
	table_expression = table_expression.merge(table_annotations, left_on = 'locusTag', right_on = 'locusTag',
		how = 'left')
	table_expression.to_csv(filename_expression, sep = "\t", index = False)

	table_full = transform_table(table_expression)
	table_full = table_full.merge(table_annotations, left_on = "locusTag", right_on = 'locusTag', how = 'left')
	table_full.to_csv(filename_full, sep = "\t", index = False)

	significant_locus_tags = get_significant_locus_tags(filename_expression)

	table_significant = table_expression[table_expression['locusTag'].isin(significant_locus_tags)]
	table_significant = transform_table(table_significant)
	table_significant = table_significant.merge(table_annotations, left_on = "locusTag", right_on = 'locusTag',
		how = 'left')
	table_significant.to_csv(filename_significant, sep = "\t", index = False)


if __name__ == "__main__":

	main()
