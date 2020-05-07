"""
	generates a table mapping `target_id` from the abundance table to `protein_id` in the reference cds fasta file.

"""

import re
from pathlib import Path
from typing import *
import argparse
import pandas
from pprint import pprint
from Bio import SeqIO


def parse_description(item: str) -> Dict[str, str]:
	regex = "\[([a-z_]+)[=](.+?)\]"
	matches = re.findall(regex, item)

	return dict(matches)


def generate_annotation_table(filename: Path) -> pandas.DataFrame:
	"""
		Returns a dict with the annotations in the reference file.
		Returns
		-------
		List[Dict[str,str]]
		- `locus_tag`
		- `protein`
		- `protein_id`
		- `location`
		- `gbkey`
		- `name`
	"""
	table = list()
	for record in SeqIO.parse(filename, 'fasta'):
		record_data = parse_description(record.description)
		record_data['name'] = record.name
		table.append(record_data)
	table = pandas.DataFrame(table)
	return table


def merge_annotations(reference_table: pandas.DataFrame, count_table: pandas.DataFrame) -> pandas.DataFrame:
	# Merge tables based on target_id and protein_id

	def clean_label(s) -> str:
		result = s.split('_cds_')[-1]
		result = result.rpartition('_')[0]
		return result

	# count_table['protein_id'] = count_table['target_id'].apply(clean_label)

	merged_table = pandas.merge(reference_table, count_table, left_on = 'name', right_on = 'target_id', how = 'left')

	# Remove unused columns so the file isn't too large.

	keep_columns = ['locus_tag', 'target_id', 'eff_length', 'est_counts']
	merged_table = merged_table[keep_columns]
	merged_table.columns = ['locusTag', 'targetId', 'length', 'counts']
	#merged_table = merged_table.set_index('locusTag')

	return merged_table


def generate_matrix(data: Dict[str, pandas.DataFrame]) -> pandas.DataFrame:
	# Extract the 'counts' value as a pandas.Series object.

	counts = list()
	for sample_name, sample_table in data.items():
		series = sample_table['counts']
		# Make sure the sample name is attached
		series.name = sample_name
		counts.append(series)

	df = pandas.DataFrame(counts).transpose()

	# Convert to integers
	for label in df.columns:
		df[label] = df[label].astype(int)
	return df


def create_parser(args: List[str] = None) -> argparse.Namespace:
	parser = argparse.ArgumentParser()
	parser.add_argument(
		"--reference",
		help = "The reference cds fasta file that was used with kallisto.",
		type = Path
	)
	parser.add_argument(
		"--folder",
		help = "Folder with kallisto results.",
		type = Path
	)

	parser.add_argument(
		"--filename",
		help = "The consolidated table.",
		type = Path
	)
	parser.add_argument(
		"--annotations",
		help = "The annotations.tsv table.",
		type = Path
	)


	if args:
		args = [str(i) for i in args]
		args = parser.parse_args(args)
	else:
		args = parser.parse_args()
	return args

def generate_experimental_design_matrix(samples:List[str])->pandas.DataFrame:
	table = list()
	for label in samples:
		strain, replicate = label.split('_')
		row = {
			'strain': strain,
			'replicate': replicate,
			'sampleId': label
		}
		table.append(row)

	df = pandas.DataFrame(table).set_index('sampleId')
	df['type'] = df['strain'].apply(lambda s: 'ppc' if (s == 'S3C' or s == 'T41T') else 'tils')
	return df


def main(args = None):
	args = create_parser(args)
	filename_matrix = args.filename.with_suffix('.matrix.tsv')
	filename_design = args.filename.with_suffix('.design.tsv')
	filename_annotations = args.annotations

	if args.folder.is_file():
		filenames = [args.table]
	else:
		filenames = list(args.folder.glob("**/abundance.tsv"))

	#annotation_table = generate_annotation_table(args.reference)
	annotation_table = pandas.read_csv(filename_annotations, sep = "\t")
	plasmid_table = annotation_table[annotation_table['chromosome'] == 'Plasmid1']
	plasmid_locus_tags = plasmid_table['locusTag'].tolist()

	annotation_table = generate_annotation_table(args.reference)

	kallisto_tables = dict()
	for filename in filenames:
		sample_name = filename.parent.name
		strain, replicate = sample_name.split('_')
		count_table = pandas.read_csv(filename, sep = "\t")

		fulltable = merge_annotations(annotation_table, count_table)
		fulltable['sampleId'] = sample_name
		fulltable['strain'] = strain
		fulltable['replicate'] = replicate
		fulltable = fulltable.set_index('locusTag')
		kallisto_tables[sample_name] = fulltable

	abundance_table = pandas.concat(kallisto_tables.values())
	abundance_table.info()
	#abundance_table = abundance_table[~abundance_table['locusTag'].isin(plasmid_locus_tags)]
	abundance_table = abundance_table.loc[[i for i in abundance_table.index if i not in plasmid_locus_tags]]
	#abundance_table = abundance_table.set_index('locusTag')
	abundance_table.to_csv(args.filename, sep = "\t", index = True)

	matrix = generate_matrix(kallisto_tables)
	matrix.to_csv(filename_matrix, sep = "\t", index = True)

	table_design = generate_experimental_design_matrix(matrix.columns)
	table_design.to_csv(filename_design, sep = "\t")

if __name__ == "__main__":
	project_folder = Path.home() / "storage" / "projects" / "tils" / "RNAseq" / "differential_expression_trimmed"
	project_folder = Path.home() / "storage" / "projects" / "eisha_rna"
	debug_args = [
		"--reference", project_folder / "data" / "reference" / "GCA_000203955.1_ASM20395v1_cds_from_genomic.fna",
		"--folder", project_folder / "results" / "kallisto_results",
		"--filename", project_folder / "data" / "abundance.all.tsv",
		"--annotations", project_folder / "data" / "annotations.tsv"
	]

	debug_args = None

	main(debug_args)
