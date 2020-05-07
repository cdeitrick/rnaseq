from pathlib import Path
from typing import *
import pandas
from loguru import logger
import argparse
def read_table(filename:Path)->pandas.DataFrame:

	table = pandas.read_csv(filename, sep = "\t")

	sample_id = filename.parent.name
	print(sample_id)
	try:
		strain, replicate = sample_id.split('_')
	except ValueError:
		replicate = sample_id[-1]
		strain = sample_id[:-1]

	logger.debug(f"{strain}, {replicate}")

	table['strain'] = strain
	table['replicate'] = replicate
	table['sampleId'] = sample_id

	return table




def consolidate_kallisto_results():
	""" Combines all of the abundance tables into a single table."""
	project_folder = Path.home() / "storage" / "projects" / "tils" / "RNAseq" / "differential_expression"
	results_folder = project_folder /"results" / "kallisto-results"
	results_folder = Path("/media/cld100/FA86364B863608A1/Users/cld100/Storage/projects/eisha_rna/RNA_seq/kallisto/")


	table = list()
	filenames = list(results_folder.glob("**/*.genes.tsv"))

	from pprint import pprint
	pprint(filenames)
	for filename in filenames:
		t = read_table(filename)
		table.append(t)

	df = pandas.concat(table)
	print(df.head().to_string())

	df.to_csv("eisha.abundance.genes.all.tsv", sep = "\t", index = False)



if __name__ == "__main__":
	consolidate_kallisto_results()