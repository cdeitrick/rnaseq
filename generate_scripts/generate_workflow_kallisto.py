"""
	Generates a bash script to run all samples through kallisto
"""
from pathlib import Path
from typing import *
import argparse


def get_command(filename: Path, reference: Path, output_folder: Path) -> str:
	# f"kallisto quant -i {reference} --single -l 180 -s 20 -o {output}"
	#print(filename.name, filename.stem.split('_')[:2])
	sample, replicate = filename.stem.split('_')[:2]
	sample_id = f"{sample}_{replicate}"
	command = [
		"kallisto", "quant",
		"-i", reference,
		"--single",
		"-l 180",
		"-s", "20",
		"-o", output_folder / sample_id,
		filename
	]

	command = [str(i) for i in command]
	return " ".join(command)


def get_index_command(reference: Path, output: Path) ->str:
	# kallisto index --index=HI2424.cds.index reference/GCA_000203955.1_ASM20395v1_cds_from_genomic.fna
	command = [
		"kallisto", "index",
		f"--index={output}",
		reference
	]
	command = [str(i) for i in command]
	return " ".join(command)

def create_parser(args: List[str] = None) -> argparse.Namespace:
	parser = argparse.ArgumentParser()

	parser.add_argument(
		"--samples",
		help = "The folder with the sample fastqs.",
		type = Path,
	)

	parser.add_argument(
		"--filename",
		help = "Path of the bash script that will be created.",
		type = Path
	)
	parser.add_argument(
		"--reference",
		help = "A cds feature table from the ncbi.",
		type = Path
	)
	parser.add_argument(
		"--output",
		help = "The output folder to send the kallisto results to.",
		type = Path
	)
	if args:
		args = [str(i) for i in args]
	args = parser.parse_args(args)
	return args


def main():
	# /media/cld100/FA86364B863608A1/Users/cld100/Storage/projects/eisha_rna/RNA_seq/
	project_folder = Path.home() / "storage" / "projects" / "tils" / "RNAseq"
	args = [
		"--samples", project_folder / "trimmed",
		"--filename", project_folder / "kallisto_workflow.sh",
		"--reference", project_folder / "data" / "reference" / "GCA_000203955.1_ASM20395v1_cds_from_genomic.fna",
		"--output", project_folder / "results" / "results_kallisto"
	]
	args = None
	args = create_parser(args)

	if not args.output.exists():
		args.output.mkdir()

	assert args.reference.exists()


	filenames = list(args.samples.glob("**/*.fastq"))
	filenames += list(args.samples.glob("**/*.fastq.gz"))
	filenames = [i for i in filenames if ('R1' in i.name or 'R2' in i.name)]


	index_filename = args.reference.with_suffix('.index')
	index_command = get_index_command(args.reference, index_filename)
	kallisto_commands = [index_command] + [get_command(f, index_filename, args.output) for f in filenames]

	content = "\n".join(kallisto_commands)

	args.filename.write_text(content)


if __name__ == "__main__":
	main()
