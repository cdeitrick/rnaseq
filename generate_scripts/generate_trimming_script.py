from pathlib import Path
from typing import *

def get_command(filename:Path, output:Path):
	adaptors_filename = "/home/cld100/Documents/github/workflows/pipelines/resources/adapters.fa"
	command = [
		"trimmomatic", "SE", "-phred33",
		filename, output,
		f"ILLUMINACLIP:{adaptors_filename}:2:30:10",
		f"LEADING:3",
		f"TRAILING:3",
		f"SLIDINGWINDOW:4:15",
		f"MINLEN:36"
	]

	command = [str(i) for i in command]
	command = " ".join(command)
	return command

def checkdir(path:Path)->Path:
	path = Path(path)
	if not path.exists():
		path.mkdir()
	return path

def main():
	project_folder = Path.home() / "storage" / "projects" / "tils" / "RNAseq"
	sample_folder = project_folder / "samples"
	results_folder = checkdir(project_folder / "trimmed")
	filename_script = project_folder / "trimming_workflow.sh"

	filenames = sample_folder.glob("**/*.fastq.gz")
	commands = list()
	for filename in filenames:
		sample_folder = filename.parent.name
		name = filename.stem + '.fastq'
		output_folder = checkdir(results_folder / sample_folder)
		print(output_folder.exists(), output_folder)
		output_filename = output_folder / name

		command = get_command(filename, output_filename)
		commands.append(command)

	content = "\n".join(commands)
	filename_script.write_text(content)

if __name__ == "__main__":
	main()