 
library("tximport")
library("readr")
library("tximportData")
library("rjson")
library("DESeq2")
library("dplyr")
library("rhdf5")
"
	The `samples` table should have metadata about the individual samples.
		Columns:
		- pop
		- center
		- assay
		- sample
		- experiment
		- run
		- condition
"

install_dependencies <- function() {
	# Installs the dependencies needed for the analysis.
	if (!requireNamespace("BiocManager", quietly = TRUE))
		install.packages("BiocManager")
	
	BiocManager::install("tximport")
	BiocManager::install("tximportData")
	#BiocManager::install("readr")
	#BiocManager::install("rjson")
	BiocManager::install("rDESeq2")
	BiocManager::install("rhdf5")
	
}


get_filenames_tables_count <- function(folder) {
	regex <- paste(folder, "**/*.tsv", sep = "/")
	filenames <- Sys.glob(regex)
	sample_names <- get_sample_names(filenames)
	names(filenames) <- sample_names
	
	return(filenames)
}

get_sample_name <- function(filename) {

	parts <- strsplit(filename, "/")[[1]]
	sample_name <- parts[length(parts)-1]
	return(sample_name)
}

get_sample_names <- function(filenames) {
	sample_names <- c()
	for (filename in filenames){
		sample_name <- get_sample_name(filename[[1]])
		sample_names <- append(sample_names, c(sample_name))
	}
	return(sample_names)
}

read_gene_map <- function(filename) {
	# Reads in the transcript id to gene table. Implemented as a standalone function in case the table needs to be modified (cleaned) before returning it.
	# This currently uses `locus_tag` as the geneid, but that may change later.
	table <- read.csv(filename, sep = "\t")
	# Remove any extra columns
	reduced_table = select(table, "target_id", "locus_tag")
	# The columns should be named "TXNAME" and "GENEID" so that tximport can parse it.
	names(reduced_table) <- c("TXNAME", "GENEID")
	return(reduced_table)
}

read_sample_table <- function(filename) {
	table <- read.csv(filename, sep = "\t")

	return(table)
}

generate_sample_table <- function(filenames_abundance) {
	# 	samples[,c("pop","center","run","condition")]
	# x <- data.frame("SN" = 1:2, "Age" = c(21,15), "Name" = c("John","Dora"))
	# Get the condition from the filename.
	condition <- names(filenames_abundance)
	condition <- gsub("_1", "", condition)
	condition <- gsub("_2", "", condition)
	condition <- gsub("_3", "", condition)
	df <- data.frame(
		"pop" = "Population1",
		"center" = "Pitt",
		"run" = names(filenames_abundance),
		"condition" = condition
		)
	rownames(df) <- df$run
	
	return(df)
}


prepare_data <- function(folder_abundance, filename_gene_map) {
	
	filenames_abundance <- get_filenames_tables_count(folder_abundance)
	
	sample_table <- generate_sample_table(filenames_abundance)
	table_genemap = read_gene_map(filename_gene_map)
	
	txi <- tximport(filenames_abundance, type="kallisto", tx2gene=table_genemap)
	
	ddsTxi <- DESeqDataSetFromTximport(
		txi, 
		colData = sample_table, 
		design = ~ condition
	)
	
	ddsTxi$condition <- relevel(ddsTxi$condition, ref = "WT")
	
	return(ddsTxi)
}

generate_output_filename <- function(output_folder, left, right) {
	basename <- paste(left, right, sep = "_")
	basename <- paste(basename, "csv", sep = ".")
	filename <- paste(output_folder, basename, sep = "/")
	return(filename)
}

generate_treatment_keys <- function(elements, pairwise=FALSE, reference="WT"){
	if (pairwise){
		treatment_matrix <- combn(elements, 2)
		treatments <- data.frame("left" = treatment_matrix[1,], "right" = treatment_matrix[2,])
	}
	else{
		treatments <- data.frame("left" = reference, "right" = elements)
	}
	return(treatments)
}
sort_by_reference <- function(array, reference) {
	# Orders the array so that the reference value is the first element.
	# Mainly used so that the output files put the reference label on the left.
	array_1 <- array[array != reference]
	array_1 <- c(c(reference), array_1)
	return (array_1)
}
run_deseq <- function(ddsTxi, reference, output_folder, pairwise = FALSE) {
	if (!file.exists(output_folder)){
		dir.create(output_folder)
	}
	ddsTxi$condition <- relevel(ddsTxi$condition, ref = reference)
	dds <- DESeq(ddsTxi)
	
	conditions <- as.character(unique(ddsTxi$condition))

	conditions <- sort_by_reference(conditions, reference)

	treatments <- generate_treatment_keys(conditions, pairwise)

	for (treatment_index in 1:dim(treatments)[1]) {
		left <- treatments[treatment_index, 1]
		right <- treatments[treatment_index, 2]
		
		left <- toString(left)
		right <- toString(right)
		if (left == right){
			next
		}
		message <- paste("Analyzing:", left, right, sep = " ")
		print(message)
		treatment_result <- results(dds, contrast = c("condition", left, right))
		filename <- generate_output_filename(output_folder, left, right)
		write.csv(treatment_result, filename)
	}
}
#############################################################################################################################################################
########################################################## Entry Point for the script #######################################################################
#############################################################################################################################################################

# Need the folder where kallisto generated all of the `abundance.tsv` tables for each sample.
folder_kallisto_output <- "/media/cld100/FA86364B863608A1/Users/cld100/Storage/projects/tils/RNAseq/differential_expression/results_kallisto"
# Where to put the output tables.
folder_output <- "/media/cld100/FA86364B863608A1/Users/cld100/Storage/projects/tils/RNAseq/differential_expression/results_deseq"
# This table maps transcript ids from kallisto to gene identifiers based on the same reference file used with kallisto.
filename_gene_map <- "/media/cld100/FA86364B863608A1/Users/cld100/Storage/projects/tils/RNAseq/differential_expression/gene_map.tsv"
reference <- "WT"
# Whether to compare all-vs-all (pairwise == TRUE) or just WT vs everthing else (pairwise == False)	
pairwise <- FALSE


prepared_data = prepare_data(folder_kallisto_output, filename_gene_map)
run_deseq(prepared_data,reference, folder_output, pairwise = pairwise)

