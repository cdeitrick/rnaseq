	library("DESeq2")
	
	read_table_counts <- function(filename){
		table_counts <- read.csv(filename, sep = "\t", row.names = 'locusTag')
		# Need to convert the table counts to integers because they may be float with zeros after the decimal point.	
		#table_counts[] <- lapply(table_counts, function(x) as.numeric(x))
	
		return(table_counts)
	}
	
	read_table_annotations <- function(filename){
		table_annotations <- read.csv(filename, sep = "\t")
	
		row.names(table_annotations) <- table_annotations$locusTag
		
		return(table_annotations)
	}
	
	read_table_design <- function(filename) {
		table_design <- read.csv(filename, sep = "\t")
		row.names(table_design) <- table_design$sampleId
		
		return(table_design)
	}
	get_sample_labels <- function(labels){
		#split_labels <- strsplit(labels, "_")
		split_labels <-c()
		
		for (label in labels){
			split_label <- strsplit(label, "_")
			value = split_label <- split_label[[1]][[1]]
			
			split_labels <- c(split_labels, value)
		}
		split_labels <- unique(split_labels)
		return(split_labels)
	}
	save_deseq_results <- function(results, filename) {
		write.table(as.data.frame(results), filename, sep = "\t", row.names = TRUE)
	}
	run_deseq_pairwise <- function(deseq_object, output_folder, labels) {
		for (left in labels) {
			for (right in labels){
				if (left!=right){
					deseq_results_strain <- results(deseq_object, contrast=c("strain",left, right))
					
					filename = paste(output_folder,"/", left,".",right, ".differentialexpression.tsv", sep = "")
					print(filename)
					save_deseq_results(deseq_results_strain, filename)
				}
			}
		}
	}
	run_deseq <- function(table_counts, table_design, table_annotations, output_folder){
		print("\tGenerating the Deseq matrix object...")
		deseq_matrix <- DESeqDataSetFromMatrix(
			countData = table_counts,
			colData = table_design,
			design = ~ strain	 + strain:type# Formula summarizing the experimental design
		)
		
		# Make sure the WT is considered the reference sample.
		deseq_matrix$strain <- relevel(deseq_matrix$strain, ref = "WT")
		#deseq_matrix$ppc <- relevel(deseq_matrix$group, ref = 'group')
		print("\tGenerating the deseq object...")
		deseq_object <- DESeq(deseq_matrix)
	
		#strain_order <- c("WT","A244T","N274Y","N445K", "P421L", "S3C", "T41T", "tRNA")
		#strain_order <- get_sample_labels(table_design$strain)
		print("Running pairwise workflow...")
		deseq_results <- run_deseq_pairwise(deseq_object, output_folder, table_design$strain)
		print("\tGenerating results...")
	
		filename_deseq_results <- paste(output_folder, "ppc_vs_tils.differentialexpression.tsv", sep  = "/")
		print("\tSaving results...")
		save_deseq_results(deseq_results, filename_deseq_results)
		# Write the deseq results to a table.
		newlist <- list(deseq_object, deseq_results)
	}
	
	
	
	generate_figure_clusters <- function(deseq_object, filename) {
		deseq_object <- estimateSizeFactors (deseq_object)
		deseq_rlog <- rlog ( deseq_object , blind = TRUE )
		rlog.norm.counts <- assay (deseq_rlog)
		distance.m_rlog <- as.dist(1 - cor( rlog.norm.counts, method = "pearson"))
	}
	
	generate_figure_pca <- function(deseq_object) {
		deseq_object <- estimateSizeFactors (deseq_object)
		deseq_rlog <- rlog ( deseq_object , blind = TRUE )
		P <- plotPCA(deseq_rlog, intgroup = "strain")
		P + theme_bw() + ggtitle ( " Rlog transformed counts " )
	}
	
	generate_report <- function(deseq_object, table_counts, table_design){
		library("ReportingTools")
		
		desReport <- HTMLReport(
			shortName = 'RNAseq_analysis_with_DESeq',
			title = 'RNA-seq analysis of differential expression using DESeq',
			reportDirectory = "./reports"
		)
		publish(
			deseq_object,
			desReport,
			countTable=table_counts, 
			pvalueCutoff=0.05,
			factor = table_design$strain,
			#conditions = table_annotations$strain,
			#conditions=conditions,
			#.modifyDF=makeDESeqDF,
			expName="deseq",
			reportDir="./reports"
		)
		finish(desReport)
	}
	
	run_pca_interactive <- function(deseq_object, table_design){
		library(pcaExplorer)
		pcaExplorer(dds = deseq_object, coldata = table_design)
	}
	
	main <- function(){
		project_folder <- "/home/cld100/storage/projects/tils/RNAseq/differential_expression_trimmed"
		#project_folder <- "/home/cld100/storage/projects/eisha_rna/"
		data_folder <- paste(project_folder, "data", sep = "/")
		results_folder <- paste(project_folder, "results", sep = "/")
		
		
		filename_counts <- paste(data_folder, "abundance.all.matrix.tsv", sep = "/")
		filename_design <- paste(data_folder, "abundance.all.design.tsv", sep = "/")
		filename_annotations <- paste(data_folder, "annotations.tsv", sep = "/")
		
		
			# Specify where the deseq tables should be saved
		folder_deseq <- paste(results_folder, "deseq_results", sep = "/")
		
		folder_figures <- paste(project_folder, "figures", sep = "/")
		filename_clusters <- paste(folder_figures, "clusters.histogram.png", sep = "/")
		
		dir.create(file.path(folder_deseq), showWarnings = FALSE)
		
		print("Reading in the source tables...")
		# Matrix of read counts
		print("Reading in the count matrix...")
		table_counts <- read_table_counts(filename_counts)
		print("Reading the design table...")
		table_design <- read_table_design(filename_design)
		print("Reading the annotations table...")
		table_annotations <- read_table_annotations(filename_annotations)
		print("Running deseq...")
		deseq_tables <- run_deseq(table_counts, table_design, table_annotations, folder_deseq)
		deseq_object <- deseq_tables[[1]]
		deseq_results <- deseq_tables[[2]]
		print("Generating the cluster image...")
		generate_figure_clusters(deseq_object, filename_clusters)
		
		#resLFC <- lfcShrink(deseq_object, coef="strain_A244T_vs_WT", type="apeglm")
		print("Generating report...")
		generate_report(deseq_object, table_counts, table_design)
		print("Generating PCA figure...")
		#generate_figure_pca(deseq_object)
		
		newlist <- list(table_counts, table_design, table_annotations, deseq_object, deseq_results)
		return(newlist)
	}
								
	results <- main()
	table_annotations <- results[[3]]
	deseq_object <- results[[4]]
	table_counts <- results[[1]]
	table_design <- results[[2]]
	
	#library(pcaExplorer)
	#pcaExplorer(dds = deseq_object, coldata = table_design)
	
	
