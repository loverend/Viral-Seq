library(tidyverse)
`%notin%` <- Negate(`%in%`)
library(data.table)
library(ggplot2)
library(ggforce)
library(Gviz)
library(GenomicRanges)
library(GenomicAlignments)
library(foreach)
library(doParallel)
library(data.table)
library(edgeR)
library(limma)
library(Glimma)
library(org.Mm.eg.db)
library(gplots)
library(RColorBrewer)
library(cowplot)

#outputdir <- outputdir <- '/well/immune-rep/shared/MISEQ/VIRAL_SEQ/VT_FIRST/BULK_CMV/'

get_counts <- function(outputdir, type_run){
	paths <- list.files(paste0(outputdir, "/ViRNA_SEQ_UNIQUE_Mapping_Analysis"), recursive=TRUE, full.name=TRUE)
	paths <- grep("COUNTS_SUMMARY", paths, value=TRUE)
	#create a list of dataframes 
	df_list <- lapply(paths, fread, header = FALSE, sep="\t")
	# identify the sample they came from 
	d <- basename(paths) 
	d <- gsub(".txt", "", d)
	d <- gsub("VIRNAseq_COUNTS_SUMMARY_", "", d)
	d <- gsub("_RNA-Seq_R1_001", "", d)
	d <- gsub("Homo_sapiens", "HS", d)
	d <- gsub("Aspergillus_fumigatus", "AF", d)
	d <- gsub("Human_betaherpesvirus_5", "CMV", d)
	d <- gsub("Donor", "", d)
	d <- gsub("plus", "", d)
	s <- str_split_fixed(d, "_", 3)
	d <- s[,3]

	# Name each dataframe with the run and filename
	names(df_list) <- d
	# Create combined dataframe  
	df <- rbindlist(df_list, idcol=TRUE)
	df$V1=NULL 
	colnames(df) <- c("Sample", "Geneid", "Chr", "Start", "End", "Strand", "Length", "count", "VirusFilter")
	## Annotate Human chromosomes counts  
	human_chromosomes <- paste(0:22)
	human_chromosomes <- c(human_chromosomes, "X", "Y")
	df$Chr[df$Chr %in% human_chromosomes] <- paste0("Chr", df$Chr[df$Chr %in% human_chromosomes])
	human_chromosomes <- paste0("Chr", human_chromosomes)
	viral_genomes <- unique(df$Chr[df$Chr %like% "NC_"])
	
	## Make wide format
	df_wide <- dcast(df, Sample ~ Chr, value.var = "count")
	viral_counts <- data.frame(df_wide)
	rownames(viral_counts) <- df_wide$Sample 
	viral_counts$Sample <- NULL
	lib_size <- data.frame(rowSums(viral_counts))
	colnames(lib_size) <- "Counts"
	lib_size$Logcounts <- log10(lib_size$Counts)
	
	human_counts <- viral_counts[,c(!colnames(viral_counts) %like% "NC_")] 
	viral_counts_virus <- viral_counts[, c(colnames(viral_counts) %like% "NC_")] 
	### THIS IS DATA FRAME OF FINAL COUNTS (post FEATURECOUNTS)
	write.table(viral_counts, paste0(outputdir, "FINAL_COUNTS_WIDE_", type_run, ".txt"), sep="\t")
	write.table(viral_counts_virus, paste0(outputdir, "FINAL_COUNTS_WIDE_VIRUS_", type_run, ".txt"), sep="\t")
	write.table(human_counts, paste0(outputdir, "FINAL_COUNTS_WIDE_HUMAN_", type_run, ".txt"), sep="\t")
	
	### Normalise the counts
	df_normalise <- data.frame(t(df_wide))
	colnames(df_normalise) <- df_normalise[1,]
	df_normalise <- df_normalise[-1,]
	mx_normalise <- as.matrix(df_normalise)
	storage.mode(mx_normalise) <- "numeric"
	countdata <- mx_normalise
	y <- DGEList(mx_normalise)
	myCPM <- cpm(mx_normalise)
	human_myCPM <- myCPM[!rownames(myCPM) %like% "NC_",] 
	viral_myCPM <- myCPM[rownames(myCPM) %like% "NC_",] 
	

	## Lets look at the normalised counts 
	pdf(paste0(outputdir, "CPM_summary_", type_run, ".pdf"), width=15, height=17)
	plot_form <- reshape2::melt(myCPM)
	colnames(plot_form) <- c("Var", "Sample", "CPM")
	plot_formv <- plot_form[plot_form$Var %in% viral_genomes,]
	plot_form2 <- reshape2::melt(mx_normalise)
	colnames(plot_form2) <- c("Var1", "Sample1", "RAW_COUNT")
	plot_form2v <- plot_form2[plot_form2$Var %in% viral_genomes,]
	## Include only viruses which were detected
	s <- aggregate(plot_form2v$RAW_COUNT, by=list(plot_form2v$Var1), FUN=sum)
	s <- s[s$x > 0,]
	found_viruses <- as.character(s$Group.1)
	new_plot <- cbind(plot_form, plot_form2)
	new_plotv <- cbind(plot_formv, plot_form2v)
	new_plotv <- new_plotv[new_plotv$Var %in% found_viruses,]
	p1 <- ggplot(new_plot, aes(color=Sample, x=RAW_COUNT, y=CPM)) + geom_point()  + theme_classic() + xlab("Raw Counts") + ylab("Counts Per Million")+ggtitle(paste0("Human + Virus: ", type_run)) +labs(colour="Sample")  
	p2 <- ggplot(new_plotv, aes(color=Var, x=RAW_COUNT, y=CPM)) + geom_point()  + theme_classic() + xlab("Raw Counts") + ylab("Counts Per Million")+ggtitle(paste0("Virus: ", type_run)) +labs(colour="Virus")
	p3 <- ggplot(lib_size, aes(x=Counts)) + geom_histogram(fill="White", colour="black", bins = 50) + theme_bw() +xlab("Library Size") +ylab("Count")+ggtitle(paste0("Human + Virus: ", type_run))
	plot(plot_grid(p1, p2, p3, ncol=1))
	dev.off()
	
	myCPM <- t(myCPM)
	human_myCPM  <- t(human_myCPM)
	viral_myCPM <- t(viral_myCPM)
	
	write.table(myCPM, paste0(outputdir, "FINAL_CPM_WIDE_", type_run, ".txt"), sep="\t")
	write.table(human_myCPM, paste0(outputdir, "FINAL_CPM_WIDE_HUMAN_", type_run, ".txt"), sep="\t")
	write.table(viral_myCPM, paste0(outputdir, "FINAL_CPM_WIDE_VIRUS_", type_run, ".txt"), sep="\t")
	
}
