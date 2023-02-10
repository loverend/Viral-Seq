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

#outputdir <- outputdir <- '/gpfs2/well/immune-rep/shared/MISEQ/VIRAL_SEQ/HUMAN_FIRST/BULK_CMV_UNMAPPED/'
#outputdir <- '/well/immune-rep/shared/MISEQ/VIRAL_SEQ/SEPSIS_GAINS_1/'

get_qc <- function(outputdir, type_run, infection=NA){
	paths <- list.files(paste0(outputdir, "/ViRNA_SEQ_UNIQUE_Mapping_Analysis"), recursive=TRUE, full.name=TRUE)
	paths <- grep("unfiltered", paths, value=TRUE)
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
	d <- gsub("_R1_001", "", d)
	d <- gsub("_RNA-Seq", "", d)
	d <- gsub(".Unmapped.out", "", d)
	d <- gsub("QC_unfiltered", "", d)
	if(outputdir %like% "CMV"){
		s <- str_split_fixed(d, "_", 3)
		d <- s[,3]
	}

	# Name each dataframe with the run and filename
	names(df_list) <- d
	# Create combined dataframe  
	df <- rbindlist(df_list, idcol=TRUE, fill=TRUE)
	colnames(df) <- c("Sample", "Virus", "N_reads", "N_unique_reads", "Percent_uniquely_mapped", "Mean_read_quality", "Sd_read_quality", "A", "C", "G", "T", "Sequence_entropy", "Spatial_distribution", "Longest_contig", "DUST_score", "Percent_high_quality_reads", "Coverage_per_read_Ratio", "Mean_Read_Length", "Coverage_Percent_1000bp", "genome", "PassedFiltering")
	
	df <- df[df$Virus!="N_mapped_reads",]
	df <- df[!df$N_unique_reads=="0",]
	
	df$VIRUS_USE <- "NONE" 
	if(!is.na(infection)){
		if(infection=="CMV"){
			df$VIRUS_USE[df$Sample %like% "CMV"] <- "POSITIVE"
			df$VIRUS_USE[!df$Sample %like% "CMV"] <- "NEGATIVE"
		}
		
		if(infection=="EBV"){
			df$VIRUS_USE[df$Sample %like% "MUC2258" | df$Sample %like% "MUC2266" | df$Sample %like% "MUC2274"] <- "NEGATIVE"
			df$VIRUS_USE[!df$VIRUS_USE %like% "NEGATIVE"] <- "POSITIVE"
		}
	}
	
	
	numeric_cols <- sapply(df, Hmisc::all.is.numeric)
	numeric_cols <- names(numeric_cols[numeric_cols=="TRUE"])
	df <- data.frame(df)
	df[numeric_cols] <- sapply(df[numeric_cols], as.numeric)

	## Annotate Human chromosomes counts  
	### THIS IS DATA FRAME OF FINAL COUNTS (post FEATURECOUNTS)
	write.table(df, paste0(outputdir, "QC_measures_", type_run, ".txt"), sep="\t")
	
	## Plot summary Per Virus 
	pdf(paste0(outputdir, "QC_summary_", type_run, ".pdf"), width=15, height=15)
	for (i in 1:length(unique(df$Virus))){
		QC_result <- df[df$Virus==unique(df$Virus)[i],]
		Mean_mapping_length <- 101
		minread <- 50
		name_virus <- Virus_database$Virus_name[Virus_database$viral_genome==unique(df$Virus)[i]]
		p1 <- ggplot(QC_result, aes(color=VIRUS_USE, x=N_reads, y=(N_unique_reads))) + geom_point()  + theme_classic() + xlab("Number of Mapped Reads") + ylab("Number Uniquely Mapped Reads")  + ggtitle(name_virus)+ expand_limits(x = 0, y = 0)+ guides(color="none")
		p2 <- ggplot(QC_result, aes(color=VIRUS_USE, x=N_unique_reads, y=(Spatial_distribution*100))) + geom_point()  + theme_classic() + xlab("Number of Uniquely Mapped Reads") + ylab("% Mapped genome") + ylim(0, 100) + ggtitle(name_virus)+ guides(color="none")
		p3 <- ggplot(QC_result, aes(color=VIRUS_USE, x=Sequence_entropy, y=(Spatial_distribution*100))) + geom_point()  + theme_classic() + xlab("Sequence Entropy") + ylab("% Mapped genome")  + ylim(0, 100) + ggtitle(name_virus)+ guides(color="none")
		p4 <- ggplot(QC_result, aes(color=VIRUS_USE, x=Longest_contig, y=DUST_score)) + geom_point()  + theme_classic() + xlab("Longest Contig (nt)") + ylab("DUST Score") + ggtitle(name_virus)+ guides(color="none")
		p5 <- ggplot(QC_result, aes(color=VIRUS_USE, x=Mean_read_quality, y=Sd_read_quality)) + geom_point()  + theme_classic() + xlab("Mean Read Quality") + ylab("SD Read Quality") + ggtitle(name_virus)+ guides(color="none")
		p6 <- ggplot(QC_result, aes(color=VIRUS_USE, x=Longest_contig, y=Coverage_per_read_Ratio)) + geom_point()  + theme_classic() + xlab("Longest Contig") + ylab("Coverage per Read Ratio") + ggtitle(name_virus) + guides(color="none") +geom_vline(xintercept=100, col="blue") +geom_vline(xintercept=200, col="green")
		p7 <- ggplot(QC_result, aes(color=VIRUS_USE, x=(Spatial_distribution*100), y=Coverage_per_read_Ratio)) + geom_point() + theme_classic() + xlab("% Mapped genome") + ylab("Coverage per Read Ratio") + ggtitle(name_virus) + xlim(0, 100) + guides(color="none")
		p8 <- ggplot(QC_result, aes(color=VIRUS_USE, y=Coverage_Percent_1000bp, x=N_unique_reads)) + geom_point()  + theme_classic() + ylab("% Coverage for 1000bp with highest coverage") + xlab("Number Uniquely Mapped Reads") + ggtitle(name_virus)+ guides(color="none")
		if(length(rownames(QC_result))>=1){
			QC_result <- QC_result[QC_result$N_unique_reads>=1, ]
			if(length(rownames(QC_result))>1){
				QC_scale <- data.frame(QC_result)
				rownames(QC_scale) <- QC_scale$Sample
				QC_scale <- QC_scale[, c(3:17, 19)]
				cols_discard <- apply(QC_scale, 2, function(x) length(unique(x)))
				cols_discard <- names(cols_discard[cols_discard<2])
				cols_discard2 <- lapply(lapply(QC_scale, is.na), sum)
				cols_discard2 <- names(cols_discard2)[cols_discard2 >0]
				cols_discard <- unique(c(cols_discard, cols_discard2))
				QC_scale <- QC_scale[, !colnames(QC_scale) %in% cols_discard]
				QC_scale <- sapply(QC_scale, as.numeric, 2)
				rownames(QC_scale) <-  QC_result$Sample
				QC_scaled <- scale(QC_scale)
				QC_pca <- prcomp(QC_scaled, center = TRUE)
				eigencoordinates <- data.frame(QC_pca$x)
				if(!is.na(infection)){
					if(infection=="CMV"){
						eigencoordinates$CMV <- NA 
						eigencoordinates$CMV[rownames(eigencoordinates) %like% "CMV"] <- "POSITIVE"
						eigencoordinates$CMV[!rownames(eigencoordinates) %like% "CMV"] <- "NEGATIVE"
					}
					if(infection=="EBV"){
						eigencoordinates$CMV <- NA 
						eigencoordinates$CMV[rownames(eigencoordinates) %like% "MUC2258" | rownames(eigencoordinates) %like% "MUC2266" | rownames(eigencoordinates) %like% "MUC2274"] <- "NEGATIVE"
						eigencoordinates$CMV[!rownames(eigencoordinates) %like% "MUC2258" & !rownames(eigencoordinates) %like% "MUC2266" & !rownames(eigencoordinates) %like% "MUC2274"] <- "POSITIVE"
					}
				}
				if(is.na(infection)){
					eigencoordinates$CMV <- "NONE"
				}
				
				p9 <- ggplot(eigencoordinates, aes( x=PC1, y=PC2, color=CMV)) + geom_point()  + theme_classic() + ylab("PC2") + xlab("PC1") + ggtitle(name_virus)+ guides(color="none")
				plot(plot_grid( p1, p2, p3, p4, p5, p6, p7, p8, p9,  ncol=3, labels = "AUTO"))
				} else {
			plot(plot_grid( p1, p2, p3, p4, p5, p6, p7, p8,  ncol=3, labels = "AUTO"))
			}
		}
	}
	dev.off()	
}
