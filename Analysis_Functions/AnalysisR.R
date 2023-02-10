
module use -a /apps/eb/dev/ivybridge/modules/all
module load R-bundle-Bioconductor/3.11-foss-2020a-R-4.0.0
module load SAMtools/1.10-GCC-9.3.0
R

## Set up environment 
.libPaths(c("/gpfs2/well/immune-rep/users/kvi236/MORE_R_MODULES", .libPaths()))
options(repos='http://cran.ma.imperial.ac.uk/')

## Reading in sepsis samples and forming matrix 
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
library(NMF)
library("scatterplot3d") 

# Working directory
setwd('/well/immune-rep/users/kvi236/VIRUS/Sepsis_Scripts')
outputdir <- ('/well/immune-rep/shared/MISEQ/VIRAL_SEQ/SEPSIS_GAINS_1/')

## Virus Database (where all the name/nc ids are stored
Virus_database = read.delim("/gpfs2/well/immune-rep/shared/CODE/VIRAL_SEQ_Reference/NCBI_Viral_Seq_Reference.txt",header=T,sep="\t")
Virus_database$viral_genome <- Virus_database$Name_sequence
Virus_lengths <- Virus_database[,c("Virus_name", "Genome_length")]
Virus_database <- Virus_database[,c("viral_genome", "Virus_name")]
Virus_database$Virus_name <- gsub("\\|", "", Virus_database$Virus_name)
#-------------------------------------------------------------------------------------------------------------
# READ COUNTING ANALYSIS 
source('/gpfs2/well/immune-rep/shared/CODE/Viral-Seq/Analysis_Functions/CombineCounts.R')

get_counts('/well/immune-rep/shared/MISEQ/VIRAL_SEQ/VT_FIRST/BULK_CMV/', "All_Reads")
get_counts('/gpfs2/well/immune-rep/shared/MISEQ/VIRAL_SEQ/HUMAN_FIRST/BULK_CMV_UNMAPPED/', "Unmapped_Reads")

get_counts('/well/immune-rep/shared/MISEQ/VIRAL_SEQ/VT_FIRST/BULK_EBV_LCL/', "All_Reads")
get_counts('/gpfs2/well/immune-rep/shared/MISEQ/VIRAL_SEQ/HUMAN_FIRST/BULK_EBV_UNMAPPED/', "Unmapped_Reads")

get_counts('/well/immune-rep/shared/MISEQ/VIRAL_SEQ/VT_FIRST/BULK_HEPC/', "All_Reads")
get_counts('/well/immune-rep/shared/MISEQ/VIRAL_SEQ/HUMAN_FIRST/BULK_HEPC_UNMAPPED/', "Unmapped_Reads")

get_counts('/well/immune-rep/shared/MISEQ/VIRAL_SEQ/SEPSIS_GAINS_1/', "Unmapped_Reads")


## Compare CMV pre and post 
CMV_ALL <- read.delim('/well/immune-rep/shared/MISEQ/VIRAL_SEQ/VT_FIRST/BULK_CMV/FINAL_COUNTS_WIDE_VIRUS_All_Reads.txt', header=TRUE)
CMV_Unmapped <- read.delim('/well/immune-rep/shared/MISEQ/VIRAL_SEQ/HUMAN_FIRST/BULK_CMV_UNMAPPED/FINAL_COUNTS_WIDE_VIRUS_Unmapped_Reads.txt', header=TRUE)
sample_id <- str_split_fixed(rownames(CMV_ALL), "_", 3)
sample_id <- paste0("Donor", sample_id[,1], "_", sample_id[,2])
rownames(CMV_ALL) <- sample_id
CMV_ALL$Sample <- rownames(CMV_ALL)
sample_id <- str_split_fixed(rownames(CMV_Unmapped), "_", 3)
sample_id <- paste0("Donor", sample_id[,1], "_", sample_id[,2])
rownames(CMV_Unmapped) <- sample_id
CMV_Unmapped$Sample <- rownames(CMV_Unmapped)
CMV_ALL <- reshape2::melt(CMV_ALL, id.vars=c("Sample"))
colnames(CMV_ALL) <- c("Sample", "Virus", "Count_ALL")
CMV_Unmapped <- reshape2::melt(CMV_Unmapped, id.vars=c("Sample"))
colnames(CMV_Unmapped) <- c("Sample", "Virus", "Count_UNMAPPED")
CMV <- merge(CMV_Unmapped, CMV_ALL, by=c("Sample", "Virus"), all=TRUE)
CMV[is.na(CMV)] <- 0
CMV$TOTAL = (CMV$Count_UNMAPPED+ CMV$Count_ALL) 
## Detected Viruses
s <- aggregate(CMV$TOTAL, by=list(CMV$Virus), FUN=sum)
s <- s[s$x > 0,]
found_viruses <- as.character(s$Group.1)
CMV <- CMV[CMV$Virus %in% found_viruses,]
CMV <- merge(CMV, Virus_database, by.x="Virus", by.y="viral_genome")	
pdf(paste0("/well/immune-rep/shared/MISEQ/VIRAL_SEQ/CMV_comparison.pdf"), width=25, height=10)
p1 <- ggplot(CMV, aes(color=Virus_name, x=Count_ALL, y=Count_UNMAPPED)) + geom_point()  + theme_bw() + xlab("Counts Using All Reads + Stringent Human Mapping Paramters") + ylab("Counts Using Unmapped Reads + Less Stringent Mapping Paramters")+ggtitle(paste0("CMV Method Comparison")) +labs(colour="Virus") + geom_abline(intercept=0, slope=1, col="red") +facet_zoom(xlim=c(0,20), ylim=c(0,20)) +guides(colour=guide_legend(ncol=2))
plot(p1)
dev.off()

## Compare EBV pre and post
EBV_ALL <- read.delim('/well/immune-rep/shared/MISEQ/VIRAL_SEQ/VT_FIRST/BULK_EBV_LCL/FINAL_COUNTS_WIDE_VIRUS_All_Reads.txt', header=TRUE)
EBV_Unmapped <- read.delim('/well/immune-rep/shared/MISEQ/VIRAL_SEQ/HUMAN_FIRST/BULK_EBV_UNMAPPED/FINAL_COUNTS_WIDE_VIRUS_Unmapped_Reads.txt', header=TRUE)
EBV_ALL$Sample <- rownames(EBV_ALL)
EBV_Unmapped$Sample <- rownames(EBV_Unmapped)
EBV_ALL <- reshape2::melt(EBV_ALL, id.vars=c("Sample"))
colnames(EBV_ALL) <- c("Sample", "Virus", "Count_ALL")
EBV_Unmapped <- reshape2::melt(EBV_Unmapped, id.vars=c("Sample"))
colnames(EBV_Unmapped) <- c("Sample", "Virus", "Count_UNMAPPED")
EBV_Unmapped$Sample <- gsub(".Unmapped.out", "", EBV_Unmapped$Sample)
EBV <- merge(EBV_ALL, EBV_Unmapped, by=c("Sample", "Virus"), all = TRUE)
EBV[is.na(EBV)] <- 0
EBV$TOTAL = (EBV$Count_UNMAPPED+ EBV$Count_ALL) 
## Detected Viruses
s <- aggregate(EBV$TOTAL, by=list(EBV$Virus), FUN=sum)
s <- s[s$x > 0,]
found_viruses <- as.character(s$Group.1)
EBV <- EBV[EBV$Virus %in% found_viruses,]
EBV <- merge(EBV, Virus_database, by.x="Virus", by.y="viral_genome")	
pdf(paste0("/well/immune-rep/shared/MISEQ/VIRAL_SEQ/EBV_comparison.pdf"), width=25, height=10)
p1 <- ggplot(EBV, aes(color=Virus_name, x=Count_ALL, y=Count_UNMAPPED)) + geom_point()  + theme_bw() + xlab("Counts Using All Reads + Stringent Human Mapping Paramters") + ylab("Counts Using Unmapped Reads + Less Stringent Mapping Paramters")+ggtitle(paste0("CMV Method Comparison")) +labs(colour="Virus") + geom_abline(intercept=0, slope=1, col="red") +facet_zoom(xlim=c(0,20), ylim=c(0,20)) +guides(colour=guide_legend(ncol=2))
plot(p1)
dev.off()

## Compare HEPC pre and post
HEPC_ALL <- read.delim('/well/immune-rep/shared/MISEQ/VIRAL_SEQ/VT_FIRST/BULK_HEPC/FINAL_COUNTS_WIDE_VIRUS_All_Reads.txt', header=TRUE)
HEPC_Unmapped <- read.delim('/well/immune-rep/shared/MISEQ/VIRAL_SEQ/HUMAN_FIRST/BULK_HEPC_UNMAPPED/FINAL_COUNTS_WIDE_VIRUS_Unmapped_Reads.txt', header=TRUE)
HEPC_ALL$Sample <- rownames(HEPC_ALL)
HEPC_Unmapped$Sample <- rownames(HEPC_Unmapped)
HEPC_ALL <- reshape2::melt(HEPC_ALL, id.vars=c("Sample"))
colnames(HEPC_ALL) <- c("Sample", "Virus", "Count_ALL")
HEPC_Unmapped <- reshape2::melt(HEPC_Unmapped, id.vars=c("Sample"))
colnames(HEPC_Unmapped) <- c("Sample", "Virus", "Count_UNMAPPED")
HEPC_Unmapped$Sample <- gsub(".Unmapped.out", "", HEPC_Unmapped$Sample)
HEPC <- merge(HEPC_ALL, HEPC_Unmapped, by=c("Sample", "Virus"), all = TRUE)
HEPC[is.na(HEPC)] <- 0
HEPC$TOTAL = (HEPC$Count_UNMAPPED+ HEPC$Count_ALL) 
## Detected Viruses
s <- aggregate(HEPC$TOTAL, by=list(HEPC$Virus), FUN=sum)
s <- s[s$x > 0,]
found_viruses <- as.character(s$Group.1)
HEPC <- HEPC[HEPC$Virus %in% found_viruses,]
HEPC <- merge(HEPC, Virus_database, by.x="Virus", by.y="viral_genome")	
pdf(paste0("/well/immune-rep/shared/MISEQ/VIRAL_SEQ/HEPC_comparison.pdf"), width=30, height=15)
p1 <- ggplot(HEPC, aes(color=Virus_name, x=Count_ALL, y=Count_UNMAPPED)) + geom_point()  + theme_bw() + xlab("Counts Using All Reads + Stringent Human Mapping Paramters") + ylab("Counts Using Unmapped Reads + Less Stringent Mapping Paramters")+ggtitle(paste0("CMV Method Comparison")) +labs(colour="Virus") + geom_abline(intercept=0, slope=1, col="red") +facet_zoom(xlim=c(0,20), ylim=c(0,20)) +guides(colour=guide_legend(ncol=3))
plot(p1)
dev.off()


#-------------------------------------------------------------------------------------------------------------
# READ COUNTING ANALYSIS 
source('/gpfs2/well/immune-rep/shared/CODE/Viral-Seq/Analysis_Functions/CombineQC.R')
get_qc('/well/immune-rep/shared/MISEQ/VIRAL_SEQ/VT_FIRST/BULK_CMV/', "All_Reads", "CMV")
get_qc('/gpfs2/well/immune-rep/shared/MISEQ/VIRAL_SEQ/HUMAN_FIRST/BULK_CMV_UNMAPPED/', "Unmapped_Reads", "CMV")
get_qc('/well/immune-rep/shared/MISEQ/VIRAL_SEQ/VT_FIRST/BULK_EBV_LCL/', "All_Reads", "EBV")
get_qc('/gpfs2/well/immune-rep/shared/MISEQ/VIRAL_SEQ/HUMAN_FIRST/BULK_EBV_UNMAPPED/', "Unmapped_Reads", "EBV")
get_qc('/well/immune-rep/shared/MISEQ/VIRAL_SEQ/SEPSIS_GAINS_1/', "Unmapped_Reads")

EBV_QC_ALL <- read.delim('/well/immune-rep/shared/MISEQ/VIRAL_SEQ/VT_FIRST/BULK_EBV_LCL/QC_measures_All_Reads.txt', header=TRUE)
colnames(EBV_QC_ALL) <- paste0(colnames(EBV_QC_ALL), "_ALL")
colnames(EBV_QC_ALL)[colnames(EBV_QC_ALL)=="Sample_ALL"] <- "Sample"
colnames(EBV_QC_ALL)[colnames(EBV_QC_ALL) =="Virus_ALL"] <- "Virus"

EBV_QC_Unmapped <- read.delim('/well/immune-rep/shared/MISEQ/VIRAL_SEQ/HUMAN_FIRST/BULK_EBV_UNMAPPED/QC_measures_Unmapped_Reads.txt', header=TRUE)
colnames(EBV_QC_Unmapped) <- paste0(colnames(EBV_QC_Unmapped), "_Unmapped")
colnames(EBV_QC_Unmapped)[colnames(EBV_QC_Unmapped)=="Sample_Unmapped"] <- "Sample"
colnames(EBV_QC_Unmapped)[colnames(EBV_QC_Unmapped) =="Virus_Unmapped"] <- "Virus"
EBV_QC <- merge(EBV_QC_ALL, EBV_QC_Unmapped, by=c("Sample", "Virus"), all = TRUE)	
EBV_QC <- merge(EBV_QC, Virus_database, by.x="Virus", by.y="viral_genome")
EBV_QC[is.na(EBV_QC)]<-0
EBV_QC$VIRUS_USE_Unmapped[EBV_QC$VIRUS_USE_Unmapped==0] <- "NEGATIVE"

pdf(paste0("/well/immune-rep/shared/MISEQ/VIRAL_SEQ/EBV_Contig_comparison.pdf"), width=20, height=20)
p1 <- ggscatter(EBV_QC, y="Longest_contig_Unmapped", x="Longest_contig_ALL", add = "reg.line") +geom_point(aes(color=Virus_name), alpha=0.7)  + theme_bw() + xlab("Longest Contig Using All Reads + Stringent Human Mapping Paramters") + ylab("Longest Contig Using Unmapped Reads + Less Stringent Mapping Paramters")+ggtitle(paste0("EBV Method Comparison")) +labs(colour="Virus") + geom_abline(intercept=0, slope=1, col="red") +guides(colour=guide_legend(ncol=1)) +geom_vline(xintercept=100, col="blue")+geom_vline(xintercept=200, col="green")+geom_hline(yintercept=100, col="blue")+geom_hline(yintercept=200, col="green") + stat_cor(method = "pearson", label.x = 300, label.y =2000)+facet_zoom(xlim=c(0,200), ylim=c(0,200))
p2<- ggscatter(EBV_QC, y="Longest_contig_Unmapped", x="Longest_contig_ALL", add = "reg.line") +geom_point(aes(color=VIRUS_USE_Unmapped), alpha=0.7)  + theme_bw() + xlab("Longest Contig Using All Reads + Stringent Human Mapping Paramters") + ylab("Longest Contig Using Unmapped Reads + Less Stringent Mapping Paramters")+ggtitle(paste0("EBV Method Comparison")) +labs(colour="Virus") + geom_abline(intercept=0, slope=1, col="red") +guides(colour=guide_legend(ncol=1)) +geom_vline(xintercept=100, col="blue")+geom_vline(xintercept=200, col="green")+geom_hline(yintercept=100, col="blue")+geom_hline(yintercept=200, col="green") + stat_cor(method = "pearson", label.x = 300, label.y =2000)+facet_zoom(xlim=c(0,200), ylim=c(0,200))
plot(plot_grid(p1, p2, ncol=1))
p1 <- ggscatter(EBV_QC, y="Longest_contig_Unmapped", x="N_reads_Unmapped") +geom_point(aes(color=Virus_name), alpha=0.7)  + theme_bw() + xlab("Number of Reads (Unmapped Reads + Less Stringent Mapping Paramters)") + ylab("Longest Contig Using Unmapped Reads + Less Stringent Mapping Paramters")+ggtitle(paste0("EBV Method Comparison")) +labs(colour="Virus") +guides(colour=guide_legend(ncol=1)) +geom_vline(xintercept=5, col="blue")+geom_vline(xintercept=10, col="green")+geom_hline(yintercept=100, col="blue")+geom_hline(yintercept=200, col="green") +facet_zoom(xlim=c(0,75), ylim=c(0,300))
p2 <- ggscatter(EBV_QC, y="Longest_contig_Unmapped", x="N_reads_Unmapped") +geom_point(aes(color=VIRUS_USE_Unmapped), alpha=0.7)  + theme_bw() + xlab("Number of Reads (Unmapped Reads + Less Stringent Mapping Paramters)") + ylab("Longest Contig Using Unmapped Reads + Less Stringent Mapping Paramters")+ggtitle(paste0("EBV Method Comparison")) +labs(colour="Virus") +guides(colour=guide_legend(ncol=1)) +geom_vline(xintercept=5, col="blue")+geom_vline(xintercept=10, col="green")+geom_hline(yintercept=100, col="blue")+geom_hline(yintercept=200, col="green") +facet_zoom(xlim=c(0,75), ylim=c(0,300))
plot(plot_grid(p1, p2, ncol=1))
p1 <- ggscatter(EBV_QC, y="Longest_contig_Unmapped", x="Coverage_per_read_Ratio_Unmapped") +geom_point(aes(color=Virus_name), alpha=0.7)  + theme_bw() + xlab("Coverage per read (Unmapped Reads + Less Stringent Mapping Paramters)") + ylab("Longest Contig (Unmapped Reads + Less Stringent Mapping Paramters)")+ggtitle(paste0("EBV Method Comparison")) +labs(colour="Virus") +guides(colour=guide_legend(ncol=1)) +geom_vline(xintercept=5, col="blue")+geom_vline(xintercept=10, col="green")+geom_hline(yintercept=100, col="blue")+geom_hline(yintercept=200, col="green") +facet_zoom(xlim=c(0,1.5), ylim=c(0,300))
p2 <- ggscatter(EBV_QC, y="Longest_contig_Unmapped", x="Coverage_per_read_Ratio_Unmapped") +geom_point(aes(color=VIRUS_USE_Unmapped), alpha=0.7)  + theme_bw() + xlab("Coverage per read (Unmapped Reads + Less Stringent Mapping Paramters)") + ylab("Longest Contig (Unmapped Reads + Less Stringent Mapping Paramters)")+ggtitle(paste0("EBV Method Comparison")) +labs(colour="Virus") +guides(colour=guide_legend(ncol=1)) +geom_vline(xintercept=5, col="blue")+geom_vline(xintercept=10, col="green")+geom_hline(yintercept=100, col="blue")+geom_hline(yintercept=200, col="green") +facet_zoom(xlim=c(0,1.5), ylim=c(0,300))
plot(plot_grid(p1, p2, ncol=1))
p1 <- ggscatter(EBV_QC, y="Longest_contig_Unmapped", x="Coverage_Percent_1000bp_Unmapped") +geom_point(aes(color=Virus_name), alpha=0.7)  + theme_bw() + xlab("Coverage 1000bp Percent (Unmapped Reads + Less Stringent Mapping Paramters)") + ylab("Longest Contig (Unmapped Reads + Less Stringent Mapping Paramters)")+ggtitle(paste0("EBV Method Comparison")) +labs(colour="Virus") +guides(colour=guide_legend(ncol=1)) +geom_vline(xintercept=40, col="blue")+geom_vline(xintercept=50, col="green")+geom_hline(yintercept=100, col="blue")+geom_hline(yintercept=200, col="green") +facet_zoom(xlim=c(0,25), ylim=c(0,300))
p2 <- ggscatter(EBV_QC, y="Longest_contig_Unmapped", x="Coverage_Percent_1000bp_Unmapped") +geom_point(aes(color=VIRUS_USE_Unmapped), alpha=0.7)  + theme_bw() + xlab("Coverage 1000bp Percent (Unmapped Reads + Less Stringent Mapping Paramters)") + ylab("Longest Contig (Unmapped Reads + Less Stringent Mapping Paramters)")+ggtitle(paste0("EBV Method Comparison")) +labs(colour="Virus") +guides(colour=guide_legend(ncol=1)) +geom_vline(xintercept=40, col="blue")+geom_vline(xintercept=50, col="green")+geom_hline(yintercept=100, col="blue")+geom_hline(yintercept=200, col="green") +facet_zoom(xlim=c(0,25), ylim=c(0,300))
plot(plot_grid(p1, p2, ncol=1))
dev.off()


CMV_QC_ALL <- read.delim('/well/immune-rep/shared/MISEQ/VIRAL_SEQ/VT_FIRST/BULK_CMV/QC_measures_All_Reads.txt', header=TRUE)
colnames(CMV_QC_ALL) <- paste0(colnames(CMV_QC_ALL), "_ALL")
colnames(CMV_QC_ALL)[colnames(CMV_QC_ALL)=="Sample_ALL"] <- "Sample"
colnames(CMV_QC_ALL)[colnames(CMV_QC_ALL) =="Virus_ALL"] <- "Virus"

CMV_QC_Unmapped <- read.delim('/well/immune-rep/shared/MISEQ/VIRAL_SEQ/HUMAN_FIRST/BULK_CMV_UNMAPPED/QC_measures_Unmapped_Reads.txt', header=TRUE)
colnames(CMV_QC_Unmapped) <- paste0(colnames(CMV_QC_Unmapped), "_Unmapped")
colnames(CMV_QC_Unmapped)[colnames(CMV_QC_Unmapped)=="Sample_Unmapped"] <- "Sample"
colnames(CMV_QC_Unmapped)[colnames(CMV_QC_Unmapped) =="Virus_Unmapped"] <- "Virus"
CMV_QC <- merge(CMV_QC_ALL, CMV_QC_Unmapped, by=c("Sample", "Virus"), all = TRUE)	
CMV_QC <- merge(CMV_QC, Virus_database, by.x="Virus", by.y="viral_genome")
CMV_QC[is.na(CMV_QC)]<-0

pdf(paste0("/well/immune-rep/shared/MISEQ/VIRAL_SEQ/CMV_Contig_comparison.pdf"), width=20, height=10)
p1 <- ggscatter(CMV_QC, y="Longest_contig_Unmapped", x="Longest_contig_ALL", add = "reg.line") +geom_point(aes(color=Virus_name), alpha=0.7)  + theme_bw() + xlab("Longest Contig Using All Reads + Stringent Human Mapping Paramters") + ylab("Longest Contig Using Unmapped Reads + Less Stringent Mapping Paramters")+ggtitle(paste0("EBV Method Comparison")) +labs(colour="Virus") + geom_abline(intercept=0, slope=1, col="red") +guides(colour=guide_legend(ncol=1)) +geom_vline(xintercept=100, col="blue")+geom_vline(xintercept=200, col="green")+geom_hline(yintercept=100, col="blue")+geom_hline(yintercept=200, col="green") + stat_cor(method = "pearson", label.x = 300, label.y =300)+facet_zoom(xlim=c(0,200), ylim=c(0,200))
plot(p1)
p1 <- ggscatter(CMV_QC, y="Longest_contig_Unmapped", x="N_reads_Unmapped") +geom_point(aes(color=Virus_name), alpha=0.7)  + theme_bw() + xlab("Number of Reads (Unmapped Reads + Less Stringent Mapping Paramters)") + ylab("Longest Contig Using Unmapped Reads + Less Stringent Mapping Paramters")+ggtitle(paste0("EBV Method Comparison")) +labs(colour="Virus") +guides(colour=guide_legend(ncol=1)) +geom_vline(xintercept=5, col="blue")+geom_vline(xintercept=10, col="green")+geom_hline(yintercept=100, col="blue")+geom_hline(yintercept=200, col="green") +facet_zoom(xlim=c(0,75), ylim=c(0,300))
plot(p1)
p1 <- ggscatter(CMV_QC, y="Longest_contig_Unmapped", x="Coverage_per_read_Ratio_Unmapped") +geom_point(aes(color=Virus_name), alpha=0.7)  + theme_bw() + xlab("Coverage per read (Unmapped Reads + Less Stringent Mapping Paramters)") + ylab("Longest Contig (Unmapped Reads + Less Stringent Mapping Paramters)")+ggtitle(paste0("EBV Method Comparison")) +labs(colour="Virus") +guides(colour=guide_legend(ncol=1)) +geom_vline(xintercept=5, col="blue")+geom_vline(xintercept=10, col="green")+geom_hline(yintercept=100, col="blue")+geom_hline(yintercept=200, col="green") +facet_zoom(xlim=c(0,1.5), ylim=c(0,300))
plot(p1)
p1 <- ggscatter(CMV_QC, y="Longest_contig_Unmapped", x="Coverage_Percent_1000bp_Unmapped") +geom_point(aes(color=Virus_name), alpha=0.7)  + theme_bw() + xlab("Coverage 1000bp Percent (Unmapped Reads + Less Stringent Mapping Paramters)") + ylab("Longest Contig (Unmapped Reads + Less Stringent Mapping Paramters)")+ggtitle(paste0("EBV Method Comparison")) +labs(colour="Virus") +guides(colour=guide_legend(ncol=1)) +geom_vline(xintercept=40, col="blue")+geom_vline(xintercept=50, col="green")+geom_hline(yintercept=100, col="blue")+geom_hline(yintercept=200, col="green") +facet_zoom(xlim=c(0,25), ylim=c(0,300))
plot(p1)
dev.off()



#+facet_zoom(xlim=c(0,200), ylim=c(0,200))
#-------------------------------------------------------------------------------------------------------------
# Reading in samples and making matrix for QC_Results 
# Reading in samples and making matrix for unique reads post feature counts:
paths <- list.files(paste0(outputdir, "/ViRNA_SEQ_UNIQUE_Mapping_Analysis/"), recursive=TRUE, full.name=TRUE)
paths <- grep("QC_un", paths, value=TRUE)
df_list <- lapply(paths, fread, sep="\t")
# identify the sample they came from 
d <- str_split(paths, "Reports/") 
d <- sapply(d, "[[", 2)  
d <- gsub(".Unmapped.outQC_unfiltered.txt", "", d)
# Name each dataframe with the run and filename
names(df_list) <- d
# Create combined dataframe  
## We need to remove samples where no viruses were detected 
df_list <- Filter(function(x) dim(x)[1] > 0, df_list)
df <- rbindlist(df_list, idcol=TRUE)
viral_qc <- df
### THIS iS DATA FRAME OF QC (post FEATURECOUNTS)
write.table(viral_qc, paste0(outputdir, "FINAL_QC_ALL.txt"), sep="\t", row.names=FALSE)
#-------------------------------------------------------------------------------------------------------------
## Summary Coverage for each virus 
paths <- list.files(paste0(outputdir, "/ViRNA_SEQ_UNIQUE_Mapping_Analysis"), recursive=TRUE, full.name=TRUE)
paths <- grep("_genome_coverage.txt", paths, value=TRUE)
path_df <- data.frame(paths)
d <- str_split(paths, "Viral_BAM_files/") 
d <- sapply(d, "[[", 2)  
d <- gsub("_genome_coverage.txt", "", d)
path_df$Virus <- d
d <- str_split(paths, "ViRNA_SEQ_UNIQUE_Mapping_Analysis/") 
d <- sapply(d, "[[", 2)  
d <- str_split(d, "/Viral_BAM_files") 
d <- sapply(d, "[[", 1) 
d <- gsub(".Unmapped.out", "", d)
path_df$Sample <- d
coverage_list <- c()
for(i in 1:length(unique(path_df$Virus))){
	paths_read <- path_df$paths[path_df$Virus==unique(path_df$Virus)[i]]
	paths_sample <- path_df$Sample[path_df$Virus==unique(path_df$Virus)[i]]
	df_list <- lapply(paths_read, fread, sep="\t")
	df <- data.frame(do.call(cbind, df_list))
	df <- df[, !colnames(df) %like% "V1"]
	write.table(df, paste0(outputdir, "Coverage_all_samples_", unique(path_df$Virus)[i], ".txt"), sep="\t", row.names=FALSE)
	df <- list(df)
	names(df) <- unique(path_df$Virus)[i]
	coverage_list <- c(coverage_list, df)
	}
#-------------------------------------------------------------------------------------------------------------
## Plot QC results
pdf(paste0(outputdir, "QC_summary.pdf"), width=15, height=15)
for (i in 1:length(unique(viral_qc$V1))){
	QC_result <- viral_qc[viral_qc$V1==unique(viral_qc$V1)[i],]
	Mean_mapping_length <- 101
	minread <- 50
	name_virus <- Virus_database$Virus_name[Virus_database$viral_genome==unique(viral_qc$V1)[i]]
	p1 <- ggplot(QC_result, aes(color=.id, x=N_reads, y=(N_unique_reads))) + geom_point()  + theme_classic() + xlab("Number of Mapped Reads") + ylab("Number Uniquely Mapped Reads")  + ggtitle(name_virus)+ expand_limits(x = 0, y = 0)+ guides(color="none")
	p2 <- ggplot(QC_result, aes(color=.id, x=N_unique_reads, y=(Spatial_distribution*100))) + geom_point()  + theme_classic() + xlab("Number of Uniquely Mapped Reads") + ylab("% Mapped genome") + ylim(0, 100) + ggtitle(name_virus)+ guides(color="none")
	p3 <- ggplot(QC_result, aes(color=.id, x=Sequence_entropy, y=(Spatial_distribution*100))) + geom_point()  + theme_classic() + xlab("Sequence Entropy") + ylab("% Mapped genome")  + ylim(0, 100) + ggtitle(name_virus)+ guides(color="none")
	p4 <- ggplot(QC_result, aes(color=.id, x=Longest_contig, y=DUST_score)) + geom_point()  + theme_classic() + xlab("Longest Contig (nt)") + ylab("DUST Score") + ggtitle(name_virus)+ guides(color="none")
	p5 <- ggplot(QC_result, aes(color=.id, x=Mean_read_quality, y=Sd_read_quality)) + geom_point()  + theme_classic() + xlab("Mean Read Quality") + ylab("SD Read Quality") + ggtitle(name_virus)+ guides(color="none")
	p6 <- ggplot(QC_result, aes(color=.id, x=Longest_contig, y=Coverage_per_read_Ratio)) + geom_point()  + theme_classic() + xlab("Longest Contig") + ylab("Coverage per Read Ratio") + ggtitle(name_virus) + guides(color="none")
	p7 <- ggplot(QC_result, aes(color=.id, x=(Spatial_distribution*100), y=Coverage_per_read_Ratio)) + geom_point() + theme_classic() + xlab("% Mapped genome") + ylab("Coverage per Read Ratio") + ggtitle(name_virus) + xlim(0, 100) + guides(color="none")
	p8 <- ggplot(QC_result, aes(color=.id, y=Coverage_Percent_1000bp, x=N_unique_reads)) + geom_point()  + theme_classic() + ylab("% Coverage for 1000bp with highest coverage") + xlab("Number Uniquely Mapped Reads") + ggtitle(name_virus)+ guides(color="none")
	if(length(rownames(QC_result))>=1){
		QC_result <- QC_result[QC_result$N_unique_reads>=1, ]
		if(length(rownames(QC_result))>1){
			QC_scale <- data.frame(QC_result)
			rownames(QC_scale) <- QC_scale$.id
			QC_scale <- QC_scale[, c(3:17)]
			cols_discard <- apply(QC_scale, 2, function(x) length(unique(x)))
			cols_discard <- names(cols_discard[cols_discard<2])
			QC_scale <- QC_scale[, !colnames(QC_scale) %in% cols_discard]
			QC_scaled <- scale(QC_scale)
			QC_pca <- prcomp(QC_scaled, center = TRUE)
			eigencoordinates <- data.frame(QC_pca$x)
			p9 <- ggplot(eigencoordinates, aes( x=PC1, y=PC2)) + geom_point()  + theme_classic() + ylab("PC2") + xlab("PC1") + ggtitle(name_virus)+ guides(color="none")
			plot(plot_grid( p1, p2, p3, p4, p5, p6, p7, p8, p9,  ncol=3, labels = "AUTO"))
			} else {
		plot(plot_grid( p1, p2, p3, p4, p5, p6, p7, p8,  ncol=3, labels = "AUTO"))
		}
	}
}
dev.off()
#-------------------------------------------------------------------------------------------------------------









## Ready to explore!!!
viruses <- grep("refseq|NC", colnames(final_counts), value=TRUE)
cases_less500 <- c()
cases_less5000 <- c()
cases_more5000 <- c()

for(i in 1:length(viruses)){
	column_name <- viruses[i]
	column <- as.numeric(final_counts[,column_name])
	p <- sum(column)
	if(p >0 & p<500){
	cases_less500 <- c(cases_less500, column_name)
	}
	if (p>=500 & p<5000){
	cases_less5000 <- c(cases_less5000, column_name)
	}
	if(p>=5000){
	cases_more5000 <- c(cases_more5000, column_name)
	}
	}

cases <-c(cases_less500, cases_less5000, cases_more5000)

## Plotting Results of Unique analysis 
cases_1 <- cases
# subset the dataframe
final_counts_postive <- final_counts[, cases_1]
final_counts_postive$Sample <- rownames(final_counts_postive)
final_counts_postive[cases_1] <- sapply(final_counts_postive[cases_1],as.numeric)
## reshape long format
data_long <- gather(final_counts_postive, key = viral_genome, value = count, -Sample)
data_long$logcount <- log10(data_long$count)
data_long <- merge(data_long, Virus_database, by="viral_genome")
data_long$Mapping <- "Multi"
multi_mapping_counts <- data_long
## Now we can do a trial plot!
pdf("Plots/Viral_hits_sepsis_multi.pdf", width=25, height=10)
p <- ggplot(data_long, aes(x=Virus_name, y=count)) + geom_boxplot() + theme_bw()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.title = element_text(color = "blue", size = 100)) + stat_summary(fun=mean, geom="point", shape=23, size=2, col="red")+ facet_zoom(ylim = c(0, 550)) +xlab("Virus") + ylab("Total Counts") + ggtitle("ViRNAseq Multi Mapping Analysis")
plot(p)
dev.off()

## ALL POSSIBLE VIRUSES 
all_identified_viruses <- cases 
write.table(multi_mapping_counts, "Counts/Multi_mapping_counts.txt", sep="\t")
write.table(all_identified_viruses , "Counts/Multi_mapping_identified_viruses.txt", sep="\t")

##---------------------------------------------------------------------------------------------
## Repeat for uniquely mapping viruses. 

# Reading in samples and making matrix 
# FOR UNIQUE READS 
paths <- list.dirs('/well/immune-rep/users/kvi236/VIRUS/Sepsis_Bulk/ViRNA_SEQ_UNIQUE_Mapping_Analysis')
paths <- grep("Reports", paths, value=TRUE)

## Sequentially read in and merge data
for(i in 1:length(paths)){
      path <- paths[i] 
	  counts_file <- paste0(path, "/")
	  files <- list.files(counts_file)
	  files <- grep("COUNTS_SUMMARY", files, value =TRUE)
	  name <- unlist(str_split(files, "[_]"))
	  name <- name[[4]]
	  name <- unlist(str_split(name, "[.]"))
	  name <- name[[1]]
	  x <- read.delim((paste0(counts_file, files)), header=TRUE)
	  x <- x[, c("Chr", "count", "VirusFilter")]
	  colnames(x) <- c("Genome", name, paste0(name, "_Filtering"))
	  ## Set up dataframe 
	  if (i ==1){
		Viral_counts <- x
	  } else {	  
	  Viral_counts <- merge(Viral_counts, x, by="Genome")
	  }
	}
	
counts <- grep("Filtering", colnames(Viral_counts), value=TRUE)
samples<- colnames(Viral_counts)[colnames(Viral_counts) %notin% counts]
counts_only <- Viral_counts[samples]

## Rearrange to have better genome names 
human_chromosomes <- paste(0:22)
human_chromosomes <- c(human_chromosomes, "X", "Y")

for (i in 1:length(counts_only$Genome)){
	genome <- counts_only$Genome[i] 
	if (genome %in% human_chromosomes){
		counts_only$Genome[i] <- paste0("Chr", genome)
	}
}

# dim(counts_only) should be number of samples + 1 !)
s <- as.data.frame(t(counts_only))
colnames(s) <- s[1,]
s <- s[-1,]
final_counts <- s
final_counts$sample <- rownames(final_counts)

## Ready to explore!!!
viruses <- grep("refseq|NC", colnames(final_counts), value=TRUE)
cases_less500 <- c()
cases_less5000 <- c()
cases_more5000 <- c()

for(i in 1:length(viruses)){
	column_name <- viruses[i]
	column <- as.numeric(final_counts[,column_name])
	p <- sum(column)
	if(p >0 & p<500){
	cases_less500 <- c(cases_less500, column_name)
	}
	if (p>=500 & p<5000){
	cases_less5000 <- c(cases_less5000, column_name)
	}
	if(p>=5000){
	cases_more5000 <- c(cases_more5000, column_name)
	}
	}

cases <-c(cases_less500, cases_less5000, cases_more5000)
## ALL POSSIBLE VIRUSES 
all_identified_unique_viruses <- cases 

## Check to see whether all_viruses encapsulates those with just no unique reads
all_identified_viruses[all_identified_viruses %notin% all_identified_unique_viruses]
## Below should be 0 
all_identified_unique_viruses[all_identified_unique_viruses %notin% all_identified_viruses]
## Yes!

## Plotting Results of Unique analysis using all the viruses (e.g. all_identified_viruses NOT just unique viruses.  
cases_1 <- all_identified_viruses
# subset the dataframe
final_counts_postive <- final_counts[, cases_1]
final_counts_postive$Sample <- rownames(final_counts_postive)
final_counts_postive[cases_1] <- sapply(final_counts_postive[cases_1],as.numeric)
## reshape long format
data_long <- gather(final_counts_postive, key = viral_genome, value = count, -Sample)
data_long$logcount <- log10(data_long$count)
data_long <- merge(data_long, Virus_database, by="viral_genome")
data_long$Mapping <- "Unique"
unique_mapping_counts <- data_long
## Now we can do a trial plot!
pdf("Plots/Viral_hits_sepsis_unique.pdf", width=25, height=10)
p <- ggplot(data_long, aes(x=Virus_name, y=count)) + geom_boxplot() + theme_bw()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.title = element_text(color = "blue", size = 100)) + stat_summary(fun=mean, geom="point", shape=23, size=2, col="red")+ facet_zoom(ylim = c(0, 550)) +xlab("Virus") + ylab("Total Counts") + ggtitle("ViRNAseq Unique Mapping Analysis")
plot(p)
dev.off()

## ALL POSSIBLE VIRUSES 
all_identified_unique_viruses <- cases 
write.table(unique_mapping_counts, "Counts/Unique_mapping_counts.txt", sep="\t")
write.table(all_identified_unique_viruses , "Counts/Unique_mapping_identified_viruses.txt", sep="\t")

## MERGING DATA 
## CHeck that we have the right subset
## Both should be 0! 
unique_mapping_counts[unique_mapping_counts$viral_genome %notin% multi_mapping_counts$viral_genome,]
multi_mapping_counts[multi_mapping_counts$viral_genome %notin% unique_mapping_counts$viral_genome,]

## Now we combine the two datasets which both include the same number of viruses
## Combine the data 
all_counts <- rbind(unique_mapping_counts, multi_mapping_counts)
pdf("Plots/Viral_hits_sepsis_all.pdf", width=25, height=10)
p <- ggplot(all_counts, aes(x=Virus_name, y=count)) + geom_boxplot() + theme_bw()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.title = element_text(color = "blue", size = 100)) + stat_summary(fun=mean, geom="point", shape=23, size=2, col="red")+  facet_wrap(~ Mapping) +xlab("Virus") + ylab("Total Counts") + ggtitle("ViRNAseq Multi Mapping Analysis")
plot(p)
dev.off()

# Save Counts
write.table(all_counts, "Counts/All_counts.txt", sep="\t")


## Look at any viruses with lots of multi-mapping reads  e.g. difference between unique and multimapping is high
df = data.frame(sample_id = character(), viral_genome <- character(), count_unique <- numeric(), count_multi <- numeric(), Difference <- numeric())
colnames(df) <- c("sample_id", "viral_genome", "Count_Unique", "Count_Multi")
unique_mapping_counts$count <- as.numeric(unique_mapping_counts$count)
multi_mapping_counts$count <- as.numeric(multi_mapping_counts$count)
for(i in 1:length(unique_mapping_counts[,1])){
		sample_id <- unique_mapping_counts$Sample[i]
		sample_count <- unique_mapping_counts$count[i]
		genome_unique <- unique_mapping_counts$Virus_name[i]
		sample_count_multi <- multi_mapping_counts$count[multi_mapping_counts$Virus_name ==genome_unique & multi_mapping_counts$Sample == sample_id]
		row <- c(sample_id, genome_unique, sample_count, sample_count_multi)
        df <- rbind(df, row)
		}
	
colnames(df) <- c("sample_id", "viral_genome", "Count_Unique", "Count_Multi")
df$Count_Unique <- as.numeric(df$Count_Unique)
df$Count_Multi <- as.numeric(df$Count_Multi)
df$Difference <- df$Count_Multi - df$Count_Unique
df$Percentage_Unique <- df$Count_Unique / df$Count_Multi * 100
df$Percentage_Multi <- df$Difference / df$Count_Multi * 100

#Save df
write.table(df, "Counts/Counts_DF.txt", sep="\t")


## Plot the differences to identify any odd viruses 
pdf("Plots/Viral_hits_sepsis_difference.pdf", width=25, height=10)
p <- ggplot(df, aes(x=viral_genome, y=Percentage_Unique)) + geom_boxplot() + theme_bw()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.title = element_text(color = "blue", size = 100)) +xlab("Virus") + ylab("Percentage of Primary Mapped Reads which are Unique") + ggtitle("ViRNAseq Unique Mapping Analysis")
plot(p)
p <- ggplot(df, aes(x=viral_genome, y=Percentage_Multi)) + geom_boxplot() + theme_bw()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.title = element_text(color = "blue", size = 100)) +xlab("Virus") + ylab("Percentage of Primary Mapped Reads which are Non-Unique") + ggtitle("ViRNAseq Unique Mapping Analysis")
plot(p)
p <- ggplot(df, aes(x=viral_genome, y=Difference)) + geom_boxplot() + theme_bw()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.title = element_text(color = "blue", size = 100)) + facet_zoom(ylim = c(0, 550)) +xlab("Virus") + ylab("Difference in Counts (Multi vs Unique Mapping)") + ggtitle("ViRNAseq Unique Mapping Analysis")
plot(p)
dev.off()

#------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------
## READING IN DATA SO DON'T HAVE TO REDO!!
multi_mapping_counts <- read.delim("Counts/Unique_mapping_counts.txt", sep="\t")
all_identified_unique_viruses <- read.delim("Counts/Unique_mapping_identified_viruses.txt", sep="\t")
unique_mapping_counts <- read.delim("Counts/Multi_mapping_counts.txt", sep="\t")
all_identified_viruses <- read.delim("Counts/Multi_mapping_identified_viruses.txt", sep="\t")
all_counts <- read.delim("Counts/All_counts.txt", sep="\t")
df <- read.delim("Counts/Counts_DF.txt", sep="\t")

#------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------

## Look at positives in dataset 
cases_df <- data.frame(viral_name = character(), number_cases = numeric(), percentage = numeric())
for(i in 1:length(all_identified_viruses[,1])){
	virus <- all_identified_viruses[i,]
	virus_name <- Virus_database$Virus_name[Virus_database$viral_genome==virus]
	subset2 <- df[df$viral_genome==virus_name,]
	count <- c()
    count <- length(subset2$Count_Unique[subset2$Count_Unique >=1])
	percent_postitive <- count / 903 * 100
	row <- c(virus_name, count, percent_postitive)
	cases_df <- rbind(row, cases_df)
}
colnames(cases_df) <- c("viral_name", "cases_1", "percent_1")

cases_df_10 <- data.frame(viral_name = character(), number_cases = numeric(), percentage = numeric())
for(i in 1:length(all_identified_viruses[,1])){
	virus <- all_identified_viruses[i,]
	virus_name <- Virus_database$Virus_name[Virus_database$viral_genome==virus]
	subset2 <- df[df$viral_genome==virus_name,]
	count <- c()
    count <- length(subset2$Count_Unique[subset2$Count_Unique >=10])
	percent_postitive <- count / 903 * 100
	row <- c(virus_name, count, percent_postitive)
	cases_df_10 <- rbind(row, cases_df_10)
}
colnames(cases_df_10) <- c("viral_name", "cases_10", "percent_10")

cases_df_30<- data.frame(viral_name = character(), number_cases = numeric(), percentage = numeric())
for(i in 1:length(all_identified_viruses[,1])){
	virus <- all_identified_viruses[i,]
	virus_name <- Virus_database$Virus_name[Virus_database$viral_genome==virus]
	subset2 <- df[df$viral_genome==virus_name,]
	count <- c()
    count <- length(subset2$Count_Unique[subset2$Count_Unique >=30])
	percent_postitive <- count / 903 * 100
	row <- c(virus_name, count, percent_postitive)
	cases_df_30 <- rbind(row, cases_df_30)
}
colnames(cases_df_30) <- c("viral_name", "cases_30", "percent_30")

cases_df_50<- data.frame(viral_name = character(), number_cases = numeric(), percentage = numeric())
for(i in 1:length(all_identified_viruses[,1])){
	virus <- all_identified_viruses[i,]
	virus_name <- Virus_database$Virus_name[Virus_database$viral_genome==virus]
	subset2 <- df[df$viral_genome==virus_name,]
	count <- c()
    count <- length(subset2$Count_Unique[subset2$Count_Unique >=50])
	percent_postitive <- count / 903 * 100
	row <- c(virus_name, count, percent_postitive)
	cases_df_50 <- rbind(row, cases_df_50)
}
colnames(cases_df_50) <- c("viral_name", "cases_50", "percent_50")


complete <- merge(cases_df_10, cases_df, by="viral_name")
complete <- merge(complete, cases_df_30, by="viral_name")
complete <- merge(complete, cases_df_50, by="viral_name")


complete_percentage <- complete[ ,c("viral_name", "percent_1", "percent_10", "percent_30", "percent_50")]
new_format <- gather(complete_percentage, key = min_reads, value = percentage, -viral_name)
new_format$percentage <- as.numeric(new_format$percentage)
## Plot the differences to identify any odd viruses 
pdf("Plots/Viral_incidence.pdf", width=30, height=10)
p <- ggplot(new_format, aes(x=viral_name, y=percentage, fill=min_reads)) + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.15) + theme_bw()+ facet_zoom(ylim = c(0, 15)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.title = element_text(color = "blue")) +xlab("Virus") + ylab("Percentage Positive") + ggtitle("ViRNAseq Unique Mapping Analysis") + labs(fill="Min reads for positive id") + scale_fill_discrete(labels = c("1", "10", "30", "50"))
plot(p)
dev.off()

#------------------------------------------------
##Save and Read 
write.table(complete, 'Counts/viral_incidence.txt', sep='\t')
write.table(new_format, 'Counts/viral_incidence_new_format.txt', sep='\t')
complete<- read.delim('Counts/viral_incidence.txt', sep='\t')
new_format <- read.delim('Counts/viral_incidence_new_format.txt', sep='\t')
#----------------------------------------------------

## ASSIGN SRS membership 
srs <- read.delim('/gpfs2/well/immune-rep/users/kvi236/VIRUS/Sepsis_Scripts/Sample_Sheets/SRS_model_all_samples_864.txt', sep="\t")
key <- read.delim('/gpfs2/well/immune-rep/users/kvi236/VIRUS/Sepsis_Scripts/Sample_Sheets/Sample_key.txt', sep="\t")

## Add srs information
colnames(df)[colnames(df)=="sample_id"]  <- "SangerSampleID"
df <- merge(df, key, by="SangerSampleID", all.x=TRUE)
df_2 <- df
df_2$SampleID <- df_2$GAinSID
df_2 <- merge(df_2, srs, by="SampleID", , all.x=TRUE)


## Plot the differences to identify any odd viruses 
pdf("Plots/Viral_cases_SRS.pdf", width=20, height=10)
p <- ggplot(df_2, aes(x=viral_genome, y=Count_Unique))+ geom_boxplot() + theme_bw()+ facet_wrap(~SRSModel) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.title = element_text(color = "blue")) +xlab("Virus") + ylab("Percentage Positive") + ggtitle("ViRNAseq Unique Mapping Analysis")
plot(p)
dev.off()

##-------------------------------------------------
write.table(df_2, 'Counts/df_srs.txt', sep='\t')
df <- read.delim('Counts/df_srs.txt', sep='\t')
##-------------------------------------------------

## Investigating Genome coverage on all 903 samples and on an individual basis'

##Run parralised function on single and multi-mapping 
source('/gpfs2/well/immune-rep/users/kvi236/VIRUS/Sepsis_Scripts/auxillary_functions.R')
directory <- '/well/immune-rep/users/kvi236/VIRUS/Sepsis_Bulk/ViRNA_SEQ_UNIQUE_Mapping_Analysis'
for(i in 1:length(all_identified_viruses[,1])){
    virus_use <- all_identified_viruses$x[i]
	map_viral_coverage_unique(virus_use, directory)
	}
on.exit(stopCluster(c1))

directory <- '/well/immune-rep/users/kvi236/VIRUS/Sepsis_Bulk/ViRNA_SEQ_Multi_Mapping_Analysis'
for(i in 1:length(all_identified_viruses[,1])){
    virus_use <- all_identified_viruses$x[i]
	map_viral_coverage_multi(virus_use, directory)
	}
on.exit(stopCluster(c1))


##---------------------------------------------------------------------
##---------------------------------------------------------------------
##---------------------------------------------------------------------
## Investigating Coverage 
## Per individual 
paths <- list.files('/gpfs2/well/immune-rep/users/kvi236/VIRUS/Sepsis_Scripts/Plots/Coverage_unique', full.names=TRUE)
paths <- grep("*individual_coverage.txt", paths, value=TRUE)

trial <- data.frame()
for(i in 1:length(paths)){
	path <- paths[i]
	coverage_individual <- read.delim(path, sep='\t')
	a <- unlist(str_split(path, "/Coverage_unique/"))[2]
	a <- unlist(str_split(a, "_unique_individual_coverage.txt"))[1]
	virus_id <- a
	Virus_eng <- Virus_database$Virus_name[Virus_database$viral_genome == virus_id]
    coverage_individual$viral_genome <- Virus_eng
	coverage_individual$Virus_id <- virus_id
	new <- merge(coverage_individual, df, by=c("SangerSampleID", "viral_genome"))
	trial <- rbind(new, trial)
}
results <- trial[, c("SangerSampleID", "GAinSID", "viral_genome", "Virus_id", "SRSModel.x", "Count_Unique", "Count_Multi", "Difference", "Percentage_Unique", "Percentage_Multi", "Coverage_bp", "Percent_covered")] 
colnames(results) <- c("SangerSampleID", "GAinSID", "viral_genome", "Virus_id", "SRSModel", "Count_Unique", "Count_Multi", "Difference", "Percentage_Unique", "Percentage_Multi", "Coverage_bp_unique", "Percent_covered_unique")
results$SRSModel <- as.factor(results$SRSModel)

paths <- list.files('/gpfs2/well/immune-rep/users/kvi236/VIRUS/Sepsis_Scripts/Plots/Coverage_multi', full.names=TRUE)
paths <- grep("*individual_coverage.txt", paths, value=TRUE)

trial_multi <- data.frame()
for(i in 1:length(paths)){
	path <- paths[i]
	coverage_individual <- read.delim(path, sep='\t')
	a <- unlist(str_split(path, "/Coverage_multi/"))[2]
	a <- unlist(str_split(a, "_multi_individual_coverage.txt"))[1]
	virus_id <- a
	Virus_eng <- Virus_database$Virus_name[Virus_database$viral_genome == virus_id]
    coverage_individual$viral_genome <- Virus_eng
	coverage_individual$Virus_id <- virus_id
	new <- merge(coverage_individual, df, by=c("SangerSampleID", "viral_genome"))
	trial_multi <- rbind(new, trial_multi)
}
results_multi <- trial_multi[, c("SangerSampleID", "GAinSID", "viral_genome", "Virus_id", "SRSModel.x", "Count_Unique", "Count_Multi", "Difference", "Percentage_Unique", "Percentage_Multi", "Coverage_bp", "Percent_covered")] 
colnames(results_multi) <- c("SangerSampleID", "GAinSID", "viral_genome", "Virus_id", "SRSModel", "Count_Unique", "Count_Multi", "Difference", "Percentage_Unique", "Percentage_Multi", "Coverage_bp_multi", "Percent_covered_multi")
results_multi$SRSModel <- as.factor(results_multi$SRSModel)


new <- merge(results_unique, results_multi, by=c("SangerSampleID", "viral_genome"))
new <- new[, c("SangerSampleID", "GAinSID.x", "viral_genome", "Virus_id.x", "SRSModel.x", "Count_Unique.x", "Count_Multi.x", "Difference.x", "Percentage_Unique.x", "Percentage_Multi.x", "Coverage_bp_unique", "Percent_covered_unique","Coverage_bp_multi", "Percent_covered_multi" )] 
colnames(new) <- c("SangerSampleID", "GAinSID", "viral_genome", "Virus_id", "SRSModel", "Count_Unique", "Count_Multi", "Difference", "Percentage_Unique", "Percentage_Multi", "Coverage_bp_unique", "Percent_covered_unique","Coverage_bp_multi", "Percent_covered_multi" )
results_all <- new
results_all <- merge(results_all, key, by = "SangerSampleID")
colnames(results_all)[17] <- "Sample_timepoint"
colnames(results_all)[3] <- "GAinSID"

##-------------------------------------------------------------------------------------------------------
write.table(results, '/gpfs2/well/immune-rep/users/kvi236/VIRUS/Sepsis_Scripts/results_unique_longformat.txt', sep='\t')
write.table(results_multi, '/gpfs2/well/immune-rep/users/kvi236/VIRUS/Sepsis_Scripts/results_multi_longformat.txt', sep='\t')
write.table(results_all, '/gpfs2/well/immune-rep/users/kvi236/VIRUS/Sepsis_Scripts/results_all_longformat.txt', sep='\t')
results_unique <- read.delim('/gpfs2/well/immune-rep/users/kvi236/VIRUS/Sepsis_Scripts/results_unique_longformat.txt', sep='\t')
results_multi <- read.delim('/gpfs2/well/immune-rep/users/kvi236/VIRUS/Sepsis_Scripts/results_multi_longformat.txt', sep='\t')
results_all <- read.delim('/gpfs2/well/immune-rep/users/kvi236/VIRUS/Sepsis_Scripts/results_all_longformat.txt', sep='\t')
##-----------------------------------------------------------------------------------------------------


## Investigating Coverage 
## Looking into in more detail overal coverage of viruses in all 903 samples  
paths <- list.files('/gpfs2/well/immune-rep/users/kvi236/VIRUS/Sepsis_Scripts/Plots/Coverage_unique', full.names=TRUE)
paths <- grep("*_unique_read_alignment", paths, value=TRUE)

investigate <- data.frame()
for(i in 1:length(paths)){
	path <- paths[i]
	coverage_NC <- read.delim(path, sep='\t')
	genome_length <- length(coverage_NC[,1])
	coverage_bases_0 <- rowSums(coverage_NC)
	coverage_positive <- length(coverage_bases_0[coverage_bases_0 > 0])
	coverage_percent <- coverage_positive / genome_length * 100
	a <- unlist(str_split(path, "/"))[11]
	a <- unlist(str_split(a, "_unique"))[1]
	virus_id <- a
	Virus_eng <- Virus_database$Virus_name[Virus_database$viral_genome == virus_id]
    new_row0 <- data.frame(virus_id, Virus_eng, coverage_percent)
	colnames(new_row0) <- c("virus_id", "viral_name","percent_genome_covered")
    try <- merge(complete, new_row0, by="viral_name")
	investigate <- rbind(try, investigate)
}


paths <- list.files('/gpfs2/well/immune-rep/users/kvi236/VIRUS/Sepsis_Scripts/Plots/Coverage_multi', full.names=TRUE)
paths <- grep("*_multi_read_alignment", paths, value=TRUE)

investigate_multi <- data.frame()
for(i in 1:length(paths)){
	path <- paths[i]
	coverage_NC <- read.delim(path, sep='\t')
	genome_length <- length(coverage_NC[,1])
	coverage_bases_0 <- rowSums(coverage_NC)
	coverage_positive <- length(coverage_bases_0[coverage_bases_0 > 0])
	coverage_percent <- coverage_positive / genome_length * 100
	a <- unlist(str_split(path, "/"))[11]
	a <- unlist(str_split(a, "_multi"))[1]
	virus_id <- a
	Virus_eng <- Virus_database$Virus_name[Virus_database$viral_genome == virus_id]
    new_row0 <- data.frame(virus_id, Virus_eng, coverage_percent)
	colnames(new_row0) <- c("virus_id", "viral_name","percent_genome_covered")
    try <- merge(complete, new_row0, by="viral_name")
	investigate_multi <- rbind(try, investigate_multi)
}


whole_genome_coverage <- merge(investigate, investigate_multi, by='viral_name')
whole_genome_coverage <- whole_genome_coverage[,c('viral_name', 'cases_10.x', 'percent_10.x', 'cases_1.x', 'percent_1.x', 'cases_30.x', 'percent_30.x', 'cases_50.x', 'percent_50.x', 'percent_genome_covered.x', 'percent_genome_covered.y')]
colnames(whole_genome_coverage) <- c('viral_name', 'cases_10_unique', 'percent_10_unique', 'cases_1_unique', 'percent_1_unique', 'cases_30_unique', 'percent_30_unique', 'cases_50_unique', 'percent_50_unique', 'percent_genome_covered_unique', 'percent_genome_covered_multi')

#-------------------------------------------------------------------------------------------------------
write.table(whole_genome_coverage, '/gpfs2/well/immune-rep/users/kvi236/VIRUS/Sepsis_Scripts/results_coverage_shortformat.txt', sep='\t')
whole_genome_coverage <- read.delim('/gpfs2/well/immune-rep/users/kvi236/VIRUS/Sepsis_Scripts/results_coverage_shortformat.txt', sep='\t')
##-----------------------------------------------------------------------------------------------------
# ALL SUMMARY TABLES MADE 
# ------------------------------------------------------------------
##-------------------------------------------------------------------
##-------------------------------------------------------------------
##-------------------------------------------------------------------
## PART TWO 
## ANALYSIS 
library(ggforce)
library(ggExtra)
## Can we begin to seperate out cases / non cases

## One read length (calculate average across data set for mapped reads) - mapping is the same for unique and multi so 
## don't need to run twice!
paths <- list.dirs('/well/immune-rep/users/kvi236/VIRUS/Sepsis_Bulk/ViRNA_SEQ_Multi_Mapping_Analysis')
paths <- grep("Reports", paths, value=TRUE)
source('/gpfs2/well/immune-rep/users/kvi236/VIRUS/Sepsis_Scripts/auxillary_functions.R')
results$SRSModel <- as.factor(results$SRSModel)
mapped_read_length <- c()
for(i in 1:length(paths)){
	report_path <- paths[i]
	log_paths <- list.files(report_path, full.names=TRUE)
	log_path <- grep("Log.final.out", log_paths, value=TRUE)
	s <- Extraction_Log_final(log_path)
	length_mapped_reads <- s$Length_vector[1]
	mapped_read_length <- c(mapped_read_length,length_mapped_reads)
}
average_mapped_read_length <- mean(mapped_read_length)

## coverage per each genome corresponding to one mapped read 
Virus_lengths$percent_one_read <- average_mapped_read_length  / Virus_lengths$Genome_length *100
Virus_lengths$percent_two_read <- (average_mapped_read_length *2)  / Virus_lengths$Genome_length *100
colnames(Virus_lengths)[1] <- "viral_genome"
hlines_df <- Virus_lengths[,c("viral_genome","percent_one_read","percent_two_read" )]
results_all <- merge(results_all, hlines_df, by="viral_genome")

results_all$SRSModel <- as.factor(results_all$SRSModel)
## Plot read counts against percentage genome covered 
pdf("Plots/coverage_counts_unique.pdf", width=10, height=10)
for(i in 1:23){
	p<-ggplot(results_all, aes(x=Count_Unique, y=Percent_covered_unique, col=SRSModel)) + geom_point() + facet_wrap_paginate(~viral_genome, ncol = 2, nrow = 2, page = i, scales = "free") +xlab("Unique Mapping Read Count") +ylab("Unique Mapping Percent Genome Covered")+ geom_smooth(method=lm, aes(col=NULL)) + geom_rug(col=rgb(.5,0,0,alpha=.2)) +theme_classic()  +  geom_hline(aes(yintercept = percent_one_read), color="red")+  geom_hline(aes(yintercept = percent_two_read), color="green")
	plot(p)
}
for(i in 1:23){
	p<-ggplot(results_all, aes(x=Count_Multi, y=Percent_covered_multi, col=SRSModel)) + geom_point() + facet_wrap_paginate(~viral_genome, ncol = 2, nrow = 2, page = i, scales = "free") +xlab("Multi Mapping Read Count") +ylab("Multi Mapping Percent Genome Covered")+ geom_smooth(method=lm, aes(col=NULL)) + geom_rug(col=rgb(.5,0,0,alpha=.2)) +theme_classic()  +  geom_hline(aes(yintercept = percent_one_read), color="red")+  geom_hline(aes(yintercept = percent_two_read), color="green")
	plot(p)
}
dev.off()

pdf("Plots/counts_uniquevsmulti.pdf", width=10, height=10)
for(i in 1:23){
	p<-ggplot(results, aes(x=Count_Unique, y=Count_Multi, col=SRSModel)) + geom_point() + facet_wrap_paginate(~viral_genome, ncol = 2, nrow = 2, page = i, scales = "free") +xlab("Unique Mapping Read Count") +ylab("Multi Mapping Read Count")+ geom_smooth(method=lm, aes(col=NULL)) + geom_rug(col=rgb(.5,0,0,alpha=.2)) +theme_classic()  
	plot(p)
}
dev.off()
##-----------------------------------------------------------
##---------------------------------------------------------

# Viral Families
pappilloma <- unique(grep("papillomavirus", whole_genome_coverage$viral_name, value=TRUE))
hepatitis  <- unique(grep("Hepatitis", whole_genome_coverage$viral_name, value=TRUE))
retrovirus <- unique(grep("retrovirus", whole_genome_coverage$viral_name, value=TRUE))
torque_teno <- unique(grep("teno", whole_genome_coverage$viral_name, value=TRUE))
herpes <- unique(grep("herpes", whole_genome_coverage$viral_name, value=TRUE))

#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
##cases dataframe = new_format
library(gridExtra)
pdf('Plots/Herpes_SUMMARY.pdf', width=7.5, height=5)
p1<-ggplot(whole_genome_coverage[whole_genome_coverage$viral_name %in% herpes,] , aes(x=viral_name, y=percent_genome_covered_unique, fill=percent_1_unique)) + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + geom_hline(yintercept=5)+ labs(fill="Percentage Positive") + xlab("Herpes Virus") + ylab("Percent genome covered by unique reads")
p2<-ggplot(whole_genome_coverage[whole_genome_coverage$viral_name %in% herpes,], aes(x=viral_name, y=percent_1_unique)) + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +xlab("Herpes Virus") + ylab("Positive Individuals (1 or more reads) %")
p3 <-ggplot(results_all[results_all$viral_genome %in% herpes,], aes(x=viral_genome, y=Count_Unique, color=SRSModel)) + 
  geom_boxplot() + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + labs(fill="SRS Model") +xlab("Herpes Virus") + ylab("Number of Unique Reads") + facet_zoom(ylim = c(0, 200))
grid.arrange(p1, p2, ncol=2)
plot(p3)
dev.off()

pdf('Plots/hepatitis_SUMMARY.pdf', width=7.5, height=5)
p1<-ggplot(whole_genome_coverage[whole_genome_coverage$viral_name %in% hepatitis,] , aes(x=viral_name, y=percent_genome_covered_unique, fill=percent_1_unique)) + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + geom_hline(yintercept=5)+ labs(fill="Percentage Positive") + xlab("Hepatitis Virus") + ylab("Percent genome covered by unique reads")
p2<-ggplot(whole_genome_coverage[whole_genome_coverage$viral_name %in% hepatitis,], aes(x=viral_name, y=percent_1_unique)) + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +xlab("Hepatitis Virus") + ylab("Positive Individuals (1 or more reads) %")
p3 <-ggplot(results_all[results_all$viral_genome %in% hepatitis,], aes(x=viral_genome, y=Count_Unique, color=SRSModel)) + 
  geom_boxplot() + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + labs(fill="SRS Model") +xlab("Hepatitis Virus") + ylab("Number of Unique Reads") + facet_zoom(ylim = c(0, 20))
grid.arrange(p1, p2, ncol=2)
plot(p3)
dev.off()

pdf('Plots/pappilloma_SUMMARY.pdf', width=7.5, height=5)
p1<-ggplot(whole_genome_coverage[whole_genome_coverage$viral_name %in% pappilloma,] , aes(x=viral_name, y=percent_genome_covered_unique, fill=percent_1_unique)) + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + geom_hline(yintercept=5)+ labs(fill="Percentage Positive") + xlab("Pappilloma Virus") + ylab("Percent genome covered by unique reads")
p2<-ggplot(whole_genome_coverage[whole_genome_coverage$viral_name %in% pappilloma,], aes(x=viral_name, y=percent_1_unique)) + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +xlab("Pappilloma Virus") + ylab("Positive Individuals (1 or more reads) %")
p3 <-ggplot(results_all[results_all$viral_genome %in% pappilloma,], aes(x=viral_genome, y=Count_Unique, color=SRSModel)) + 
  geom_boxplot() + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + labs(fill="SRS Model") +xlab("Pappilloma Virus") + ylab("Number of Unique Reads") + facet_zoom(ylim = c(0, 20))
grid.arrange(p1, p2, ncol=2)
plot(p3)
dev.off()

pdf('Plots/torque_teno_SUMMARY.pdf', width=7.5, height=5)
p1<-ggplot(whole_genome_coverage[whole_genome_coverage$viral_name %in% torque_teno,] , aes(x=viral_name, y=percent_genome_covered_unique, fill=percent_1_unique)) + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + geom_hline(yintercept=5)+ labs(fill="Percentage Positive") + xlab("Torque-teno Virus") + ylab("Percent genome covered by unique reads")
p2<-ggplot(whole_genome_coverage[whole_genome_coverage$viral_name %in% torque_teno,], aes(x=viral_name, y=percent_1_unique)) + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +xlab("Torque-teno Virus") + ylab("Positive Individuals (1 or more reads) %")
p3 <-ggplot(results_all[results_all$viral_genome %in% torque_teno,], aes(x=viral_genome, y=Count_Unique, color=SRSModel)) + 
  geom_boxplot() + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + labs(fill="SRS Model") +xlab("Torque-teno Virus Virus") + ylab("Number of Unique Reads") + facet_zoom(ylim = c(0, 20))
grid.arrange(p1, p2, ncol=2)
plot(p3)
dev.off()

##---------------------------------------------------------
##Looking at EBV adding Cyndis Data

ddPCR <- read.delim('/gpfs2/well/immune-rep/users/kvi236/VIRUS/Sepsis_Scripts/Sample_Sheets/EBV_reactivated_samples_ddpcr.csv', sep=",")
plasma_EBV <- read.delim('/gpfs2/well/immune-rep/users/kvi236/VIRUS/Sepsis_Scripts/Sample_Sheets/EBV_reactivated_samples.csv', sep=",")
ebv_analysis_plasma <- merge(results_all, plasma_EBV, by="SampleID")
ebv_analysis_plasma <- ebv_analysis_plasma[ebv_analysis_plasma$viral_genome=="Human gammaherpesvirus 4",]
ebv_analysis_ddpcr <- merge(results_all, ddPCR, by="SampleID")
ebv_analysis_ddpcr <- ebv_analysis_ddpcr[ebv_analysis_ddpcr$viral_genome=="Human gammaherpesvirus 4",]
ebv_analysis_ddpcr$EBV_positive_ddpcr <- as.factor(ebv_analysis_ddpcr$EBV_positive_ddpcr)
ebv_analysis_plasma$EBV_positive <- as.factor(ebv_analysis_plasma$EBV_positive)


pdf('Plots/ebv_ddPCR.pdf', width=7.5, height=5)
p1<-ggplot(ebv_analysis_ddpcr, aes(x=Count_Multi, y=ddPCR, color=EBV_positive_ddpcr, shape=SRSModel)) +
+   geom_point() + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + labs(fill="ddEBV Positive") + xlab("RNAseq PBMC MultiMapping Read Count") + ylab("PCR Read quantificication") + facet_zoom(xlim = c(0, 100), ylim = c(0, 500))+ geom_rug(col=rgb(.5,0,0,alpha=.2))
plot(p1)
p1<-ggplot(ebv_analysis_plasma, aes(x=Count_Multi, y=EBV_n_reads_all, color=EBV_positive, shape=SRSModel)) +
+   geom_point() + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  + xlab("MultiMapping Read Count") + ylab("Plasma Read Count") + facet_zoom(xlim = c(0, 100), ylim = c(0, 500)) + geom_rug(col=rgb(.5,0,0,alpha=.2))
plot(p1)
p1<-ggplot(ebv_analysis_plasma, aes(x=EBV_positive, y=Count_Multi, fill=SRSModel)) +
+   geom_boxplot() + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  + xlab("EBV Diagnosis from Plasma") + ylab("RNAseq PBMC MultiMapping Read Count") + facet_zoom(ylim = c(0, 100))
plot(p1)
p1<-ggplot(ebv_analysis_ddpcr, aes(x=EBV_positive_ddpcr, y=Count_Multi, fill=SRSModel)) +
+   geom_boxplot() + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  + xlab("EBV Diagnosis from ddPCR") + ylab("RNAseq PBMC MultiMapping Read Count") + facet_zoom(ylim = c(0, 100))
plot(p1)
dev.off()


write.table(ebv_analysis_ddpcr, '/gpfs2/well/immune-rep/users/kvi236/VIRUS/Sepsis_Scripts/Counts/ebv_analysis_ddpcr.txt', sep='\t')
write.table(ebv_analysis_plasma, '/gpfs2/well/immune-rep/users/kvi236/VIRUS/Sepsis_Scripts/Counts/ebv_analysis_plasma.txt', sep='\t')














multi_mapping_counts_pca <- multi_mapping_counts[,c("Virus_name", "Sample", "count")]

srs <- read.delim('/gpfs2/well/immune-rep/users/kvi236/VIRUS/Sepsis_Scripts/Sample_Sheets/SRS_model_all_samples_864.txt', sep="\t")
key <- read.delim('/gpfs2/well/immune-rep/users/kvi236/VIRUS/Sepsis_Scripts/Sample_Sheets/Sample_key.txt', sep="\t")

## Add srs information
colnames(multi_mapping_counts_pca)[colnames(multi_mapping_counts_pca)=="Sample"]  <- "SangerSampleID"
multi_mapping_counts_pca <- merge(multi_mapping_counts_pca, key, by="SangerSampleID", all.x=TRUE)
multi_mapping_counts_pca$SampleID <- multi_mapping_counts_pca$GAinSID
multi_mapping_counts_pca <- merge(multi_mapping_counts_pca, srs, by="SampleID", , all.x=TRUE)
multi_mapping_counts_pca$count <- as.numeric(multi_mapping_counts_pca$count)
multi_mapping_counts_pca <- multi_mapping_counts_pca[, c("SangerSampleID", "Virus_name", "count", "SRSModel")

data_wide <- spread(multi_mapping_counts_pca, Virus_name, count)
data_wide_1 <- data_wide[, c(-1, -2)]
data_wide_1  <- mutate_all(data_wide_1, function(x) as.numeric(x))
data_wide_1$Sample <- data_wide$Sample
data_wide_1$SRSModel <- as.factor(data_wide$SRSModel)
for_pca <- data_wide_1[,c(-91, -92)]
data_wide.pca <- prcomp(for_pca[ , which(apply(for_pca, 2, var) != 0)] , center = TRUE, scale=TRUE)
summary(data_wide.pca)

library(ggbiplot)
pdf("pca.pdf", width=20, height=20)
autoplot(data_wide.pca, data = data_wide_1, colour = 'SRSModel', label = TRUE, loadings = TRUE, loadings.colour = 'blue', loadings.label = TRUE, loadings.label.size = 3)
dev.off()



#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


paths <- list.files('/gpfs2/well/immune-rep/users/kvi236/VIRUS/Sepsis_Scripts/Plots/Coverage', full.names=TRUE)
paths <- grep(".txt", paths, value=TRUE)
path <- '/gpfs2/well/immune-rep/users/kvi236/VIRUS/Sepsis_Scripts/Plots/Coverage/NC_009334.1_multi.txt'
investigate <- data.frame()
for(i in 1:length(paths)){
	path <- paths[i]
	coverage_NC <- read.delim(path, sep='\t')
	genome_length <- length(coverage_NC[,1])
	
	coverage_bases_0 <- rowSums(coverage_NC)
	coverage_positive <- length(coverage_bases_0[coverage_bases_0 > 0])
	coverage_percent <- coverage_positive / genome_length * 100
	

	a <- unlist(str_split(path, "/"))[11]
	a <- unlist(str_split(a, "_multi"))[1]
	virus_id <- a
	Virus_eng <- Virus_database$Virus_name[Virus_database$viral_genome == virus_id]
    new_row0 <- data.frame(virus_id, Virus_eng, coverage_percent)
	colnames(new_row0) <- c("virus_id", "viral_name","percent_genome_covered")
    try <- merge(complete, new_row0, by="viral_name")
	investigate <- rbind(try, investigate)
}


lapply(coverage_NC, FUN = sum()


per_sample_coverage <- data.frame(sapply(coverage_NC, function(x) count(x>0)))
colnames(per_sample_coverage) <- "BP_covered"
per_sample_coverage$percent_covered <- per_sample_coverage$BP_covered / genome_length

per_sample_coverage_positive <- per_sample_coverage[per_sample_coverage$BP_covered >0,]
pdf("plottyplott.pdf")
ggplot(per_sample_coverage, aes(y=BP_covered)) + geom_boxplot() + theme_minimal()  + stat_summary(fun=mean, geom="point", shape=23, size=2, col="red")
ggplot(per_sample_coverage, aes(y=percent_covered)) + geom_boxplot()+ theme_minimal() + stat_summary(fun=mean, geom="point", shape=23, size=2, col="red")
ggplot(per_sample_coverage_positive, aes(y=BP_covered)) + geom_boxplot()+ theme_minimal() + stat_summary(fun=mean, geom="point", shape=23, size=2, col="red")
ggplot(per_sample_coverage_positive, aes(y=percent_covered)) + geom_boxplot()+ theme_minimal() + stat_summary(fun=mean, geom="point", shape=23, size=2, col="red")
dev.off()
