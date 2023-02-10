#-----------------------------------------------------------------
Extraction_Log_final = function(path_to_Log_file) {
  #Loading the final log file from STAR
  Log_file =read.delim(path_to_Log_file,header = F,sep="\t")
  Name_variable = as.character(Log_file$V1)
  Value_variable = as.character(Log_file$V2)
  
  #Extracting the informations about the mapping quality 
  Uniquely_mapped_percent = Value_variable[Name_variable=="                        Uniquely mapped reads % |"]
  Multiple_mapped_percent = Value_variable[Name_variable=="             % of reads mapped to multiple loci |"]
  Unmapped_too_short_percent = Value_variable[Name_variable=="                 % of reads unmapped: too short |"]
  Unmapped_mismatch_percent = Value_variable[Name_variable=="       % of reads unmapped: too many mismatches |"]
  Unmapped_other_percent = Value_variable[Name_variable=="                     % of reads unmapped: other |"]
  
  remove_percent = function(x) {
    l = nchar(x)
    x = substr(x,1,l-1)
    x=as.numeric(x)
    return(x)
  }  
  
  Uniquely_mapped_percent=remove_percent(Uniquely_mapped_percent)
  Multiple_mapped_percent=remove_percent(Multiple_mapped_percent)
  Unmapped_too_short_percent=remove_percent(Unmapped_too_short_percent)
  Unmapped_mismatch_percent = remove_percent(Unmapped_mismatch_percent)
  Unmapped_other_percent = remove_percent(Unmapped_other_percent)
  Total_unmapped = Unmapped_too_short_percent + Unmapped_mismatch_percent + Unmapped_other_percent
  
  Matrix_mapping = matrix(c(Uniquely_mapped_percent,Multiple_mapped_percent,Total_unmapped),nrow = 3)
  
  #####Looking at length of mapping, deletion and insertion
  
  Mean_mapped_length = as.numeric(Value_variable[Name_variable=="                          Average mapped length |"])
  Mean_deletion_length = as.numeric(Value_variable[Name_variable=="                        Deletion average length |"])
  Mean_insertion_length = as.numeric(Value_variable[Name_variable=="                       Insertion average length |"])
  #####Looking at rates of  of mismatch, insertion and deletion
  
  Mismatch_rate= remove_percent(Value_variable[Name_variable=="                      Mismatch rate per base, % |"])
  Deletion_rate = remove_percent(Value_variable[Name_variable=="                         Deletion rate per base |"])
  Insertion_rate = remove_percent(Value_variable[Name_variable=="                        Insertion rate per base |"])
  
  List_elements = list(Mapping_result= Matrix_mapping , 
                       Length_vector = c(Mean_mapped_length,Mean_insertion_length,Mean_deletion_length),
                       Rate_vector = c(Mismatch_rate,Insertion_rate,Deletion_rate) )
  
  return(List_elements)
  }


##PARALISED COVERAGE FUNCION
map_viral_coverage_multi <- function(virus_id, directory){
	cl <- 7
	paths <- list.dirs(directory)
    paths <- grep("Reports", paths, value=TRUE)
	registerDoParallel(cl)
	z <- virus_id
	q <- foreach(i = 1:length(paths), .combine=cbind, .packages=c('GenomicAlignments', 'GenomicRanges', 'tidyverse')) %dopar% {
			bam_path <- list.files(paths[i], full.names=TRUE)
			bam_path <- grep("for_feature_counts.bam", bam_path, value=TRUE)
			print(paste0('Reading in file ', i, ' of 903'))
			x <- readGAlignments(bam_path, param=ScanBamParam(what=scanBamWhat()))
			# Subset to the reads analysed (multi)
			x= x[x@elementMetadata$flag !=355 & x@elementMetadata$flag != 419 & x@elementMetadata$flag != 403 & x@elementMetadata$flag != 339 ]
			xcov <- coverage(x)
			xnum <- as.numeric(xcov[[z]])
			l <- paths[i]
			a <- unlist(str_split(l, "/ViRNA_SEQ_Multi_Mapping_Analysis/"))[2]
			a <- unlist(str_split(a, ".Unmapped.out/Reports"))[1]
			xnum <- c(a, xnum)
				}
	colnames(q) <- q[1,]
	q <- q[-1,]
	coverage_df <- data.frame(q)
	bp <- as.numeric(rownames(coverage_df))
	virus <- z 
	Virus_eng <- Virus_database$Virus_name[Virus_database$viral_genome == virus]
	pdf_name <- paste0('Plots/Coverage_multi/', virus, "_multi.pdf")
	print("Making matplot")
	pdf(pdf_name, width=15, height=15)
	matplot(x=bp, y = coverage_df, type = 'l', lty = 1, xlab=paste0(virus, " / ", Virus_eng), ylab="Read Coverage", xaxt="none", main=paste0("Read coverage of ", z, " in GAinS RNAseq"), cex.main=3) 
	axis(side = 1, at=seq(0, length(bp), by=round((length(bp)/20))), cex=30)
	dev.off()
	write.table(coverage_df, paste0('Plots/Coverage_multi/', virus, '_multi_read_alignment.txt'), sep='\t')
	coverage_function <- function(x){
	  length(x[x>0])
	}
	srs <- read.delim('/gpfs2/well/immune-rep/users/kvi236/VIRUS/Sepsis_Scripts/Sample_Sheets/SRS_model_all_samples_864.txt', sep="\t")
	key <- read.delim('/gpfs2/well/immune-rep/users/kvi236/VIRUS/Sepsis_Scripts/Sample_Sheets/Sample_key.txt', sep="\t")
	print("Merging Dataframes")
	individual_coverage <- data.frame(apply(coverage_df, 2,  coverage_function))
	colnames(individual_coverage) <- "Coverage_bp"
	individual_coverage$Percent_covered <- individual_coverage$Coverage_bp / length(bp) * 100
	individual_coverage$SangerSampleID<- rownames(individual_coverage)
	individual_coverage <- merge(key, individual_coverage, by="SangerSampleID", all.y = TRUE)
	colnames(individual_coverage)[2] <- "SampleID"
	individual_coverage <- merge(srs, individual_coverage, by="SampleID", , all.y = TRUE)
	individual_coverage$SRSModel <- factor(individual_coverage$SRSModel)
	write.table(individual_coverage, paste0('Plots/Coverage_multi/', virus, '_multi_individual_coverage.txt'), sep='\t')
	pdf_name <- paste0('Plots/Coverage_multi/', virus, "_multi_coverage_srs.pdf")
	print("Making Coverage Plot")
	pdf(pdf_name, width=15, height=15)
	p1<- ggplot(individual_coverage, aes(x=SRSModel, y=Percent_covered, color=SRSModel)) + geom_boxplot() + theme_classic() +ggtitle(paste0("Distribution of Genome Coverage in GAinS RNAseq for ", Virus_eng)) +xlab("SRS assignment") +ylab("Percent Genome covered %") 
	plot(p1)
	dev.off()
	gc()
}


map_viral_coverage_unique <- function(virus_id, directory){
	cl <- 7
	registerDoParallel(cl)
	paths <- list.dirs(directory)
    paths <- grep("Reports", paths, value=TRUE)
	z <- virus_id
	q <- foreach(i = 1:length(paths), .combine=cbind, .packages=c('GenomicAlignments', 'GenomicRanges', 'tidyverse')) %dopar% {
			bam_path <- list.files(paths[i], full.names=TRUE)
			bam_path <- grep("for_feature_counts.bam", bam_path, value=TRUE)
			print(paste0('Reading in file ', i, ' of 903'))
			x <- readGAlignments(bam_path, param=ScanBamParam(what=scanBamWhat()))
			# Subset to the reads analysed (multi)
			x= x[x@elementMetadata$flag !=355 & x@elementMetadata$flag != 419 & x@elementMetadata$flag != 403 & x@elementMetadata$flag != 339 ]
			x= x[x@elementMetadata$mapq==255]
			xcov <- coverage(x)
			xnum <- as.numeric(xcov[[z]])
			l <- paths[i]
			a <- unlist(str_split(l, "/ViRNA_SEQ_UNIQUE_Mapping_Analysis/"))[2]
			a <- unlist(str_split(a, ".Unmapped.out/Reports"))[1]
			xnum <- c(a, xnum)
				}
	colnames(q) <- q[1,]
	q <- q[-1,]
	coverage_df <- data.frame(q)
	bp <- as.numeric(rownames(coverage_df))
	virus <- z 
	Virus_eng <- Virus_database$Virus_name[Virus_database$viral_genome == virus]
	pdf_name <- paste0('Plots/Coverage_unique/', virus, "_unique.pdf")
	print("Making matplot")
	pdf(pdf_name, width=15, height=15)
	matplot(x=bp, y = coverage_df, type = 'l', lty = 1, xlab=paste0(virus, " / ", Virus_eng), ylab="Read Coverage", xaxt="none", main=paste0("Read coverage of ", z, " in GAinS RNAseq"), cex.main=3) 
	axis(side = 1, at=seq(0, length(bp), by=round((length(bp)/20))), cex=30)
	dev.off()
	write.table(coverage_df, paste0('Plots/Coverage_unique/', virus, '_unique_read_alignment.txt'), sep='\t')
	coverage_function <- function(x){
	  length(x[x>0])
	}
	srs <- read.delim('/gpfs2/well/immune-rep/users/kvi236/VIRUS/Sepsis_Scripts/Sample_Sheets/SRS_model_all_samples_864.txt', sep="\t")
	key <- read.delim('/gpfs2/well/immune-rep/users/kvi236/VIRUS/Sepsis_Scripts/Sample_Sheets/Sample_key.txt', sep="\t")
	print("Merging Dataframes")
	individual_coverage <- data.frame(apply(coverage_df, 2,  coverage_function))
	colnames(individual_coverage) <- "Coverage_bp"
	individual_coverage$Percent_covered <- individual_coverage$Coverage_bp / length(bp) * 100
	individual_coverage$SangerSampleID<- rownames(individual_coverage)
	individual_coverage <- merge(key, individual_coverage, by="SangerSampleID", all.y = TRUE)
	colnames(individual_coverage)[2] <- "SampleID"
	individual_coverage <- merge(srs, individual_coverage, by="SampleID", , all.y = TRUE)
	individual_coverage$SRSModel <- factor(individual_coverage$SRSModel)
	write.table(individual_coverage, paste0('Plots/Coverage_unique/', virus, '_unique_individual_coverage.txt'), sep='\t')
	pdf_name <- paste0('Plots/Coverage_unique/', virus, "_unique_coverage_srs.pdf")
	print("Making Coverage Plot")
	pdf(pdf_name, width=15, height=15)
	p1<- ggplot(individual_coverage, aes(x=SRSModel, y=Percent_covered, color=SRSModel)) + geom_boxplot() + theme_classic() +ggtitle(paste0("Distribution of Genome Coverage in GAinS RNAseq for ", Virus_eng)) +xlab("SRS assignment") +ylab("Percent Genome covered %") 
	plot(p1)
	dev.off()
	gc()
}



