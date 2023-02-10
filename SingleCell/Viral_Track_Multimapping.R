## ---------------------------
## Script name: VIRAL TRACK: Viral_Scanning: Module 1 - MAPPING
## Function: Map Single Cell Virals from individual FASTQ file using STAR.
## Author: Pierre Bost (as used in Viral TRACK paper). Updated by Lauren Overend (LEO) - Wellcome Trust Centre for Human Genetics
##
## Date Created: 2020 - June 
##
## Email: lauren.overend@oriel.ox.ac.uk
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------
if(nzchar(system.file(package = "optparse"))==FALSE){
  stop("optparse not installed. Terminating. \n")
}
suppressMessages(library(optparse))
parser <- OptionParser()
option_list <- list( 
  make_option(c("-n", "--nThreadmap"), action="store", default=8, type="integer", help="runThreadN for Star Mapping. Note will also be used as threads for Feature Counts [default]"),
  make_option(c("-o", "--outputdir"), action="store", default='/gpfs2/well/immune-rep/users/kvi236/LCL_VIRALTRACK', type="character", help="Path to output directory"),
  make_option(c("-i", "--indexgenome"), action="store", type="character", default="/well/immune-rep/users/kvi236/References/VIRAL_TRACK_REFERENCE_BUILD_273a", help="Path to VIRAL TRACK reference genome [default]"),
  make_option(c("-s", "--nThreadsort"), action="store", type="integer", default=1, help="outBAMsortingThreadN for STAR Mapping [default] - usually < runThreadN"),
  make_option(c("-m", "--minreads"), action="store", type="integer", default=50, help="Minimum number of mapped viral reads [default]"),
  make_option(c("-b", "--bins"), action="store", type="integer", default=50, help="outBAMsortingBinsN for STAR Mapping [default]"),
  make_option(c("-f", "--fastq"), action="store", type="character", default = '/well/immune-rep/shared/10X_GENOMICS/EBV_LCLS/FASTQ/SRR8427168/DOWNSAMPLES/sub100k.fa', help="Path to input FASTQ file [default]"),
  make_option(c("-r", "--runname"), action="store", type="character", default="Viral_Track", help="Run Name [default]"),
  make_option(c("-v", "--viralannotation"), action="store", type="character", default="/gpfs2/well/immune-rep/users/kvi236/References/Updated_VirusSite_Reference.txt", help="Path to VirusSite annotation file [default]"),
  make_option(c("-a", "--auxfunctions"), action="store", type="character", default="/gpfs2/well/immune-rep/users/kvi236/ViralTrackProgram/RPipeline/Viral-Track/AuxillaryFunctions/auxillary_viral_track_functions.R", help="Path to ViralTrack Auxillary Functions [default]"),
  make_option(c("-g", "--gtffile"), action="store", type="character", default="FALSE", help="Path to GTF file. If no GTF file exists use FALSE and it will be created [default]")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, print_help_and_exit = TRUE, args = commandArgs(trailingOnly = TRUE) )

if (is.null(opt$fastq)){
  stop("No FASTq file provided. Terminating. \n", call.=FALSE)
}
if (is.null(opt$outputdir)){
  stop("Output directory not present: must be specified. Terminating. \n", call.=FALSE)
}
if (opt$runname == "Viral_Track"){
  warning("Using Default Run-Name: Reconsider using Custom Run-Name. \n")
}
if (!dir.exists(opt$outputdir)) {
  warning("Output directory does not exist! Creating it! \n")
  dir.create(opt$outputdir)
}
if (!is.numeric(opt$nThreadmap) | opt$nThreadmap < 1 ) {
  stop("STAR Mapping Threads are Incorrect. Terminating. \n")
}
if (!is.numeric(opt$nThreadsort)) {
  stop("STAR BAM Sorting Threads are incorrect. Terminating. \n")
}
if (opt$nThreadsort > 1 ) {
  warning("STAR BAM Sorting Threads are higher than recommended (Default: 1). Consider adjusting paramaters.  \n")
}
if (!dir.exists(opt$indexgenome)) {
  stop("Index Genome Directory Does Not Exist. Terminating. \n")
}
if (!file.exists(opt$viralannotation)) {
  stop("VirusSite Database File Does Not Exist. Terminating. \n")
}
if (opt$nThreadsort > opt$nThreadmap ) {
  warning("Bam Sorting Threads exceeds Bam Mapping Threads. May be incompatible with memory request. \n")
}
## Check all packages are installed. 
if(nzchar(system.file(package = "Biostrings"))==FALSE){
  stop("Biostrings not installed. Terminating. \n")
} 
if(nzchar(system.file(package = "ShortRead"))==FALSE){
  stop("ShortRead not installed. Terminating. \n")
}
if(nzchar(system.file(package = "doParallel"))==FALSE){
  stop("doParallel not installed. Terminating. \n")
}
if(nzchar(system.file(package = "GenomicAlignments"))==FALSE){
  stop("GenomicAlignments not installed. Terminating. \n")
}
if(nzchar(system.file(package = "Matrix"))==FALSE){
  stop("Matrix not installed. Terminating. \n")
}
## Check Validiamety of Input FASTQ File
List_target_path = c()
if (!is.null(opt$fastq)) {
  if(file.exists(opt$fastq)){
    cat("FASTQ File present. \n") 
    if(any(grepl(".fa|.fq|.fasta", opt$fastq))==TRUE){
      cat("FASTQ File Type is Valid. \n")
      List_target_path = opt$fastq
    } else {
      stop("Fastq File Provided is Not of Type '.fasta/.fq/.fa'. Terminating. \n")
    }
  } else {
    stop("FASTQ Provided but Path is Invalid. Terminating. \n")
  }
} else {
   stop("No FASTQ File Provided. Terminating. \n")
}

## Load Required Libraries
suppressMessages(library(Biostrings))
suppressMessages(library(ShortRead))
suppressMessages(library(doParallel))
suppressMessages(library(GenomicAlignments))
suppressMessages(library(Matrix))

## Notin Function:
`%notin%` <- Negate(`%in%`)

## Setting up log.file: 
name <- unlist(strsplit(opt$fastq,"/",fixed = T))
sample_name <- name[length(name)]
sample_name = gsub('.fastq|.fa|.fq|.gz','',sample_name) 
log <-  paste0(opt$outputdir, "/ViralTrack_MultiMapping_", sample_name, ".log")

## Checking the parameters values
cat("ViralTrack_Scanning.2.0 by Lauren Overend & Pierre Bost \n", file=log, append=TRUE)
cat(paste0("Run Name: ", opt$runname, "\n"), file=log, append=TRUE)
start_time_1 <- Sys.time()
cat(paste0("Start time: ", start_time_1, "\n"), file=log, append=TRUE)
cat(paste0("Output directory: ", opt$outputdir, "\n"), file=log, append=TRUE)
cat(paste0("FASTQ File: ", opt$fastq, "\n"), file=log, append=TRUE)
cat("----------------------------------------------\n", file=log, append=TRUE)
cat("Input Parameters Read: Commencing Viral Track. \n", file=log, append = TRUE)
cat(paste0("--nThreadmap: ", opt$nThreadmap, "\n"), file=log, append = TRUE)
cat(paste0("--outputdir: ", opt$outputdir, "\n"), file=log, append = TRUE)
cat(paste0("--indexgenome: ", opt$indexgenome, "\n"), file=log, append = TRUE)
cat(paste0("--nThreadsort: ", opt$nThreadsort, "\n"), file=log, append = TRUE)
cat(paste0("--minreads: ", opt$minreads, "\n"), file=log, append = TRUE)
cat(paste0("--bins: ", opt$bins, "\n"), file=log, append = TRUE)
cat(paste0("--fastq: ", opt$fastq, "\n"), file=log, append = TRUE)
cat(paste0("--runname: ", opt$runname, "\n"), file=log, append = TRUE)
cat(paste0("--viralannotation: ", opt$viralannotation, "\n"), file=log, append = TRUE)
cat(paste0("--auxfunctions: ", opt$gtffile, "\n"), file=log, append = TRUE)
cat(paste0("--gtffile: ", opt$auxfunctions, "\n"), file=log, append = TRUE)
cat("----------------------------------------------\n", file=log, append=TRUE)

## Renaming Paramaters to original variable names as in ViralTrack1.0 and sourcing auxillary funtions
N_thread = opt$nThreadmap 
N_thread_sort = opt$nThreadsort
N_bins = opt$bins
Output_directory = paste0(opt$outputdir, "/Multi_Mapping_Analysis")
dir.create(Output_directory)
Name_run = opt$runname    
Index_genome = opt$indexgenome 
Minimal_read_mapped = as.numeric(opt$minreads)
Viral_annotation_file = opt$viralannotation
path_to_gtf = opt$gtffile
source(opt$auxfunctions) 
## ------------------------------------------------------------------------------------
## ------------------------------------------------------------------------------------
## Mapping: 
## Setting Up Directory: 
name_target = unlist(base::strsplit(opt$fastq,"/",fixed = T))
name_target = name_target[length(name_target)]
name_target = gsub('/','',name_target)
is_gz_file = any(grepl(pattern = ".gz",name_target))
name_target = gsub('.fastq|.fa|.fq|.gz','',name_target)   #Cleaning the name to get the original Amplification batch number
temp_output_dir = paste(Output_directory, "/", name_target,sep = "")
dir.create(temp_output_dir)

## ------------------------------------------------------------------------------------
## Making GTF 
cat("----------------------------------------------\n", file=log, append=TRUE)
if(path_to_gtf ==FALSE){
  path_to_gtf  = paste0(temp_output_dir, "/referencegtf.gtf")
  cat("NO GTF FILE PROVIDED: CREATING GTF FROM PROVIDED REFERENCE. \n", file=log, append = TRUE)
  gtf_maker(Index_genome, path_to_gtf)
  cat("GTF FILE GENERATED. \n", file=log, append = TRUE)
  cat("----------------------------------------------\n", file=log, append=TRUE)
} 

## ------------------------------------------------------------------------------------
## MAPPING 
cat(paste0("Mapping: ",name_target,".fastq file \n"), file=log, append = TRUE)
start_time <- Sys.time()
cat(paste0("Start time: ", start_time, "\n"), file=log, append=TRUE)
name_prefix = paste0(temp_output_dir, "/", name_target)
  
#We construct a complex command 
STAR_mapping_command = paste("STAR --runThreadN ",N_thread," --outBAMsortingThreadN ",N_thread_sort," --outBAMsortingBinsN ", N_bins, " --genomeDir ",Index_genome," --readFilesIn ", opt$fastq," --outSAMattributes NH HI AS nM NM XS ",
                               "--outFileNamePrefix ",name_prefix," --outSAMtype BAM SortedByCoordinate --twopassMode Basic ",
                               "--outFilterMatchNmin 35 --outMultimapperOrder Random --runRNGseed 1 --outFilterScoreMinOverLread 0.6 --outFilterMatchNminOverLread 0.6 > ", name_prefix, "_STAR_MAPPING.log", sep="")

#If the file is in the format .gz then we need to add an additional paramter :
if (is_gz_file==TRUE) {
    STAR_mapping_command = paste(STAR_mapping_command, "--readFilesCommand zcat", sep=" ") #Allowing to read .gz files
    }
## Launching STAR COMMAND
system(STAR_mapping_command)
end_time <- Sys.time()
## DOCUMENTING 
cat(paste0("Is File Gzipped: ", is_gz_file, "\n"), file=log, append=TRUE)
cat(paste0("End time: ", end_time, "\n"), file=log, append=TRUE)
cat(paste0("Mapping: ",name_target,".fastq DONE! \n"), file=log, append = TRUE)
cat("----------------------------------------------\n", file=log, append=TRUE)

### ------------------------------------------------------------------------------------
### ------------------------------------------------------------------------------------

cat("Commencing QC analysis: \n", file=log, append = TRUE)
start_time <- Sys.time()
cat(paste0("QC START TIME: ", start_time, "\n"), file=log, append=TRUE)
Virus_database = read.delim(Viral_annotation_file,header=T,sep="\t")
cat("\t 1. ViralSite annotation File loaded Sucessfully. \n", file=log, append = TRUE)

# 1. Indexing the STAR Aligned.sortedByCoord.out.bam file:
cat("\t 2. Indexing BAM: \n", file=log, append = TRUE)
temp_sorted_bam = paste0(temp_output_dir, "/", name_target, "Aligned.sortedByCoord.out.bam")#[1]
#To begin with : the ordered .BAM file need to indexed
SAMtools_indexing_command = paste("samtools index",temp_sorted_bam)
system(SAMtools_indexing_command)
cat(paste0("\t Indexing of the bam file for ",name_target," is done. \n"), file=log, append = TRUE)

# 2. Compute The number of Mapped Reads For Each Chromosome/virus
cat("\t 3. Computing The number of Mapped Reads For Each Chromosome/Virus. \n", file=log, append = TRUE)
temp_chromosome_count_path = paste(temp_output_dir,"/Count_chromosomes.txt",sep = "")
SAMtools_chromosome_count_command = paste("samtools idxstats",temp_sorted_bam,">",temp_chromosome_count_path)
system(SAMtools_chromosome_count_command)
cat(paste0("\t Indexing and Stat file for ",name_target,"is done. \n"), file=log, append = TRUE)


### ---------------------------------------------------------------------------------------------------------------
# Creating VIRAL BAM FILES 
# Loading Count_Chromosome.txt and Filtering 
temp_chromosome_count = read.table(temp_chromosome_count_path,header = F,row.names = 1)
colnames(temp_chromosome_count) = c("Chromosome_length","Mapped_reads","Unknown")

# Filtering Table removing Human sequences and the removing viruses with less than a given threshold of reads:
virus <- grep("refseq|NC_", rownames(temp_chromosome_count), value=TRUE)
temp_chromosome_count <- temp_chromosome_count[rownames(temp_chromosome_count) %in% virus, ]
temp_chromosome_count = temp_chromosome_count[temp_chromosome_count$Mapped_reads>=Minimal_read_mapped,]
no_viruses <- length(temp_chromosome_count[,1])
cat(paste0("\t 4. Identified ", no_viruses, " viruses with at least ", Minimal_read_mapped, " Reads in input Fastq. \n"), file=log, append = TRUE)


# We first create a sub-directory to export the sam files corresponding to each virus
newdir <- paste(temp_output_dir,"/Viral_BAM_files",sep = "")
dir.create(newdir) 
# Extracting Viral Read BAM file per virus into new Directory: 
cat("\t 5. Extracting BAMs files to new directories. \n", file=log, append = TRUE)
# Identify Number of Parrallel threads to use for Bam sorting depending on the input number of threads and the number of identified viruses
if (length(temp_chromosome_count[,1]) <= N_thread) {
  if(length(temp_chromosome_count[,1])==0){
    N_threadR <- 1
  } else {
    N_threadR <- length(temp_chromosome_count[,1])
  } 
  } else {N_threadR <- N_thread 
}
cat("\t 5.a  Setting up Parrallel Environment in R using: ", N_threadR, " Threads. \n", file=log, append = TRUE)
cl =makeCluster(N_threadR)
registerDoParallel(cl)
cat("\t Done. \n", file=log, append = TRUE)
# Create One Bam File Per Virus: 
foreach(i=rownames(temp_chromosome_count)) %dopar% {
  temp_export_bam_command = paste("samtools view -b ",temp_sorted_bam," \'",i,"\'"," > \'",temp_output_dir,"/Viral_BAM_files/",i,".bam\'",sep = "")
  system(temp_export_bam_command)
}
cat("\t 5.b  Viral BAM Files Extracted. \n", file=log, append = TRUE)

### ---------------------------------------------------------------------------------------------------------------
# Creating HUMAN BAM FILES
temp_chromosome_count_human = read.table(temp_chromosome_count_path,header = F,row.names = 1)
colnames(temp_chromosome_count_human) = c("Chromosome_length","Mapped_reads","Unknown")
virus_sequence_names <- grep("refseq|NC_", rownames(temp_chromosome_count_human), value=TRUE)
human_chromosomes <- rownames(temp_chromosome_count_human)[rownames(temp_chromosome_count_human) %notin% virus_sequence_names]
human_chromosomes <- human_chromosomes[human_chromosomes!="*"]
temp_chromosome_count_human = temp_chromosome_count_human[rownames(temp_chromosome_count_human)%in% human_chromosomes ,]
dir.create(paste(temp_output_dir,"/HUMAN_BAM_files",sep = ""))
foreach(i=rownames(temp_chromosome_count_human)) %dopar% {
  temp_export_bam_command = paste("samtools view -b ",temp_sorted_bam," \'",i,"\'"," > \'",temp_output_dir,"/HUMAN_BAM_files/",i,".bam\'",sep = "")
  system(temp_export_bam_command)
}
cat("\t 5.c Human BAM Files Extracted. \n", file=log, append = TRUE)
stopImplicitCluster()
cat("\t 5.d Terminating Parrallel Environment. \n", file=log, append = TRUE)
### ---------------------------------------------------------------------------------------------------------------
### ---------------------------------------------------------------------------------------------------------------

cat("\t 7. Calculating QC Metrics. \n", file=log, append = TRUE)
## Generation of the QC report:
dir <- paste(temp_output_dir,"/Viral_BAM_files/",sep = "")
if(length(list.files(dir))==0){
  cat("\t ----------------------------------------------\n", file=log, append=TRUE)
  cat(paste0("\t NO VIRUSES IDENFIFIED WITH MIN READS MAPPED >=", Minimal_read_mapped, ". \n"), file=log, append=TRUE)
  cat("\t NO QC VIRAL METRICS WILL BE CALCULTED. \n", file=log, append=TRUE)
  cat("\t ----------------------------------------------\n", file=log, append=TRUE)
  }  

if(length(list.files(dir))>0){
  cat("\t ----------------------------------------------\n", file=log, append=TRUE)
  cat(paste0("\t ", length(list.files(dir)), " VIRUSES IDENFIFIED WITH MIN READS MAPPED >=", Minimal_read_mapped, ". \n"), file=log, append=TRUE)
  cat("\t QC VIRAL METRICS WILL BE CALCULTED. \n", file=log, append=TRUE)
  cat("\t ----------------------------------------------\n", file=log, append=TRUE)
  QC_result <- NULL
## Generation of the QC repor
  for(i in 1:length(rownames(temp_chromosome_count))) {
    z <- rownames(temp_chromosome_count)[i]
    BAM_file= readGAlignments(paste(temp_output_dir,"/Viral_BAM_files/",z,".bam",sep = ""),param = ScanBamParam(what =scanBamWhat()))
    #Let's check the diversity of the reads
	#Lets use just the reads that map uniquely to calculate the statistics!!!
    Viral_reads = unique(BAM_file@elementMetadata$seq)
    Viral_reads_contents = alphabetFrequency(Viral_reads,as.prob =T )
    Viral_reads_contents = Viral_reads_contents[,c("A","C","G","T")]
    
    #calculate percentage of nucleotides in the viral reads for Read Entropy
    # If only one reads present returns numeric vector rather than dataframe:
    if (class(Viral_reads_contents)=="numeric") {
      Viral_reads_contents = matrix(Viral_reads_contents_mean,ncol = 4)
    }
    
    Viral_reads_contents_mean =colMeans(Viral_reads_contents)
    Read_entropy = sum(-log(Viral_reads_contents_mean)*Viral_reads_contents_mean,na.rm = T)
    #Caclulate the spatial distribution of the mapped reads : how much percent of the genome is mapped ?
    Covered_genome = coverage(BAM_file)[[z]]
    Covered_genome = as.numeric(Covered_genome)
    Spatial_distribution =sum(Covered_genome>0)/length(Covered_genome)
    Covered_genome = rle(sign(Covered_genome))
    Longest_contig = max(Covered_genome$lengths[Covered_genome$values>0])
    
    ##The mean reads quality
    Reads_quality = as.character(BAM_file@elementMetadata$qual)
    Reads_quality = PhredQuality(Reads_quality)
    Reads_quality = as(Reads_quality,"IntegerList")
    Reads_quality = as.numeric(as.matrix(Reads_quality))
    Mean_read_quality = mean(Reads_quality)
    Sd_read_quality = sd(Reads_quality)
    
    ##... the number of mapped reads and unique mapped reads
    N_unique_mapped_reads = sum(BAM_file@elementMetadata$mapq==255) ##Code specific to STAR aligner.... 
    N_mapped_reads = length(BAM_file)
    Percent_uniquely_mapped = N_unique_mapped_reads/N_mapped_reads
    
    ##DUSTy score identifies low-complexity sequences, in a manner inspired by the dust implementation in BLAST
    Mean_dust_score = NA
    Percent_high_quality_reads = NA
    if ("ShortRead"%in%installed.packages()){
      DUST_score = dustyScore(BAM_file@elementMetadata$seq)
      Mean_dust_score = mean(DUST_score)
      Percent_high_quality_reads =  sum(DUST_score<500)/length(DUST_score)
    }
    
    ##Summary Statistics Per Virus 
    QC_temp = c(N_mapped_reads,N_unique_mapped_reads,Percent_uniquely_mapped,
                Mean_read_quality,Sd_read_quality,
                Viral_reads_contents_mean,Read_entropy,Spatial_distribution,Longest_contig,
                Mean_dust_score,Percent_high_quality_reads)
    QC_result = rbind(QC_result, QC_temp)
  }
} else {
  QC_result = data.frame(N_mapped_reads = numeric(),N_unique_mapped_reads=numeric(),Percent_uniquely_mapped=numeric(),
              Mean_read_quality=numeric(),Sd_read_quality=numeric(),
              A = numeric(), C= numeric(), G = numeric(), T = numeric(), Read_entropy=numeric(),Spatial_distribution=numeric(),Longest_contig=numeric(),
              Mean_dust_score=numeric(),Percent_high_quality_reads=numeric())
}
cat("\t Calculating VIRAL QC Metrics.... DONE!. \n", file=log, append = TRUE)


## ------------------------------------------------------------------------------------
## Now we perform filtering based ont the Calculaed QC statsitics: 
## Editied by LEO -- otherwise fails if only one virus is detected as produces numeric vector rather than dataframe. 
if (class(QC_result)=="numeric"){
  cat('Only one virus detected - output is not a dataframe: converting \n', file=log, append = TRUE)
  QC_result <- as.data.frame(t(as.data.frame(QC_result)))
}
colnames(QC_result) = c("N_reads","N_unique_reads","Percent_uniquely_mapped",
                        "Mean_read_quality","Sd_read_quality",
                        c("A","C","G","T"),"Sequence_entropy","Spatial_distribution","Longest_contig",
                        "DUST_score","Percent_high_quality_reads")
QC_result <- as.data.frame(QC_result)
QC_result <- sapply( QC_result, as.numeric )
QC_result <- as.data.frame(QC_result)
if(length(QC_result[,1])>0){
rownames(QC_result) = rownames(temp_chromosome_count)
} 
QC_result = QC_result[QC_result$N_unique_reads>0,]

### ------------------------------------------------------------------------------------
## Now we need to extract information on the mapping by itself to get information about the QC
## This Requires R auxillalry Functions: 
path_to_Log_file = paste0(temp_output_dir, "/", name_target, "Log.final.out")
Mapping_information = Extraction_Log_final(path_to_Log_file)
Mean_mapping_length = Mapping_information$Length_vector[1]
cat("\t 8. Filtering on QC Metrics. \n", file=log, append = TRUE)
detected_virus = rownames(QC_result[QC_result$Sequence_entropy>1.2 & QC_result$Longest_contig>3*Mean_mapping_length & QC_result$Spatial_distribution>0.05,])

## Returning QC Metrics for viruses that passed Filtering 
Filtered_QC=QC_result[detected_virus,]
filteredoutfile <- paste0(temp_output_dir, "/", name_target, "QC_filtered.txt")
unfilteredoutfile <- paste0(temp_output_dir, "/", name_target, "QC_unfiltered.txt")

cat("\t 9. Exporting QC Metrics / identified Viruses to .txt files. \n", file=log, append = TRUE)
##Exporting the tables of the QC analysis
write.table(file = unfilteredoutfile, x = QC_result, quote = F,sep = "\t")
write.table(file = filteredoutfile, x = Filtered_QC, quote = F,sep = "\t")
cat("\t Exporting QC Metrics... DONE! \n", file=log, append = TRUE)


## ------------------------------------------------------------------------------------
# Note removed splice table info as this didnt seem to be used anywhere else?
##------------------------------------------------------------------------------------

## Calculating some additional Statistics needed for Plots: 
## Additional info : % of reads mapped to viral vs host
cat("\t 10. Calculating Broad Metrics on Mapping. \n", file=log, append = TRUE)
Read_count_temp = read.table(temp_chromosome_count_path,header = F,row.names = 1)
colnames(Read_count_temp) = c("Chromosome_length","Mapped_reads","Unknown")
Read_count_temp = Read_count_temp[Read_count_temp$Mapped_reads!=0,]
Read_count_temp$Chr <- rownames(Read_count_temp)

virus_sequence_names <- grep("refseq", rownames(Read_count_temp), value=TRUE)
human_chromosomes <- rownames(Read_count_temp)[rownames(Read_count_temp) %notin% virus_sequence_names]
human_chromosomes <- human_chromosomes[human_chromosomes!="*"]

## Need to rename human chr to append chr to make it easy to grep:
for (x in human_chromosomes){
  Read_count_temp[["Chr"]] <- with(Read_count_temp, ifelse(Chr == x , paste0("chr", Chr), Chr))
}

## Set as Rownames
rownames(Read_count_temp) <- Read_count_temp$Chr
Read_count_temp$Chr <- NULL

## Calculate the percentage of human reads etc 
host_mapping_count = sum(Read_count_temp[grepl(pattern = "chr",rownames(Read_count_temp)),"Mapped_reads"])
viral_mapping_count = sum(Read_count_temp[grepl(pattern = "NC",rownames(Read_count_temp)),"Mapped_reads"])
total_mapping = viral_mapping_count + host_mapping_count
Ratio_host_virus = matrix(data = c(host_mapping_count,viral_mapping_count)/total_mapping,ncol = 1)*100
percent_viral = viral_mapping_count / total_mapping * 100 

## Number of uniquely mapped reads and other reads for the filtered virus 
Mapping_selected_virus = data.frame(Unique_mapping = (Filtered_QC$N_unique_reads),All_mapping = (Filtered_QC$N_reads),row.names = rownames(Filtered_QC))
Mapping_selected_virus = Mapping_selected_virus[order(Mapping_selected_virus$Unique_mapping,decreasing = T),]
cat("\t 10. Calculating Broad Metrics on Mapping...DONE \n", file=log, append = TRUE)
## ------------------------------------------------------------------------------------
## Plotting the Statisitics 

cat("\t 11. Plotting QC plots to pdf. \n", file=log, append = TRUE)
Mapping_selected_virus$Name_id <- rownames(Mapping_selected_virus)

if(length(Mapping_selected_virus[, 1])>0) {
  Mapping_selected_virus$Name_sequence<- c()
  for (i in 1:length(Mapping_selected_virus[, 1])){
    z <- Mapping_selected_virus[i, 3]
    z <- unlist(strsplit(z,"|",fixed = T))[2]
    Mapping_selected_virus$Name_sequence[i] <- z
  } 
} else {
  Mapping_selected_virus <- data.frame(Unique_mapping = numeric(), All_mapping = numeric(), Name_id = character(), Name_sequence = character())
}
Mapping_selected_virus <- merge(Mapping_selected_virus, Virus_database, by="Name_sequence")
rownames(Mapping_selected_virus) <- Mapping_selected_virus$Complete_segment_name

pdf_name <-paste0(temp_output_dir, "/", name_target, "_Summary.pdf")
## Open PDF file: 
pdf(pdf_name, height = 24, width = 20)
par(las=0,mfrow=c(4,3),mar=c(6,6,6,4))
Color_vector = c("lightskyblue1","orange","grey80")

# Plotting the proportion of uniquely mapped reas, unmapped etc...
barplot(Mapping_information$Mapping_result,ylim=c(0,100),xlim=c(0,5),ylab="Percentage of reads (%)",col=Color_vector,cex.lab=1.5,cex.axis = 1.5, main = "Percentage of Reads Mapped")
legend(x = 1.5,y=50,legend = c("Unmapped","Mapped to multiple loci","Uniquely mapped"),bty="n",fill = Color_vector[length(Color_vector):1],cex = 1.5)
# Size of mapping, insertion and deletion
barplot(Mapping_information$Length_vector,col="black",names.arg = c("Mapping Length","Insertion Length","Deletion length"), main = "Size of Mapped Reads",  horiz = TRUE, xlim=c(0,max(Mapping_information$Length_vector[1])*1.2),xlab="Nucleotide length",cex.lab=1.3,cex.axis = 1.3,cex.names=1.3)
#Rate of mismatch, deletion and insertion
barplot(Mapping_information$Rate_vector,col="black",names.arg = c("Mismatch rate","Insertion rate","Deletion rate"), main = "Rate of Mismatching", horiz = TRUE, xlim=c(0,max(Mapping_information$Rate_vector[1])*1.2),xlab="Rate (%)",cex.lab=1.3,cex.axis = 1.3,cex.names=1.3)
# Ratio 
Color_vector = c("darkred","grey")
barplot(Ratio_host_virus,ylim=c(0,100),xlim=c(0,5),ylab="Mapping events (%)",col=Color_vector,cex.lab=1.5,cex.axis = 1.5, main = "Ratio of Host:Viral Mapping")
legend(x = 1.5,y=50,legend = c("Viral mapping","Host mapping"),bty="n",fill = Color_vector[length(Color_vector):1],cex = 1.5)

##Get colour strings
if (length(detected_virus) > 0) {
  Color_vector= string.to.colors(factor(rownames(QC_result)%in%detected_virus),colors = c("orange","green")) ##Viral sequences that passed QC : green
}

if (length(QC_result[,1]) > 0 & length(detected_virus)==0) {
    Color_vector= string.to.colors(factor(rownames(QC_result)%in%detected_virus),colors = c("orange")) ##Viral sequences that passed QC : green
}

# Plot Number of Reads > 50 (filtering threshold)
if(length(rownames(QC_result))>0){
  plot(QC_result$N_reads,QC_result$Spatial_distribution*100,pch=21,bg=Color_vector,log="x",cex=1.5,xlab="Number of Mapped Reads",ylab="% Mapped genome",ylim=c(0,100),cex.lab=1.5,main="Viral Summary Multimapping Reads: Read Count")
  abline(h=5,lwd=2,lty=2,col="grey")
  abline(v=Minimal_read_mapped,lwd=2,lty=2,col="grey")

  # Plot Number of Unique Reads > 50 (filtering threshold)
  plot(QC_result$N_unique_reads,QC_result$Spatial_distribution*100,pch=21,bg=Color_vector,log="x",cex=1.5,xlab="Number of Uniquely Mapped Reads",ylab="% Mapped genome",ylim=c(0,100),cex.lab=1.5,main="Viral Summary Multimapping Reads: Unique Read Count")
  abline(h=5,lwd=2,lty=2,col="grey")

  #Second QC for the viral hits 
  plot(QC_result$Sequence_entropy,QC_result$Spatial_distribution*100,pch=21,bg=Color_vector, cex=1.5,xlab="Sequence Complexity",ylab="% Mapped genome",cex.lab=1.5,ylim=c(0,100),main="Viral Summary Multimapping Reads: Sequence Complexity")
  abline(h=5,lwd=2,lty=2,col="grey")
  abline(v=1.2,lwd=2,lty=2,col="grey")

  #Third QC for the viral hits 
  plot(QC_result$Longest_contig,QC_result$DUST_score,pch=21,bg=Color_vector, cex=1.5,xlab="Longest Contig (nt)",ylab="DUST Score",cex.lab=1.4,main="Viral Summary Multimapping Reads: DUST Score")
  abline(v=3*Mean_mapping_length,lwd=2,lty=2,col="grey")

#Number of reads for each filtered virus
  if (length(detected_virus) > 0) {
    par(las=0)
    barplot(Mapping_selected_virus$All_mapping[nrow(Mapping_selected_virus):1],
            col="black",horiz = T,cex.lab=1.3,cex.axis = 1.3,cex.names=1,
            xlim=c(0,max(Mapping_selected_virus$All_mapping)*1.2),xlab="Number of mapped reads",
            names.arg = rownames(Mapping_selected_virus)[nrow(Mapping_selected_virus):1])
    par(las=0)
    barplot(Mapping_selected_virus$Unique_mapping[nrow(Mapping_selected_virus):1],
            col="black",horiz = T,cex.lab=1.3,cex.axis = 1.3,cex.names=1,
            xlim=c(0,max(Mapping_selected_virus$Unique_mapping)*1.2),xlab="Number of uniquely mapped reads",
            names.arg = rownames(Mapping_selected_virus)[nrow(Mapping_selected_virus):1])
    
  }
}
dev.off() 

cat("\t 11. Plotting QC plots Done. \n", file=log, append = TRUE)
cat("\t ----------------------------------------------\n", file=log, append=TRUE)
cat("QC PDF and Summary Statistics Done and Saved.  \n", file=log, append = TRUE)
cat("----------------------------------------------\n", file=log, append=TRUE)


## ------------------------------------------------------------------------------------

cat("Merging Viral SAM files for Identified Viruses into one Viral Bam File. \n", file=log, append = TRUE)
list_BAM_files = paste0(temp_output_dir,"/Viral_BAM_files/")
selected_virus = list.files(list_BAM_files,full.names=F)
if (length(rownames(Filtered_QC)) >= 1){
  selected_virus = base::strsplit(x = selected_virus,split = ".bam")
  selected_virus = unlist(lapply(selected_virus, function(x) {x[1]}))
  
  list_BAM_files = list.files(list_BAM_files,full.names=T)
  names(list_BAM_files) = selected_virus
  
  list_BAM_files = list_BAM_files[rownames(Filtered_QC)]
  list_BAM_files = paste("\'", list_BAM_files, "\'", sep = "")
  
  Merging_BAM_commad = paste("samtools merge",paste(temp_output_dir,"/Merged_viral_mapping.bam",sep = ""),list_BAM_files)
  system(Merging_BAM_commad)
  cat("VIRUS Present post filtering: viral BAM made. \n", file=log, append = TRUE)
} else {
  cat("No VIRUS Present post filtering - no viral BAM made. \n", file=log, append = TRUE)
}


## ------------------------------------------------------------------------------------
## END of Log Script: 
cat("----------------------------------------------\n", file=log, append=TRUE)
cat("ViralTrack_Scanning.2.0 COMPLETE \n", file=log, append=TRUE)
end_time <- Sys.time()
cat(paste0("End time: ", end_time, "\n"), file=log, append=TRUE)
Duration <- start_time_1 - end_time
cat(paste0("End time: ", Duration, "\n"), file=log, append=TRUE)
cat("----------------------------------------------\n", file=log, append=TRUE)


## ------------------------------------------------------------------------------------
cat("----------------------------------------------\n", file=log, append=TRUE)
cat("ViralTrack_Demultiplexing 2.0  \n", file=log, append=TRUE)
start_time <- Sys.time()
cat("Start Time: ", start_time, "\n", file=log, append=TRUE)
cat("----------------------------------------------\n", file=log, append=TRUE)

## -----------------------------------------------------------------------------------
cat("Aggregating All Detected Viral and Human BAM files into Reads to Demultipled. \n", file=log, append=TRUE)

## First aggregating all the reads from the detected viruses:
List_bam_files =c()
## Identify which viral bam files we want to demultiplex (those for filtered viruses)
Identified_viral_fragments <- detected_virus  
## Assuming viral reads are present: 
for (segment_temp in Identified_viral_fragments) {
  List_bam_files = c(List_bam_files, paste(temp_output_dir,"/Viral_BAM_files/",segment_temp,".bam",sep = ""))
}
## Also we want to demultiplex human Reads
path_to_human <- paste0(temp_output_dir, "/HUMAN_BAM_files")
List_bam_files <- c(List_bam_files, list.files(path_to_human, recursive = TRUE, full.name=TRUE))
List_bam_files = paste("\'",List_bam_files,"\'",sep = "")
  
## Merge these bam files into a Read to demultiplex file:   
command_merge = base::paste("samtools merge ", temp_output_dir,"/Reads_to_demultiplex.bam -f ",paste(List_bam_files,collapse = " "),sep="")
system(command_merge)

cat("Aggregation Complete. \n", file=log, append=TRUE)
cat("----------------------------------------------\n", file=log, append=TRUE)
cat("Performing Feature Counts on Reads To Demultiplex.  \n", file=log, append=TRUE)
start_time_f <- Sys.time()
cat(paste0("-T: Threads for featureCounts: ", N_thread, " \n"), file=log, append=TRUE)
cat("Start Time Feature Counts: ", start_time_f, "\n", file=log, append=TRUE)
## Assigning reads to transcripts using Rsubread Featurecounts
featurecommand = paste("featureCounts -T ", N_thread, " -M --primary -t transcript -R BAM -g gene_id -a ", path_to_gtf, " -o ", temp_output_dir, "/counts.txt ", temp_output_dir, "/Reads_to_demultiplex.bam 2> ", name_prefix, "_FEATURE_COUNTS.log", sep="")
system(featurecommand)
end_time_f <- Sys.time()
cat("End Time Feature Counts: ", end_time_f, "\n", file=log, append=TRUE)
cat("Feature Counts Complete.  \n", file=log, append=TRUE)
cat("----------------------------------------------\n", file=log, append=TRUE)
  
## We now have to sort the BAM file:
cat("Sorting Bam File.  \n", file=log, append=TRUE)
command_sort =paste("samtools sort ",temp_output_dir,"/Reads_to_demultiplex.bam.featureCounts.bam -o ",temp_output_dir,"/Assigned_sorted.bam",sep = "")
system(command_sort)
cat("Done.  \n", file=log, append=TRUE)

## We now index the BAM file:
cat("Indexing Bam File.  \n", file=log, append=TRUE)
command_index =paste("samtools index ",temp_output_dir,"/Assigned_sorted.bam",sep = "")
system(command_index)
cat("Done.  \n", file=log, append=TRUE)
  
## Final command : Umi-tools command
## This will deduplicate reads 
cat("Starting UMI-tools Count.  \n", file=log, append=TRUE)
command_umi_tools = paste("umi_tools count --per-gene --gene-tag=XT --assigned-status-tag=XS --per-cell -I ",
                          temp_output_dir,"/Assigned_sorted.bam  -S ",temp_output_dir, "/Expression_table.tsv --wide-format-cell-counts > ", name_prefix, "_UMITOOLS_COUNTS.log", sep="")
system(command_umi_tools)
cat("Done.  \n", file=log, append=TRUE)
## -----------------------------------------------------------------------------------
## Converting to Matrix Format: 
cat("----------------------------------------------\n", file=log, append=TRUE)

cat("Reading Expression tsv file.  \n", file=log, append=TRUE)
# Saving matrix as Sparse Matrixform in an outs directory:
MTX_dir = paste(temp_output_dir,"/viral_filtered_features/", sep = "")
dir.create(MTX_dir)
  
# Read non-sparse features:
file_tsv = paste0(temp_output_dir, "/Expression_table.tsv")
expression_table <- read.table(file = file_tsv, sep = '\t', header = TRUE, row.names = 1)

# Notin function:
`%notin%` <- Negate(`%in%`)

# Convert to data-frame to allow for Calculation of Summary Statistics
expression_table <- as.data.frame(t(expression_table))
viral_cols <- grep("refseq", colnames(expression_table), value=TRUE)
human_cols <- colnames(expression_table)[!colnames(expression_table) %in% viral_cols]
mito_cols <- grep("MT_1", colnames(expression_table), value=TRUE)

cat("Calculating Meta Data Percentages.  \n", file=log, append=TRUE)
# Calculating meta data percentages:
# Percent human reads per cell:
cat("\t Human Read Percentage/counts Calculated.  \n", file=log, append=TRUE)
expression_table$human_chromosome_percent <- (rowSums(expression_table[, c(human_cols)]))/ (rowSums(expression_table))* 100
expression_table$human_chromosome_counts <- rowSums(expression_table[, c(human_cols)])
# Percent mito reads per cell:
cat("\t Human Mitochondrial Read Percentage Calculated.  \n", file=log, append=TRUE)
expression_table$human_mitochondrial_percent <- (expression_table[, c(mito_cols)])/ (rowSums(expression_table))* 100
expression_table$human_mitochondrial_counts <- expression_table[, c(mito_cols)]

# Percent viral read per cell - for all viruses:

if(length(viral_cols)>= 1){
  cat("VIRAL reads PRESENT. Calculating Total and Individual Viral Percentages. \n", file=log, append=TRUE)
  if(length(viral_cols)==1){
    cat("Only One Virus Present. \n", file=log, append=TRUE)
    expression_table$percent_viral_reads <- expression_table[, c(viral_cols)]/ (rowSums(expression_table[, c(viral_cols, human_cols)]))* 100
  } else {
    cat("Multiple Viruses Present. \n", file=log, append=TRUE)
    expression_table$percent_viral_reads <- (rowSums(expression_table[, c(viral_cols)]))/ (rowSums(expression_table[, c(viral_cols, human_cols)]))* 100  
  } 
} else {
  cat("No Viral Reads Present: percent_viral_reads = 0. \n", file=log, append=TRUE) 
  expression_table$percent_viral_reads <- "0"
}

# Percent viral reads per cell per virus.
if(length(viral_cols)>= 1){
  cat("Calculating Percentage per Virus. \n", file=log, append=TRUE)
  for (i in viral_cols){
    y <- strsplit(i, split='|', fixed=TRUE)
    y = unlist(y)
    name=paste0("percent_viral_reads_", y[2], y[4])
    virus <- i
    expression_table[, name] <- (expression_table[, c(virus)])/ (rowSums(expression_table[, c(viral_cols, human_cols)]))* 100
  }
} 
cat("Finished Calculating Viral Percentages. \n", file=log, append=TRUE)

## Subsetting Matrix to remove individualised human counts and leave a Summary:
expression_table <- expression_table[, colnames(expression_table) %notin% human_cols]

## All percentages calculated
## Converting Table To SPARSE MATRIX FORMAT
cat("Converting Expression Table into Sparse Format. \n", file=log, append=TRUE)
expression_table <- (t(expression_table))
expression_table <- Matrix(expression_table, sparse = TRUE)
cat("Done. \n", file=log, append=TRUE)
cat("Writing Viral_counts.mtx. \n", file=log, append=TRUE)
writeMM(obj = expression_table, file=paste0(MTX_dir, "viral_counts.mtx"))
cat("Writing barcodes.tsv. \n", file=log, append=TRUE)
cat("Writing genomes.tsv. \n", file=log, append=TRUE)
cellbarcode <- colnames(expression_table)
genome <- rownames(expression_table)
write.table(cellbarcode, file=paste0(MTX_dir, 'barcodes.tsv'), quote=FALSE, sep='\t', col.names = FALSE, row.names = FALSE)
write.table(genome, file=paste0(MTX_dir, 'genomes.tsv'), quote=FALSE, sep='\t', col.names = FALSE, row.names = FALSE)
cat("Done. \n", file=log, append=TRUE)


## Now we just move files of interest into the Report directory: 
cat("Creating Report Directory. \n", file=log, append=TRUE)
Report_dir <- paste0(MTX_dir, "Reports")
dir.create(Report_dir)
cat("Moving Useful Files into Report Directory. \n", file=log, append=TRUE)
move <- paste0("cp ", temp_output_dir, "/counts.txt.summary ", Report_dir )
system(move)
move <- paste0("cp ", pdf_name, " ", Report_dir )
system(move)
move <- paste0("cp ", filteredoutfile, " ", Report_dir )
system(move)
move <- paste0("cp ", unfilteredoutfile, " ", Report_dir )
system(move)
cat("Done. \n", file=log, append=TRUE)

## DONE!
cat("----------------------------------------------\n", file=log, append=TRUE)
cat("DEMULTIPLEXING COMPLETE. \n", file=log, append=TRUE)
end_time <- Sys.time()
cat("Demultipexing End Time ", end_time, "\n", file=log, append=TRUE)
Time_difference <- end_time - start_time
cat("Run Time ", Time_difference, "\n", file=log, append=TRUE)
cat("----------------------------------------------\n", file=log, append=TRUE)
cat("----------------------------------------------\n", file=log, append=TRUE)
cat("PIPELINE COMPLETE. \n", file=log, append=TRUE)
cat("----------------------------------------------\n", file=log, append=TRUE)