## ---------------------------
## Script name: VirusSiteAnnotationMaker
## Function: Generate VirusSite Annotation file for updated virus_site genomes.fa file
##           Compatible with ViralTrack
## Author: Lauren Overend (LEO).
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
library(seqinr)
suppressMessages(library(optparse))
parser <- OptionParser()
option_list <- list(
  make_option(c("-g", "--genomes"), action="store", type="character", default = "/gpfs2/well/immune-rep/users/kvi236/References/VIRUS_REFERENCE/NCBI/ncbi_reference.fasta", help="Path to genomes .fa Downloaded from VirusSite [default]"),
  make_option(c("-o", "--output"), action="store", type="character", default="/gpfs2/well/immune-rep/users/kvi236/References/", help="Path to output directory [default]"),
  make_option(c("-e", "--extra_files"), action="store", type="character", default="/gpfs2/well/immune-rep/users/kvi236/References/VIRUS_REFERENCE/Subset_Reference", help="List to directory containing extra filepaths")

)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, print_help_and_exit = TRUE, args = commandArgs(trailingOnly = TRUE) )

extra_files <- list.files(opt$e, full.names=TRUE)
extra_files <- grep("fasta", extra_files, value=TRUE)
## Read in the VirusSite genomes.fa file:
l <- read.fasta(opt$genomes, seqtype = "DNA", as.string = TRUE, forceDNAtolower = TRUE,set.attributes = TRUE)
annotation <- unlist(getAnnot(l))
## There are a couple of viruses I want to remove from the reference. 
## ebv type 1 "NC_009334.1"
# macine herpesvirus "NC_004812.1"
# hep e rat "NC_038504.1"
# heron hep c "NC_001486.1"


influenza_a <- grep("Influenza A", annotation, value = TRUE)
influenza_a <- unlist(strsplit(influenza_a, " "))
influenza_a <- grep(">NC", influenza_a, value = TRUE)
influenza_a <- unlist(strsplit(influenza_a, ">"))
influenza_a <- grep("NC", influenza_a, value = TRUE)

influenza_b <- grep("Influenza B", annotation, value = TRUE)
influenza_b <- unlist(strsplit(influenza_b, " "))
influenza_b <- grep(">NC", influenza_b, value = TRUE)
influenza_b <- unlist(strsplit(influenza_b, ">"))
influenza_b <- grep("NC", influenza_b, value = TRUE)

influenza_c <- grep("Influenza C", annotation, value = TRUE)
influenza_c <- unlist(strsplit(influenza_c, " "))
influenza_c <- grep(">NC", influenza_c, value = TRUE)
influenza_c <- unlist(strsplit(influenza_c, ">"))
influenza_c <- grep("NC", influenza_c, value = TRUE)

viruses_remove <- c(influenza_a, influenza_b, influenza_c, "NC_009334.1", "NC_004812.1", "NC_038504.1", "NC_001486.1")
`%notin%` <- Negate(`%in%`)

## Removed dodgy viruses which are causing issues!!!!
l_subset <- l[names(l) %notin% viruses_remove]
l <- l_subset
lengths <- getLength(l)
nc <- getName(l)


write.fasta(sequences=l_subset,names=names(l_subset), file.out="/gpfs2/well/immune-rep/users/kvi236/References/VIRUS_REFERENCE/NCBI_subset.fa", nbchar = 60)




## Make an Empty Data Frame
df <- NULL
## Fill with Relevant Columns
for (i in 1:length(l)){
  z <- unlist(getAnnot(l)[i])
  z <- unlist(strsplit(z,">",fixed = T))
  z <- z[2]
  z <- unlist(strsplit(z," ",fixed = T))
  nc <- z[1]
  z <- z[-1]
  z <- paste(z, collapse = ' ')
  q <- unlist(strsplit(z, ",", fixed=TRUE))[1]
  row <- c(nc, q, z)
  df = rbind(df, row)
}

df <- cbind(lengths, df)
df <- df[, c(2, 1, 3,4 )]
colnames(df) <- c("Name_sequence", "Genome_length", "Virus_name", "Complete_segment_name")
## Extra fasta files:


if (length(extra_files)>=1){
   df_2 <- NULL
   for (i in 1:length(extra_files)){
	path <- extra_files[i]
    u <- read.fasta(path, seqtype = "DNA", as.string = TRUE, forceDNAtolower = TRUE,set.attributes = TRUE)
	if(length(u)>1){
	for (i in 1:length(u)){
		  z <- unlist(getAnnot(u)[i])
		  length <- unlist(getLength(u)[i])
		  z <- unlist(strsplit(z," |",fixed = T))
		  name <- z[2]
		  z <- unlist(strsplit(z[1],">",fixed = T))
		  nc <- z[2]
		  row <- c(nc, length, name, name)
		  df_2 = rbind(df_2, row)
	} 
	}
	if(length(u)==1){
		annotation <- unlist(getAnnot(u))
		length <- unlist(getLength(u))
		z  <- names(u)
		z <- unlist(strsplit(z,"|",fixed = T))
		nc <- grep("NC", z, value=TRUE)
		
		x <- annotation
		f <- unlist(strsplit(x, "|", fixed=TRUE))
		f <- f[length(f)]
		q <- f
		row <- c(nc, length, q, f)
		df_2 = rbind(df_2, row)
    }
	}
	}
	
colnames(df_2) <- c("Name_sequence", "Genome_length", "Virus_name", "Complete_segment_name")
## Make Final Df by merging virus site and Extra Files 
final_df <- rbind(df, df_2)

## Save
write.table(final_df, file = paste0(opt$output, "/Updated_VirusSite_Reference_subset.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote=FALSE)
