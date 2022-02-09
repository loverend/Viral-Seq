## ---------------------------
## Script name: VirusSiteAnnotationMaker
## Function: Generate NCBI Annotation file for updated NCBI viral references.fa file
##           Compatible with Viral-Seq
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
  make_option(c("-o", "--output"), action="store", type="character", default="/well/immune-rep/shared/CODE/VIRAL_SEQ_Reference/", help="Path to output directory [default]"),
  make_option(c("-g", "--genomes"), action="store", type="character", default="/well/immune-rep/shared/CODE/VIRAL_SEQ_Reference/Viral_FA_REFERENCE", help="List to directory containing extra filepaths")

)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, print_help_and_exit = TRUE, args = commandArgs(trailingOnly = TRUE) )

extra_files <- list.files(opt$g, full.names=TRUE)
extra_files <- grep("fasta", extra_files, value=TRUE)

## Read in the NCBI reference files:
references <- list()
for(i in 1:length(extra_files)){
	l <- read.fasta(extra_files[i], seqtype = "DNA", as.string = TRUE, forceDNAtolower = TRUE,set.attributes = TRUE)
	print(length(l))
	references <- append(l, references)
}	

## names of sequences
annotation <- unlist(getAnnot(references))
## Lengths of sequences 
length <- unlist(getLength(references))


## Make an Empty Data Frame
df <- NULL
## Fill with Relevant Columns
## Must extract NCBI, name and full name
for (i in 1:length(references)){
  z <- unlist(getAnnot(references)[i])
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

## Bind to genome lengths
df <- cbind(length, df)
df <- df[, c(2, 1, 3,4 )]
colnames(df) <- c("Name_sequence", "Genome_length", "Virus_name", "Complete_segment_name")

## Save
write.table(df, file = paste0(opt$output, "/NCBI_Viral_Seq_Reference.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote=FALSE)
print("Done")

########################################################################