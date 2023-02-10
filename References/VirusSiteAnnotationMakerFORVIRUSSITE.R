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
  make_option(c("-g", "--genomes"), action="store", type="character", default = "/gpfs2/well/immune-rep/users/kvi236/References/VIRUS_REFERENCE/genomes.fasta", help="Path to genomes .fa Downloaded from VirusSite [default]"),
  make_option(c("-o", "--output"), action="store", type="character", default="/gpfs2/well/immune-rep/users/kvi236", help="Path to output directory [default]"),
  make_option(c("-e", "--extra_files"), action="store", type="character", default="/well/immune-rep/users/kvi236/References/VIRUS_REFERENCE/covid-19.fasta", help="Comma seperated list of extra fasta files")

)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, print_help_and_exit = TRUE, args = commandArgs(trailingOnly = TRUE) )




## Read in the VirusSite genomes.fa file:
l <- read.fasta(opt$genomes, seqtype = "DNA", as.string = TRUE, forceDNAtolower = TRUE,set.attributes = TRUE)
annotation <- unlist(getAnnot(l))
all_genomes <- names(l)
## Make an Empty Data Frame
df <- NULL
## Fill with Relevant Columns
for (i in 1:length(all_genomes)){
  z <- all_genomes[i]
  z <- unlist(strsplit(z,"|",fixed = T))
  nc <- z[2]
  nucleo <- unlist(strsplit(z[3],split='nt', fixed=TRUE))
  x <- annotation[i]
  f <- unlist(strsplit(x, "|", fixed=TRUE))[4]
  q <- unlist(strsplit(f, ",", fixed=TRUE))[1]
  row <- c(nc, nucleo, q, f)
  df = rbind(df, row)
}
colnames(df) <- c("Name_sequence", "Genome_length", "Virus_name", "Complete_segment_name")

## Extra fasta files:

if(length(opt$extra_files)>0){
  extra_files <- opt$extra_files
  extra_files <- unlist(strsplit(extra_files, ","))
}

if (length(extra_files)>=1){
  extra_df <- NULL
  for (i in extra_files){
    f <- read.fasta(i, seqtype = "DNA", as.string = TRUE, forceDNAtolower = TRUE,set.attributes = TRUE)
    annotation <- unlist(getAnnot(f))
    all_genomes <- names(f)
    for (i in 1:length(all_genomes)){
      z <- all_genomes[i]
      z <- unlist(strsplit(z,"|",fixed = T))
      nc <- z[2]
      nucleo <- unlist(strsplit(z[3],split='nt', fixed=TRUE))
      x <- annotation[i]
      f <- unlist(strsplit(x, "|", fixed=TRUE))[4]
      q <- unlist(strsplit(f, ",", fixed=TRUE))[1]
      row <- c(nc, nucleo, q, f)
      extra_df = rbind(extra_df, row)
    }
  }
} else {
  extra_df <- data.frame(Name_sequence = character(), Genome_length=character(), Virus_name=character(), Complete_segment_name=character())
}
colnames(extra_df) <- c("Name_sequence", "Genome_length", "Virus_name", "Complete_segment_name")
## Make Final Df by merging virus site and Extra Files 
final_df <- rbind(df, extra_df)

## Save
write.table(final_df, file = paste0(opt$output, "/Updated_VirusSite_Reference.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote=FALSE)
