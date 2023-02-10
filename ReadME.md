# Viral-Seq Single Cell / Bulk RNA seq Analysis

- Refactored by Lauren Overend
- Adapted from Viral-Track (scRNAseq) availible from: https://www.nature.com/articles/s41597-019-0116-4, Author Pierre Bost)
- Performs human and viral mapping of either bulk RNA seq (paired/unpaired) or single cell RNAseq data. 
-------------------------------------------------


## Dependencies (availible on Rescomp using the following modules)

```
module load UMI-tools/1.0.1-foss-2020a-Python-3.8.2 
module load R-bundle-Bioconductor/3.11-foss-2020a-R-4.0.0
module load SAMtools/1.10-GCC-9.3.0
module load STAR/2.7.3a-GCC-9.3.0
module load Subread/2.0.1-GCC-9.3.0
```
## R dependencies 
```
library(stringr)
library(Biostrings)
library(ShortRead)
library(doParallel)
library(GenomicAlignments)
library(Matrix)
library(ggplot2)
library(cowplot)
library(reshape2)
library(optparse)
```

# Build A Viral Reference 
- Summary: The first step consists in creating a **STAR** index that include both host and virus reference genomes.
- Viral Refererence:  
- Host genome has be downloaded from the [ensembl website](https://www.ensembl.org/info/data/ftp/index.html). 
- We include a reference downloaded from NCBI Virus filtering (as of Feb 2022, https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/), which has been filtered to include those relevent to a UK population. Search terms were human host and reference sequence. We also removed EBV type 2 and kept the gene which discriminates the two types to prevent multi-mapping issues.
- Edit Build_STAR_reference.sh with locations of files including the host/viral files 
- Run: 
```./Build_STAR_reference.sh```
