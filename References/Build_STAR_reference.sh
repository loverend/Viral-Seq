#!/bin/bash
#$ -cwd
#$ -N annotate 
#$ -q long.qc
#$ -pe shmem 8

module use -a /apps/eb/dev/ivybridge/modules/all
module load Python/3.8.2-GCCcore-9.3.0 
module load UMI-tools/1.0.1-foss-2020a-Python-3.8.2 
module load R-bundle-Bioconductor/3.11-foss-2020a-R-4.0.0
module load SAMtools/1.10-GCC-9.3.0
module load STAR/2.7.3a-GCC-9.3.0


gunzip HUMAN_GENOME/*
STAR --runThreadN 7 --runMode genomeGenerate --genomeDir VIRAL_SEQ_273a_NCBI/ --genomeFastaFiles HUMAN_Reference/*.fa Viral_FA_REFERENCE/*.fasta