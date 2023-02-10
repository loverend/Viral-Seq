# Introduction 
Viral Detection in bulk RNA sequencing data (paired or unpaired). 

# Author Contribution
 - Lauren Overend (lauren.overend@oriel.ox.ac.uk) 
 - Refactored from Viral-Track for scRNAseq (Pierre Bost https://github.com/PierreBSC/Viral-Track)

# Job Submission: 
- For paired end reads:
```
qsub -t 1:2 run_viral_track_rnaseq.sh SampleSheets/ViralSeq_SampleSheet_RepertoireSamples_2023.txt /gpfs2/well/immune-rep/shared/MISEQ/LEO_GAinS_RNASEQ_2023/ViralSeq TRUE
```

# Arguments: 
- ```-t 1:n```
- number of samples from sample sheet, e.g you could run just sample 1 with qsub -1 t or you could run 1-100 with qsub -t 1:100
- ```ViralSeq_SampleSheet_RepertoireSamples_2023.txt```
- sample sheet (see sample folder for examples). 
- ```/gpfs2/well/immune-rep/shared/MISEQ/LEO_GAinS_RNASEQ_2023/ViralSeq```
- outputdirectory 
- ```TRUE/FALSE```
- = if paired end use TRUE 

# Changing job duration: 
- Note you may need to edit the number of nodes, and long or short queue in this script: run_viral_track_rnaseq.sh
- ```-pe shmem 1 ```
- = 1 node (1/2 is good for unmapped reads). 
- For mapping all reads you may need up to 16 nodes. 
- ```-q long.qc ```
- for all reads, short.qc for unmapped reads. 

# Output
- During analaysis a log file will be output in the outputdirectory detailing the run paramaters etc 
- Output will include:
  -  a counts matrix of human and viral counts
  -  a pdf with viral QC measures including longest contig, entropy and genome coverage 
