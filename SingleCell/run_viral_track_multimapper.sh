#!/bin/bash

#$ -cwd
#$ -N ViralTrackMM
#$ -q short.qc
#$ -pe shmem 10

# Load software modules
module purge
module use -a /apps/eb/dev/ivybridge/modules/all
module load Python/3.8.2-GCCcore-9.3.0 
module load UMI-tools/1.0.1-foss-2020a-Python-3.8.2 
module load R-bundle-Bioconductor/3.11-foss-2020a-R-4.0.0
module load SAMtools/1.10-GCC-9.3.0
module load STAR/2.7.3a-GCC-9.3.0
module load Subread/2.0.1-GCC-9.3.0

# Job Arguments
SAMPLES_FILE=$1

# Task Arguments
SAMPLE=$(awk -F '\t' "{if (NR==$SGE_TASK_ID) print \$1}" $SAMPLES_FILE)      # Sample ID

# RunJob
CMD="Rscript /well/immune-rep/users/kvi236/ViralTrackProgram/RPipeline/Viral-Track/ViralTrack2.0/Viral_Track_Multimapping.R -f ${SAMPLE} -n 10"
eval "${CMD}"


# Done 
echo
echo "********************************************************"
echo "["`date`"] Done"
echo "********************************************************"
exit 0
