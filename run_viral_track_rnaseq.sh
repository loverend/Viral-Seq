#!/bin/bash
#$ -cwd
#$ -N ViralTrackSepsis
#$ -q short.qc
#$ -pe shmem 2
#$ -e COMMAND_LOGS/
#$ -o COMMAND_LOGS/

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
OUTPUTDIR=$2
PAIRED=$3
# Task Arguments
SAMPLE=$(awk -F '\t' "{if (NR==$SGE_TASK_ID) print \$1}" $SAMPLES_FILE)      # Sample ID

## Check pairing paramater
if [[ "${PAIRED}" -ne "TRUE" && "${PAIRED}" -ne "FALSE" ]]; then
echo "ERROR: UNKNOWN PAIRED PARAMETER PROVIDED..."
exit 888
fi 



echo "********************************************************"
echo "* Job Details"
echo "********************************************************"
echo "SGE job ID       : "$JOB_ID
echo "SGE task ID      : "$SGE_TASK_ID
echo "Run on host      : "`hostname`
echo "Operating system : "`uname -s`
echo "Username         : "`whoami`
echo "Started at       : "`date`
echo
echo "********************************************************"
echo "* Job Parameters"
echo "********************************************************"
echo "SAMPLE FILE     : ${SAMPLES_FILE}"
echo "SAMPLE          : ${SAMPLE}"
echo "OUTPUTDIR       : ${OUTPUTDIR}"
echo "PAIRED       	  : ${PAIRED}"
echo "********************************************************"

# RunJob
CMD="Rscript /well/immune-rep/shared/CODE/Viral-Seq/ViRNAseq_Unique.R -f ${SAMPLE} -n 2 -o ${OUTPUTDIR} -p ${PAIRED}"

echo "********************************************************"
echo "["`date`"] Running TCR/BCR"
echo "********************************************************"
echo "Command : ${CMD}"
echo "********************************************************"
echo
eval "${CMD}"

## If Job runs successfully sample will be output to a check file
NWCMD="echo ${SAMPLE} >> job_${1}_SuccessfullCompletedSamples.txt"
eval "${NWCMD}"

# Done 
echo
echo "********************************************************"
echo "["`date`"] Done"
echo "********************************************************"
exit 0
