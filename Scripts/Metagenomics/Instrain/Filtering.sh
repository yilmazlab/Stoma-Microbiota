#!/bin/sh
#!/bin/bash
#SBATCH --mail-user=bahtiyar.yilmaz@dbmr.unibe.ch
#SBATCH --mail-type=end
#SBATCH --job-name="BIN_CLASSIFICATION"
#SBATCH --nodes=1
#SBATCH --time=96:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=24
#SBATCH --partition=epyc2

##SBATCH --output=~/StomaMetagenomics/SlurmFiles
##SBATCH --error=~/StomaMetagenomics/SlurmFiles

#For array jobs
#Array job containing 10 tasks, run max 20 tasks at the same time
#SBATCH --array=1-43%20


#### Your shell commands below this line ####
#### Tutorial: https://github.com/bxlab/metaWRAP/blob/master/Usage_tutorial.md

export LC_CTYPE=en_US.UTF-8
export LC_ALL=en_US.UTF-8

source activate metawrap-env
module load vital-it

#### Your shell commands below this line ####
#### Tutorial: https://github.com/bxlab/metaWRAP/blob/master/Usage_tutorial.md

# Folders
cd /storage/homefs/terziev/StomaMetagenomics

gunzip *.gz

mkdir RAW_READS
mv *fastq RAW_READS
 ls RAW_READS

# Generating the ID list for "for loop"
cd RAW_READS
ls -1 *R1*.gz | awk -F '_' '{print $1}' | sort | uniq > ID # Create IDs
cd ..

# Installation of  bmtagger hg38
# https://github.com/bxlab/metaWRAP/blob/master/installation/database_installation.md
# Step 1: Run metaWRAP-Read_qc to trim the reads and remove human contamination
# Keep it in your home directory
# The main things this pipeline accomplishes are read trimming based on quality scores, and removal of human sequences.     
mkdir READ_QC
metawrap read_qc -1 RAW_READS/SF${SLURM_ARRAY_TASK_ID}_R1.fastq.gz -2 RAW_READS/SF${SLURM_ARRAY_TASK_ID}_R2.fastq.gz -t 24 -o READ_QC/SF${SLURM_ARRAY_TASK_ID}

mkdir CLEAN_READS
mv READ_QC/SF${SLURM_ARRAY_TASK_ID}/SF${SLURM_ARRAY_TASK_ID}_R1_val_1.fq.gz CLEAN_READS/SF${SLURM_ARRAY_TASK_ID}_1.fastq.gz
mv READ_QC/SF${SLURM_ARRAY_TASK_ID}/SF${SLURM_ARRAY_TASK_ID}_R2_val_2.fq.gz CLEAN_READS/SF${SLURM_ARRAY_TASK_ID}_2.fastq.gz