#!/bin/sh
#!/bin/bash
#SBATCH --mail-user=leven.terziev@dkf.unibe.ch
#SBATCH --mail-type=end
#SBATCH --job-name="Chip61"
#SBATCH --nodes=2
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=64G
#SBATCH --partition=all

##SBATCH --output=/path/to/outfile
##SBATCH --error=/path/to/errfile

# For array jobs
# Array job containing 100 tasks, run max 10 tasks at the same time
##SBATCH --array=1-100%10


#### Your shell commands below this line ####

## QIIME pipeline shell script ##

## Load QIIME module Version 1.8.0 ##
module load vital-it
module load UHTS/Analysis/qiime/1.8.0

## Python 2.7 should now be used automatically. UBELIX default is 2.6.6. R is required for p value calculation ##
## module load python/2.7 (deleted after .bashcr file was changed to inlcude 'export PATH=/software/bin:$PATH')
module load R/latest

## Use my custom config in order to have qiime use $TMPDIR instead of /tmp on the compute nodes
export QIIME_CONFIG_FP=~/.qiime_config

## Copy required files (Iontorrent fastq and mapping TXT) from HOME directory to local scratch on node ($TMPDIR) ##
cp Chip61.fastq $TMPDIR/file.fastq
cp Chip_61_MappingFile_3000.txt $TMPDIR/map.txt
mkdir ./Results_$SLURM_JOB_ID/

## Start of QIIME pipeline ##
## Convert fastq file into .fna and .qual files ##
convert_fastaqual_fastq.py -c fastq_to_fastaqual -f $TMPDIR/file.fastq -o $TMPDIR/fastaqual
cp -R $TMPDIR/fastaqual ./Results_$SLURM_JOB_ID

## Split libraries based on barcodes provided in map.txt ##
split_libraries.py -m $TMPDIR/map.txt -f $TMPDIR/fastaqual/file.fna -q $TMPDIR/fastaqual/file.qual -o $TMPDIR/split_output/ -b variable_length
cp -R $TMPDIR/split_output ./Results_$SLURM_JOB_ID
