#!/bin/sh
#!/bin/bash
#SBATCH --mail-user=leven.terziev@dbmr.unibe.ch
#SBATCH --mail-type=end
#SBATCH --job-name="Chip67"
#SBATCH --nodes=2
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=16G
#SBATCH --partition=all

##SBATCH --output=/path/to/outfile
##SBATCH --error=/path/to/errfile

# For array jobs
# Array job containing 100 tasks, run max 10 tasks at the same time
##SBATCH --array=1-100%10


#### Your shell commands below this line ####

## QIIME pipeline shell script ##

## Load QIIME module Version 1.9.1 ##
source activate qiime1

cd ~/Stoma_November_2018

# convert_fastaqual_fastq.py -c fastq_to_fastaqual -f Chip67.fastq -o fastaqual_Chip67
# split_libraries.py -m Chip67_map.txt -f fastaqual_Chip67/Chip67.fna -q fastaqual_Chip67/Chip67.qual -o split_output_Chip67 -b variable_length
split_libraries.py -m Chip67_map_5000.txt -f fastaqual_Chip67/Chip67.fna -q fastaqual_Chip67/Chip67.qual -o split_output_Chip67_5000 -b variable_length
