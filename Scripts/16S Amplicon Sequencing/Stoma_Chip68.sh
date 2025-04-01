#!/bin/sh
#!/bin/bash
#SBATCH --mail-user=leven.terziev@dbmr.unibe.ch
#SBATCH --mail-type=end
#SBATCH --job-name="Chip68"
#SBATCH --nodes=2
#SBATCH --time=24:00:00
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
# convert_fastaqual_fastq.py -c fastq_to_fastaqual -f Chip68.fastq -o fastaqual_Chip68
# split_libraries.py -m Chip68_map.txt -f fastaqual_Chip68/Chip68.fna -q fastaqual_Chip68/Chip68.qual -o split_output_Chip68 -b variable_length
 split_libraries.py -m Chip68_map_5000.txt -f fastaqual_Chip68/Chip68.fna -q fastaqual_Chip68/Chip68.qual -o split_output_Chip68_5000 -b variable_length
