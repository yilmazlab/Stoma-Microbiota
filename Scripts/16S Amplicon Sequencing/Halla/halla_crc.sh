#!/bin/sh
#!/bin/bash
#SBATCH --job-name="halla"
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=20G
#SBATCH --cpus-per-task=10
#SBATCH --partition=epyc2

##SBATCH --output=/path/to/outfile
##SBATCH --error=/path/to/errfile

# For array jobs
# Array job containing 100 tasks, run max 10 tasks at the same time
##SBATCH --array=1-100%10


#### Your shell commands below this line ####

## QIIME pipeline shell script ##
export LC_CTYPE=en_US.UTF-8
export LC_ALL=en_US.UTF-8
source activate halla # since we used the latest silvadatase trained classifier

halla -X CRC_Bacteria_L6.txt -Y CRC_Metadata_X_Ileostoma.txt --output CRC_Halla_L6_CRC_Metadata_X_Ileostoma_0.05v2 --header -m nmi --fdr bh -q 0.05  --missing-char  NA  -e -1
# hallagram similarity_table.txt hypotheses_tree.txt associations.txt  --similarity NMI --cmap YlOrBr  --axlabels "Microbial species" "Phenotypes"  --outfile CRC_Halla_L2_Ileostoma_0.2.pdf
