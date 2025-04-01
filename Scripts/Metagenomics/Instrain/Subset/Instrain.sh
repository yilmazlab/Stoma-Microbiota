#!/bin/sh
#!/bin/bash
#SBATCH --mail-user=bahtiyar.yilmaz@dbmr.unibe.ch
#SBATCH --mail-type=end
#SBATCH --job-name="Instrain"
#SBATCH --nodes=1
#SBATCH --time=96:00:00
#SBATCH --mem-per-cpu=10G
#SBATCH --cpus-per-task=24
#SBATCH --partition=epyc2

##SBATCH --output=/path/to/outfile
##SBATCH --error=/path/to/errfile

#For array jobs
#Array job containing 10 tasks, run max 20 tasks at the same time
### SBATCH --array=1-43%20

#### Your shell commands below this line ####
#### Tutorial: https://github.com/bxlab/metaWRAP/blob/master/Usage_tutorial.md

export LC_CTYPE=en_US.UTF-8
export LC_ALL=en_US.UTF-8
module load vital-it
module add UHTS/Aligner/bowtie2/2.3.4.1
module add SequenceAnalysis/GenePrediction/prodigal/2.6.3
module add UHTS/Analysis/Mash/2.0
module add UHTS/Analysis/mummer/4.0.0beta1
module add UHTS/Analysis/fastANI/1.1
module add UHTS/Analysis/samtools/1.10


bowtie2 -p 10 -x /storage/homefs/terziev/UHGGv1/Reps_bt2/UHGG_reps.fasta.bt2 -1 /storage/homefs/terziev/StomaMetagenomics/CLEAN_READS/SF${SLURM_ARRAY_TASK_ID}_1.fastq -2 /storage/homefs/terziev/StomaMetagenomics/CLEAN_READS/SF${SLURM_ARRAY_TASK_ID}_2.fastq > /storage/homefs/terziev/StomaMetagenomics/Instrain/SF${SLURM_ARRAY_TASK_ID}.sam
inStrain profile /storage/homefs/terziev/StomaMetagenomics/Instrain/SF${SLURM_ARRAY_TASK_ID}.sam UHGGv1/UHGG_reps.fasta -o /storage/homefs/terziev/StomaMetagenomics/Instrain/SF${SLURM_ARRAY_TASK_ID}.IS -p 24 -g UHGGv1/UHGG_reps.genes.fna -s UHGGv1/UHGG_reps.stb --database_mode
inStrain compare -i /storage/homefs/terziev/StomaMetagenomics/Instrain/SF1.IS /storage/homefs/terziev/StomaMetagenomics/Instrain/SF10.IS /storage/homefs/terziev/StomaMetagenomics/Instrain/SF11.IS /storage/homefs/terziev/StomaMetagenomics/Instrain/SF12.IS /storage/homefs/terziev/StomaMetagenomics/Instrain/SF13.IS /storage/homefs/terziev/StomaMetagenomics/Instrain/SF14.IS /storage/homefs/terziev/StomaMetagenomics/Instrain/SF15.IS /storage/homefs/terziev/StomaMetagenomics/Instrain/SF16.IS /storage/homefs/terziev/StomaMetagenomics/Instrain/SF17.IS /storage/homefs/terziev/StomaMetagenomics/Instrain/SF18.IS /storage/homefs/terziev/StomaMetagenomics/Instrain/SF19.IS /storage/homefs/terziev/StomaMetagenomics/Instrain/SF2.IS /storage/homefs/terziev/StomaMetagenomics/Instrain/SF20.IS /storage/homefs/terziev/StomaMetagenomics/Instrain/SF21.IS /storage/homefs/terziev/StomaMetagenomics/Instrain/SF22.IS /storage/homefs/terziev/StomaMetagenomics/Instrain/SF23.IS /storage/homefs/terziev/StomaMetagenomics/Instrain/SF24.IS /storage/homefs/terziev/StomaMetagenomics/Instrain/SF25.IS /storage/homefs/terziev/StomaMetagenomics/Instrain/SF26.IS /storage/homefs/terziev/StomaMetagenomics/Instrain/SF28.IS /storage/homefs/terziev/StomaMetagenomics/Instrain/SF29.IS /storage/homefs/terziev/StomaMetagenomics/Instrain/SF3.IS /storage/homefs/terziev/StomaMetagenomics/Instrain/SF30.IS /storage/homefs/terziev/StomaMetagenomics/Instrain/SF31.IS /storage/homefs/terziev/StomaMetagenomics/Instrain/SF32.IS /storage/homefs/terziev/StomaMetagenomics/Instrain/SF33.IS /storage/homefs/terziev/StomaMetagenomics/Instrain/SF34.IS /storage/homefs/terziev/StomaMetagenomics/Instrain/SF35.IS /storage/homefs/terziev/StomaMetagenomics/Instrain/SF36.IS /storage/homefs/terziev/StomaMetagenomics/Instrain/SF38.IS /storage/homefs/terziev/StomaMetagenomics/Instrain/SF39.IS /storage/homefs/terziev/StomaMetagenomics/Instrain/SF4.IS /storage/homefs/terziev/StomaMetagenomics/Instrain/SF40.IS /storage/homefs/terziev/StomaMetagenomics/Instrain/SF41.IS /storage/homefs/terziev/StomaMetagenomics/Instrain/SF42.IS /storage/homefs/terziev/StomaMetagenomics/Instrain/SF43.IS /storage/homefs/terziev/StomaMetagenomics/Instrain/SF6.IS /storage/homefs/terziev/StomaMetagenomics/Instrain/SF7.IS /storage/homefs/terziev/StomaMetagenomics/Instrain/SF8.IS /storage/homefs/terziev/StomaMetagenomics/Instrain/SF9.IS -s UHGGv1/UHGG_reps.stb -p 24 -o /storage/homefs/terziev/StomaMetagenomics/Instrain/IS.COMPARE --database_mode
