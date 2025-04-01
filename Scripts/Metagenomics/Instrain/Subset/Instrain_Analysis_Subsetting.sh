#!/bin/sh
#!/bin/bash
#SBATCH --mail-user=bahtiyar.yilmaz@dbmr.unibe.ch
#SBATCH --mail-type=end
#SBATCH --job-name="InStrain"
#SBATCH --nodes=1
#SBATCH --time=96:00:00
#SBATCH --mem-per-cpu=20G
#SBATCH --cpus-per-task=12
#SBATCH --partition=epyc2

##SBATCH --output=/path/to/outfile
##SBATCH --error=/path/to/errfile

#For array jobs
#Array job containing 10 tasks, run max 20 tasks at the same time
#SBATCH --array=1-26%20

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
module load UHTS/Analysis/sambamba/0.7.1

cd /storage/homefs/terziev/StomaMetagenomics/Instrain/

bowtie2 -p 10 -x /storage/homefs/terziev/UHGGv1/Reps_bt2/UHGG_reps.fasta.bt2 -1 /storage/homefs/terziev/StomaMetagenomics/CLEAN_READS/SF${SLURM_ARRAY_TASK_ID}_1.fastq -2 /storage/homefs/terziev/StomaMetagenomics/CLEAN_READS/SF${SLURM_ARRAY_TASK_ID}_2.fastq > /storage/homefs/terziev/StomaMetagenomics/Instrain/SF${SLURM_ARRAY_TASK_ID}.sam

for i in {1..26};
  do samtools view -hS -F 4 /storage/homefs/terziev/StomaMetagenomics/Instrain/SF${i}.sam | awk '$7 == "=" || $1 ~ /^@/ { print $0 }' | samtools sort -m 10G -o /storage/homefs/terziev/StomaMetagenomics/Instrain/SF${i}.bam -O bam
    done

  for i in {1..26}; 
  	do samtools coverage /storage/homefs/terziev/StomaMetagenomics/Instrain/SF${i}.sorted.bam -o /storage/homefs/terziev/StomaMetagenomics/Instrain/Coverage/SF${i}.coverage.txt
   	 done

### Subsetting to sample to 56x mean depth
sambamba view -f bam -t 10 --subsampling-seed=3 -s 0.330 SF2.sorted.bam -o E.coli_SF2_sub.bam
sambamba view -f bam -t 10 --subsampling-seed=3 -s 0.220 SF3.sorted.bam -o E.coli_SF3_sub.bam
sambamba view -f bam -t 10 --subsampling-seed=3 -s 0.032 SF12.sorted.bam -o E.coli_SF12_sub.bam
sambamba view -f bam -t 10 --subsampling-seed=3 -s 1.000 SF14.sorted.bam -o E.coli_SF14_sub.bam
sambamba view -f bam -t 10 --subsampling-seed=3 -s 1.000 SF17.sorted.bam -o E.coli_SF17_sub.bam
sambamba view -f bam -t 10 --subsampling-seed=3 -s 0.200 SF20.sorted.bam -o E.coli_SF20_sub.bam
sambamba view -f bam -t 10 --subsampling-seed=3 -s 0.910 SF21.sorted.bam -o E.coli_SF21_sub.bam
sambamba view -f bam -t 10 --subsampling-seed=3 -s 1.000 SF24.sorted.bam -o E.coli_SF24_sub.bam
sambamba view -f bam -t 10 --subsampling-seed=3 -s 0.140 SF25.sorted.bam -o E.coli_SF25_sub.bam
sambamba view -f bam -t 10 --subsampling-seed=3 -s 0.680 SF26.sorted.bam -o E.coli_SF26_sub.bam

### Subsetting to sample to 60x mean depth ### 
sambamba view -f bam -t 10 --subsampling-seed=3 -s 0.037 SF2.sorted.bam -o Klebsiella_SF2_sub.bam
sambamba view -f bam -t 10 --subsampling-seed=3 -s 0.028 SF3.sorted.bam -o Klebsiella_SF3_sub.bam
sambamba view -f bam -t 10 --subsampling-seed=3 -s 0.034 SF9.sorted.bam -o Klebsiella_SF9_sub.bam
sambamba view -f bam -t 10 --subsampling-seed=3 -s 0.028 SF10.sorted.bam -o Klebsiella_SF10_sub.bam
sambamba view -f bam -t 10 --subsampling-seed=3 -s 0.480 SF20.sorted.bam -o Klebsiella_SF20_sub.bam
sambamba view -f bam -t 10 --subsampling-seed=3 -s 0.950 SF21.sorted.bam -o Klebsiella_SF21_sub.bam
sambamba view -f bam -t 10 --subsampling-seed=3 -s 1.000 SF23.sorted.bam -o Klebsiella_SF23_sub.bam
sambamba view -f bam -t 10 --subsampling-seed=3 -s 1.000 SF25.sorted.bam -o Klebsiella_SF25_sub.bam
sambamba view -f bam -t 10 --subsampling-seed=3 -s 1.000 SF26.sorted.bam -o Klebsiella_SF26_sub.bam


### Checking the Coverage for each sample after subsetting ### 
 mkdir Coverage

for i in {2,3,12,14,17,20,21,24,25,26}; 
 do samtools coverage /storage/homefs/terziev/StomaMetagenomics/Instrain/E.coli_SF${i}_sub.bam -o /storage/homefs/terziev/StomaMetagenomics/Instrain/Coveragee/E.coli_SF${i}_sub.coverage.txt
 done

for i in {2,3,9,10,20,21,23,25,26}; 
 do samtools coverage /storage/homefs/terziev/StomaMetagenomics/Instrain/Klebsiella_SF${i}_sub.bam -o /storage/homefs/terziev/StomaMetagenomics/Instrain/Coverage/Klebsiella_SF${i}_sub.coverage.txt
 done

### Running the InStrain profile using subsetted bam files to characterizethe nucleotide diversity, SNSs and SNVs, linkage, genome or gene level metrics ### 
## For E.coli
inStrain profile /storage/homefs/terziev/StomaMetagenomics/Instrain/E.coli_SF${SLURM_ARRAY_TASK_ID}_sub.bam /storage/homefs/terziev/UHGGv1/UHGG_reps.fasta -o /storage/homefs/terziev/StomaMetagenomics/Instrain/E.coli_sub/E.coli_SF${SLURM_ARRAY_TASK_ID}.IS -p 24 -g /storage/homefs/terziev/UHGGv1/UHGG_reps.genes.fna -s /storage/homefs/terziev/UHGGv1/UHGG_reps.stb --database_mode
inStrain compare -i StomaMetagenomics/Instrain/E.coli_sub/E.coli_SF2.IS StomaMetagenomics/Instrain/E.coli_sub/E.coli_SF3.IS StomaMetagenomics/Instrain/E.coli_sub/E.coli_SF12.IS StomaMetagenomics/Instrain/E.coli_sub/E.coli_SF14.IS StomaMetagenomics/Instrain/E.coli_sub/E.coli_SF17.IS StomaMetagenomics/Instrain/E.coli_sub/E.coli_SF20.IS StomaMetagenomics/Instrain/E.coli_sub/E.coli_SF21.IS StomaMetagenomics/Instrain/E.coli_sub/E.coli_SF25.IS StomaMetagenomics/Instrain/E.coli_sub/E.coli_SF26.IS -s UHGGv1/UHGG_reps.stb -p 24 -o StomaMetagenomics/Instrain/E.coli_sub/IS.COMPARE --database_mode

## For Klebsiella
inStrain profile /storage/homefs/terziev/StomaMetagenomics/Instrain/Klebsiella_SF${SLURM_ARRAY_TASK_ID}_sub.bam /storage/homefs/terziev/UHGGv1/UHGG_reps.fasta -o /storage/homefs/terziev/StomaMetagenomics/Instrain/Klebsiella_sub/Klebsiella_SF${SLURM_ARRAY_TASK_ID}.IS -p 24 -g /storage/homefs/terziev/UHGGv1/UHGG_reps.genes.fna -s /storage/homefs/terziev/UHGGv1/UHGG_reps.stb --database_mode
inStrain compare -i StomaMetagenomics/Instrain/Klebsiella_sub/Klebsiella_SF2.IS StomaMetagenomics/Instrain/Klebsiella_sub/Klebsiella_SF3.IS StomaMetagenomics/Instrain/Klebsiella_sub/Klebsiella_SF9.IS StomaMetagenomics/Instrain/Klebsiella_sub/Klebsiella_SF10.IS StomaMetagenomics/Instrain/Klebsiella_sub/Klebsiella_SF20.IS StomaMetagenomics/Instrain/Klebsiella_sub/Klebsiella_SF21.IS StomaMetagenomics/Instrain/Klebsiella_sub/Klebsiella_SF23.IS StomaMetagenomics/Instrain/Klebsiella_sub/Klebsiella_SF25.IS StomaMetagenomics/Instrain/Klebsiella_sub/Klebsiella_SF26.IS -s UHGGv1/UHGG_reps.stb -p 24 -o StomaMetagenomics/Instrain/E.coli_sub/IS.COMPARE --database_mode
