#!/bin/sh
#!/bin/bash
#SBATCH --mail-user=bahtiyar.yilmaz@dbmr.unibe.ch
#SBATCH --mail-type=end
#SBATCH --job-name="MIDAS"
#SBATCH --nodes=1
#SBATCH --time=96:00:00
#SBATCH --mem-per-cpu=10G
#SBATCH --cpus-per-task=24
#SBATCH --partition=epyc2

##SBATCH --output=/path/to/outfile
##SBATCH --error=/path/to/errfile

#For array jobs
#Array job containing 10 tasks, run max 20 tasks at the same time
###### SBATCH --array=1-43%20

#### Your shell commands below this line ####
#### Tutorial: https://github.com/bxlab/metaWRAP/blob/master/Usage_tutorial.md

export LC_CTYPE=en_US.UTF-8
export LC_ALL=en_US.UTF-8

module load vital-it
module add UHTS/Aligner/bowtie2/2.3.4.1

export PYTHONPATH=$PYTHONPATH:/storage/homefs/terziev/MIDAS
export PATH="/storage/homefs/terziev/MIDAS/:$PATH"
export PATH="/storage/homefs/terziev/MIDAS/scripts/:$PATH"
export PATH="/storage/homefs/terziev/midas_database/:$PATH"

cd StomaMetagenomics
run_midas.py species MIDAS_Output/SF${SLURM_ARRAY_TASK_ID} -d /storage/homefs/terziev/MIDAS/midas_db_v1.2 -t 24 -1 CLEAN_READS/SF${SLURM_ARRAY_TASK_ID}_1.fastq  -2 CLEAN_READS/SF${SLURM_ARRAY_TASK_ID}_2.fastq
run_midas.py genes MIDAS_Output/SF${SLURM_ARRAY_TASK_ID} -d /storage/homefs/terziev/MIDAS/midas_db_v1.2 -t 24 -1 CLEAN_READS/SF${SLURM_ARRAY_TASK_ID}_1.fastq  -2 CLEAN_READS/SF${SLURM_ARRAY_TASK_ID}_2.fastq
run_midas.py snps MIDAS_Output/SF${SLURM_ARRAY_TASK_ID} -d /storage/homefs/terziev/MIDAS/midas_db_v1.2 -t 24 -1 CLEAN_READS/SF${SLURM_ARRAY_TASK_ID}_1.fastq  -2 CLEAN_READS/SF${SLURM_ARRAY_TASK_ID}_2.fastq


merge_midas.py species  /storage/homefs/terziev/StomaMetagenomics/MIDAS_Output/Merged/ -i  /storage/homefs/terziev/StomaMetagenomics/MIDAS_Output/sample_paths.txt -t file --sample_depth 5.0 -d /storage/homefs/terziev/MIDAS/midas_db_v1.2
merge_midas.py genes /storage/homefs/terziev/StomaMetagenomics/MIDAS_Output/Merged/Genes/ -d /storage/homefs/terziev/MIDAS/midas_db_v1.2  -i /storage/homefs/terziev/StomaMetagenomics/MIDAS_Output/sample_paths.txt -t file -t file --sample_depth 5.0
merge_midas.py snps /storage/homefs/terziev/StomaMetagenomics/MIDAS_Output/Merged/SPNS/ Klebsiella_pneumoniae_54788v2 ---species_id Klebsiella_pneumoniae_54788 d /storage/homefs/terziev/MIDAS/midas_db_v1.2  -i /storage/homefs/terziev/StomaMetagenomics/MIDAS_Output/sample_paths.txt -t file -t file --sample_depth 5.0 --site_depth 5 --site_prev 0.3 --all_snps --allele_freq 0.01




merge_midas.py snps /storage/homefs/terziev/StomaMetagenomics/MIDAS_Output/Merged/SPNS/Klebsiella_pneumoniae_54788v2 --species_id Klebsiella_pneumoniae_54788 -d /storage/homefs/terziev/MIDAS/midas_db_v1.2  -i /storage/homefs/terziev/StomaMetagenomics/MIDAS_Output/sample_paths.txt -t file
merge_midas.py snps /storage/homefs/terziev/StomaMetagenomics/MIDAS_Output/Merged/SPNS/Escherichia_coli_58110v2 --species_id Escherichia_coli_58110 -d /storage/homefs/terziev/MIDAS/midas_db_v1.2  -i /storage/homefs/terziev/StomaMetagenomics/MIDAS_Output/sample_paths.txt -t file
merge_midas.py snps /storage/homefs/terziev/StomaMetagenomics/MIDAS_Output/Merged/SPNS/Bacteroides_fragilis_54507v2 --species_id Bacteroides_fragilis_54507 -d /storage/homefs/terziev/MIDAS/midas_db_v1.2  -i /storage/homefs/terziev/StomaMetagenomics/MIDAS_Output/sample_paths.txt -t file
merge_midas.py snps /storage/homefs/terziev/StomaMetagenomics/MIDAS_Output/Merged/SPNS/Clostridium_clostridioforme_51842v2 --species_id Clostridium_clostridioforme_51842 -d /storage/homefs/terziev/MIDAS/midas_db_v1.2  -i /storage/homefs/terziev/StomaMetagenomics/MIDAS_Output/sample_paths.txt -t file
merge_midas.py snps /storage/homefs/terziev/StomaMetagenomics/MIDAS_Output/Merged/SPNS/Veillonella_atypica_58169v2 --species_id Veillonella_atypica_58169 -d /storage/homefs/terziev/MIDAS/midas_db_v1.2  -i /storage/homefs/terziev/StomaMetagenomics/MIDAS_Output/sample_paths.txt -t file
merge_midas.py snps /storage/homefs/terziev/StomaMetagenomics/MIDAS_Output/Merged/SPNS/Citrobacter_freundii_56148v2 --species_id Citrobacter_freundii_56148 -d /storage/homefs/terziev/MIDAS/midas_db_v1.2  -i /storage/homefs/terziev/StomaMetagenomics/MIDAS_Output/sample_paths.txt -t file
merge_midas.py snps /storage/homefs/terziev/StomaMetagenomics/MIDAS_Output/Merged/SPNS/Morganella_morganii_57363v2 --species_id Morganella_morganii_57363 -d /storage/homefs/terziev/MIDAS/midas_db_v1.2  -i /storage/homefs/terziev/StomaMetagenomics/MIDAS_Output/sample_paths.txt -t file
merge_midas.py snps /storage/homefs/terziev/StomaMetagenomics/MIDAS_Output/Merged/SPNS/Hafnia_alvei_56921v2 --species_id Hafnia_alvei_56921 -d /storage/homefs/terziev/MIDAS/midas_db_v1.2  -i /storage/homefs/terziev/StomaMetagenomics/MIDAS_Output/sample_paths.txt -t file
merge_midas.py snps /storage/homefs/terziev/StomaMetagenomics/MIDAS_Output/Merged/SPNS/Bacteroides_vulgatus_57955v2 --species_id Bacteroides_vulgatus_57955 -d /storage/homefs/terziev/MIDAS/midas_db_v1.2  -i /storage/homefs/terziev/StomaMetagenomics/MIDAS_Output/sample_paths.txt -t file
merge_midas.py snps /storage/homefs/terziev/StomaMetagenomics/MIDAS_Output/Merged/SPNS/Clostridium_butyricum_56361v2 --species_id Clostridium_butyricum_56361 -d /storage/homefs/terziev/MIDAS/midas_db_v1.2  -i /storage/homefs/terziev/StomaMetagenomics/MIDAS_Output/sample_paths.txt -t file
merge_midas.py snps /storage/homefs/terziev/StomaMetagenomics/MIDAS_Output/Merged/SPNS/Escherichia_fergusonii_56914v2 --species_id Escherichia_fergusonii_56914 -d /storage/homefs/terziev/MIDAS/midas_db_v1.2  -i /storage/homefs/terziev/StomaMetagenomics/MIDAS_Output/sample_paths.txt -t file
merge_midas.py snps /storage/homefs/terziev/StomaMetagenomics/MIDAS_Output/Merged/SPNS/Haemophilus_parainfluenzae_62468v2 --species_id Haemophilus_parainfluenzae_62468 -d /storage/homefs/terziev/MIDAS/midas_db_v1.2  -i /storage/homefs/terziev/StomaMetagenomics/MIDAS_Output/sample_paths.txt -t file
merge_midas.py snps /storage/homefs/terziev/StomaMetagenomics/MIDAS_Output/Merged/SPNS/Prevotella_timonensis_58235v2 --species_id Prevotella_timonensis_58235 -d /storage/homefs/terziev/MIDAS/midas_db_v1.2  -i /storage/homefs/terziev/StomaMetagenomics/MIDAS_Output/sample_paths.txt -t file
merge_midas.py snps /storage/homefs/terziev/StomaMetagenomics/MIDAS_Output/Merged/SPNS/Porphyromonas_asaccharolytica_56465v2 --species_id Porphyromonas_asaccharolytica_56465 -d /storage/homefs/terziev/MIDAS/midas_db_v1.2  -i /storage/homefs/terziev/StomaMetagenomics/MIDAS_Output/sample_paths.txt -t file
merge_midas.py snps /storage/homefs/terziev/StomaMetagenomics/MIDAS_Output/Merged/SPNS/Porphyromonas_bennonis_58999v2 --species_id Porphyromonas_bennonis_58999 -d /storage/homefs/terziev/MIDAS/midas_db_v1.2  -i /storage/homefs/terziev/StomaMetagenomics/MIDAS_Output/sample_paths.txt -t file
merge_midas.py snps /storage/homefs/terziev/StomaMetagenomics/MIDAS_Output/Merged/SPNS/Citrobacter_braakii_56022v2 --species_id Citrobacter_braakii_56022 -d /storage/homefs/terziev/MIDAS/midas_db_v1.2  -i /storage/homefs/terziev/StomaMetagenomics/MIDAS_Output/sample_paths.txt -t file
merge_midas.py snps /storage/homefs/terziev/StomaMetagenomics/MIDAS_Output/Merged/SPNS/Citrobacter_amalonaticus_54491v2 --species_id Citrobacter_amalonaticus_54491 -d /storage/homefs/terziev/MIDAS/midas_db_v1.2  -i /storage/homefs/terziev/StomaMetagenomics/MIDAS_Output/sample_paths.txt -t file
merge_midas.py snps /storage/homefs/terziev/StomaMetagenomics/MIDAS_Output/Merged/SPNS/Klebsiella_oxytoca_56762v2 --species_id Klebsiella_oxytoca_56762 -d /storage/homefs/terziev/MIDAS/midas_db_v1.2  -i /storage/homefs/terziev/StomaMetagenomics/MIDAS_Output/sample_paths.txt -t file
merge_midas.py snps /storage/homefs/terziev/StomaMetagenomics/MIDAS_Output/Merged/SPNS/Haemophilus_parainfluenzae_57123v2 --species_id Haemophilus_parainfluenzae_57123 -d /storage/homefs/terziev/MIDAS/midas_db_v1.2  -i /storage/homefs/terziev/StomaMetagenomics/MIDAS_Output/sample_paths.txt -t file
merge_midas.py snps /storage/homefs/terziev/StomaMetagenomics/MIDAS_Output/Merged/SPNS/Prevotella_buccalis_62098v2 --species_id Prevotella_buccalis_62098 -d /storage/homefs/terziev/MIDAS/midas_db_v1.2  -i /storage/homefs/terziev/StomaMetagenomics/MIDAS_Output/sample_paths.txt -t file
merge_midas.py snps /storage/homefs/terziev/StomaMetagenomics/MIDAS_Output/Merged/SPNS/Haemophilus_parainfluenzae_62356v2 --species_id Haemophilus_parainfluenzae_62356 -d /storage/homefs/terziev/MIDAS/midas_db_v1.2  -i /storage/homefs/terziev/StomaMetagenomics/MIDAS_Output/sample_paths.txt -t file
merge_midas.py snps /storage/homefs/terziev/StomaMetagenomics/MIDAS_Output/Merged/SPNS/Streptococcus_salivarius_58024v2 --species_id Streptococcus_salivarius_58024-d /storage/homefs/terziev/MIDAS/midas_db_v1.2  -i /storage/homefs/terziev/StomaMetagenomics/MIDAS_Output/sample_paths.txt -t file
merge_midas.py snps /storage/homefs/terziev/StomaMetagenomics/MIDAS_Output/Merged/SPNS/Clostridium_sp_60465v2 --species_id Clostridium_sp_60465 -d /storage/homefs/terziev/MIDAS/midas_db_v1.2  -i /storage/homefs/terziev/StomaMetagenomics/MIDAS_Output/sample_paths.txt -t file
merge_midas.py snps /storage/homefs/terziev/StomaMetagenomics/MIDAS_Output/Merged/SPNS/Veillonella_dispar_61763v2 --species_id Veillonella_dispar_61763 -d /storage/homefs/terziev/MIDAS/midas_db_v1.2  -i /storage/homefs/terziev/StomaMetagenomics/MIDAS_Output/sample_paths.txt -t file
merge_midas.py snps /storage/homefs/terziev/StomaMetagenomics/MIDAS_Output/Merged/SPNS/Alistipes_onderdonkii_55464v2 --species_id Alistipes_onderdonkii_55464 -d /storage/homefs/terziev/MIDAS/midas_db_v1.2  -i /storage/homefs/terziev/StomaMetagenomics/MIDAS_Output/sample_paths.txt -t file
merge_midas.py snps /storage/homefs/terziev/StomaMetagenomics/MIDAS_Output/Merged/SPNS/Anaerostipes_hadrus_55206v2 --species_id Anaerostipes_hadrus_55206 -d /storage/homefs/terziev/MIDAS/midas_db_v1.2  -i /storage/homefs/terziev/StomaMetagenomics/MIDAS_Output/sample_paths.txt -t file
merge_midas.py snps /storage/homefs/terziev/StomaMetagenomics/MIDAS_Output/Merged/SPNS/Ruminococcus_gnavus_57638v2 --species_id Ruminococcus_gnavus_57638 -d /storage/homefs/terziev/MIDAS/midas_db_v1.2  -i /storage/homefs/terziev/StomaMetagenomics/MIDAS_Output/sample_paths.txt -t file
merge_midas.py snps /storage/homefs/terziev/StomaMetagenomics/MIDAS_Output/Merged/SPNS/Enterobacter_cloacae_58148v2 --species_id Enterobacter_cloacae_58148 -d /storage/homefs/terziev/MIDAS/midas_db_v1.2  -i /storage/homefs/terziev/StomaMetagenomics/MIDAS_Output/sample_paths.txt -t file
merge_midas.py snps /storage/homefs/terziev/StomaMetagenomics/MIDAS_Output/Merged/SPNS/Citrobacter_sp_58515v2 --species_id Citrobacter_sp_58515 -d /storage/homefs/terziev/MIDAS/midas_db_v1.2  -i /storage/homefs/terziev/StomaMetagenomics/MIDAS_Output/sample_paths.txt -t file
merge_midas.py snps /storage/homefs/terziev/StomaMetagenomics/MIDAS_Output/Merged/SPNS/Porphyromonas_uenonis_49360v2 --species_id Porphyromonas_uenonis_49360 -d /storage/homefs/terziev/MIDAS/midas_db_v1.2  -i /storage/homefs/terziev/StomaMetagenomics/MIDAS_Output/sample_paths.txt -t file
merge_midas.py snps /storage/homefs/terziev/StomaMetagenomics/MIDAS_Output/Merged/SPNS/Ruminococcus_torques_62045v2 --species_id Ruminococcus_torques_62045 -d /storage/homefs/terziev/MIDAS/midas_db_v1.2  -i /storage/homefs/terziev/StomaMetagenomics/MIDAS_Output/sample_paths.txt -t file
merge_midas.py snps /storage/homefs/terziev/StomaMetagenomics/MIDAS_Output/Merged/SPNS/Prevotella_copri_61740v2 --species_id Prevotella_copri_61740 -d /storage/homefs/terziev/MIDAS/midas_db_v1.2  -i /storage/homefs/terziev/StomaMetagenomics/MIDAS_Output/sample_paths.txt -t file
merge_midas.py snps /storage/homefs/terziev/StomaMetagenomics/MIDAS_Output/Merged/SPNS/Bacteroides_caccae_53434v2 --species_id Bacteroides_caccae_53434 -d /storage/homefs/terziev/MIDAS/midas_db_v1.2  -i /storage/homefs/terziev/StomaMetagenomics/MIDAS_Output/sample_paths.txt -t file
merge_midas.py snps /storage/homefs/terziev/StomaMetagenomics/MIDAS_Output/Merged/SPNS/Streptococcus_vestibularis_56030v2 --species_id Streptococcus_vestibularis_56030 -d /storage/homefs/terziev/MIDAS/midas_db_v1.2  -i /storage/homefs/terziev/StomaMetagenomics/MIDAS_Output/sample_paths.txt -t file


# Strain Tracking
## Step 1
### These scripts will allow you to identify rare SNPs that discriminate individual strains and to track these SNPs between hosts to elucidate transmission patterns.
# strain_tracking.py id_markers --indir /storage/homefs/terziev/StomaMetagenomics/MIDAS_Output/Merged/SPNS/Escherichia_coli_58110 --out /storage/homefs/terziev/StomaMetagenomics/MIDAS_Output/Merged/SPNS/58110.markers --samples SF2,SF3,SF5

## Step 2
## Track rare SNPs between samples and determine transmission
# strain_tracking.py track_markers --indir /storage/homefs/terziev/StomaMetagenomics/MIDAS_Output/Merged/SPNS/Escherichia_coli_58110 --out /storage/homefs/terziev/StomaMetagenomics/MIDAS_Output/Merged/SPNS/58110.marker_sharing.txt --markers /storage/homefs/terziev/StomaMetagenomics/MIDAS_Output/Merged/SPNS/58110.markers

# call_consensus.py /storage/homefs/terziev/StomaMetagenomics/MIDAS_Output/Merged/SPNS/Escherichia_coli_58110/ --out /storage/homefs/terziev/StomaMetagenomics/MIDAS_Output/Merged/SPNS/58110.Consensus.txt
