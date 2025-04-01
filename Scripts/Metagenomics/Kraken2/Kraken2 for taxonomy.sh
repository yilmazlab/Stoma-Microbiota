#!/bin/sh
#!/bin/bash
#SBATCH --mail-user=bahtiyar.yilmaz@dbmr.unibe.ch
#SBATCH --mail-type=end
#SBATCH --job-name="Kraken2"
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=24
#SBATCH --partition=epyc2

##SBATCH --output=~/StomaMetagenomics/SlurmFiles
##SBATCH --error=~/StomaMetagenomics/SlurmFiles


export LC_CTYPE=en_US.UTF-8
export LC_ALL=en_US.UTF-8
source activate kraken2
### We built a database for Plant first
kraken2-build --download-taxonomy --db kraken2_plant_db
kraken2-build --download-library  plant --db kraken2_plant_db --threads 24
kraken2-build --build --db kraken2_plant_db --no-masking

### Then we build a Custom Database for the rest. Somehow plant DB cannot be included into this.
### It gave error and hence we kept it separately. 
kraken2-build --build --db CustomDB
kraken2-build --build --db CustomDB/library
kraken2-build --download-library bacteria  --db CustomDB
kraken2-build --download-library plasmid  --db CustomDB
kraken2-build --download-library archaea  --db CustomDB
kraken2-build --download-library viral --db CustomDB
kraken2-build --download-library human  --db CustomDB
kraken2-build --download-library fungi --db CustomDB
kraken2-build --download-library protozoa --db CustomDB
kraken2-build --build --max-db-size 96000000000 --db CustomDB
kraken2-inspect --db /storage/research/dbmr_yilmaz_gut_evo/CustomDB/ | head 10

# wget https://genome-idx.s3.amazonaws.com/kraken/k2_pluspfp_20210127.tar.gz #Alternatively one can use this library too


# Folders
cd /storage/homefs/terziev/StomaMetagenomics

# Running Kraken for Plant
kraken2 --use-names --db /storage/homefs/terziev/kraken2_plant_db  --threads 24 --minimum-base-quality 0  --confidence 0.2 --report-zero-counts --report /storage/homefs/terziev/StomaMetagenomics/kraken_analysis/SF${SLURM_ARRAY_TASK_ID}_plant.kreport --paired /storage/research/dbmr_yilmaz_gut_evo/Stoma/CLEAN_READS/SF${SLURM_ARRAY_TASK_ID}_R1.fastq  /storage/research/dbmr_yilmaz_gut_evo/Stoma/CLEAN_READS/SF${SLURM_ARRAY_TASK_ID}_R2.fastq  >  /storage/homefs/terziev/StomaMetagenomics/kraken_analysis/SF${SLURM_ARRAY_TASK_ID}_plant.kraken

# Then we run for the rest of it.
kraken2 --use-names --db /storage/homefs/terziev/CustomDB  --threads 24 --minimum-base-quality 0  --confidence 0.2 --report-zero-counts --report /storage/homefs/terziev/StomaMetagenomics/kraken_analysis/SF${SLURM_ARRAY_TASK_ID}_plant.kreport --paired /storage/research/dbmr_yilmaz_gut_evo/Stoma/CLEAN_READS/SF${SLURM_ARRAY_TASK_ID}_R1.fastq  /storage/research/dbmr_yilmaz_gut_evo/Stoma/CLEAN_READS/SF${SLURM_ARRAY_TASK_ID}_R2.fastq  >  /storage/homefs/terziev/StomaMetagenomics/kraken_analysis/SF${SLURM_ARRAY_TASK_ID}_Custom.kraken

source deactivate kraken2