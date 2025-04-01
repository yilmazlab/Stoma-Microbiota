#!/bin/sh
#!/bin/bash
#SBATCH --mail-user=bahtiyar.yilmaz@dbmr.unibe.ch
#SBATCH --mail-type=end
#SBATCH --job-name="SGV"
#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=10G
#SBATCH --cpus-per-task=5
#SBATCH --partition=epyc2

##SBATCH --output=/path/to/outfile
##SBATCH --error=/path/to/errfile

#For array jobs
#Array job containing 10 tasks, run max 20 tasks at the same time
###SBATCH --array=1-43%20

export LC_CTYPE=en_US.UTF-8
export LC_ALL=en_US.UTF-8

cd /storage/research/dbmr_yilmaz_gut_evo/
source activate SGV

# wget https://zenodo.org/record/3237975/files/DataFiles.tar.gz
# gzip -dc DataFiles.tar.gz > DataFiles
# git clone https://github.com/segalab/SGVFinder.git
# mv DataFiles SGVFinder/
# cd /storage/research/dbmr_yilmaz_gut_evo/SGVFinder/src/cy_ext/
# python setup.py build_ext

# conda create --name SGV python=2.7.8
# source activate SGV
# pip install numpy==1.14.2
# pip install scipy==1.1.0
# pip install bokeh==0.12.6
# pip install pandas==0.23.4
# pip install ujson==1.35
# pip install biopython==1.68
# pip install Cython==0.29.30
# gzip -9 /storage/research/dbmr_yilmaz_gut_evo/Stoma/RAW_READS/P1/SF1_R1.fastq

# Run ICRA
mkdir -p Stoma/icra
python /storage/research/dbmr_yilmaz_gut_evo/SGVFinder/src/ICRA_cmd.py Stoma/icra /storage/homefs/terziev/StomaMetagenomics/CLEAN_READS/Used/SF${SLURM_ARRAY_TASK_ID} --pe

# Run SGVF_PerFile
mkdir -p Stoma/SGVF_PerFile

python /storage/research/dbmr_yilmaz_gut_evo/SGVFinder/src/SGVF_PerFile_cmd.py /storage/research/dbmr_yilmaz_gut_evo/Stoma/icra/SF${SLURM_ARRAY_TASK_ID}.jsdel /storage/research/dbmr_yilmaz_gut_evo/Stoma/SGVF_PerFile/input.map.jsdel 100 --x_coverage 0.01 --rate_param 10

for i in {1..31}; 
	do python /storage/research/dbmr_yilmaz_gut_evo/SGVFinder/src/SGVF_PerFile_cmd.py /storage/research/dbmr_yilmaz_gut_evo/Stoma/icra/SF${i}.jsdel /storage/research/dbmr_yilmaz_gut_evo/Stoma/SGVF_PerFile/SF${i}.map.jsdel 100 --x_coverage 0.01 --rate_param 10
done

# python /storage/research/dbmr_yilmaz_gut_evo/SGVFinder/src/SGVF_PerFile_cmd.py /storage/research/dbmr_yilmaz_gut_evo/Stoma/icra/SF25.jsdel /storage/research/dbmr_yilmaz_gut_evo/Stoma/SGVF_PerFile/SF25.map.jsdel 100 --x_coverage 0.01 --rate_param 10

mkdir -p Stoma/results
mkdir -p Stoma/results/html
python /storage/research/dbmr_yilmaz_gut_evo/SGVFinder/src/SGVF_cmd.py Stoma/input_glob_string3/*.jsdel /storage/research/dbmr_yilmaz_gut_evo/Stoma/results/dsgv.csv /storage/research/dbmr_yilmaz_gut_evo/Stoma/results/vsgv.csv --min_samp_cutoff 5 --x_coverage 0.01 --rate_param 10 --browser_path /storage/research/dbmr_yilmaz_gut_evo/Stoma/results/html --csv_output

conda activate sgvfinder
mkdir -p /data/umcg-dwang/project/SV_test/s02.SVs_pipeline_cohort.full.res/test
mkdir -p /data/umcg-dwang/project/SV_test/s02.SVs_pipeline_cohort.full.res/test/test.html
/data/umcg-dwang/opt/Anaconda3/envs/sgvfinder/bin/python /data/umcg-dwang/software/SGVFinder/src/SGVF_cmd.py "/data/umcg-dwang/project/SV_test/input_jsdel/*.jsdel" /data/umcg-dwang/project/SV_test/s02.SVs_pipeline_cohort.full.res/test/test.dsgv.csv /data/umcg-dwang/project/SV_test/s02.SVs_pipeline_cohort.full.res/test/test.vsgv.csv --min_samp_cutoff 10 --x_coverage 0.01 --rate_param 10 --browser_path /data/umcg-dwang/project/SV_test/s02.SVs_pipeline_cohort.full.res/test/test.html --csv_output

cd /storage/research/dbmr_yilmaz_gut_evo/Stoma
 mkdir -p SV_test/s02.SVs_pipeline_cohort.full.res/test
 mkdir -p SV_test/s02.SVs_pipeline_cohort.full.res/test/test.html

mkdir -p /data/umcg-dwang/project/SV_test/s02.SVs_pipeline_cohort.full.res/test
mkdir -p /data/umcg-dwang/project/SV_test/s02.SVs_pipeline_cohort.full.res/test/test.html

python /storage/research/dbmr_yilmaz_gut_evo/SGVFinder/src/SGVF_cmd.py /storage/research/dbmr_yilmaz_gut_evo/Stoma/SV_test/input_jsdel/SV_test/input_jsdel/*.jsdel /storage/research/dbmr_yilmaz_gut_evo/Stoma/SV_test/s02.SVs_pipeline_cohort.full.res/test/test.dsgv.csv /storage/research/dbmr_yilmaz_gut_evo/Stoma/SV_test/s02.SVs_pipeline_cohort.full.res/test/test.vsgv.csv --min_samp_cutoff 10 --x_coverage 0.01 --rate_param 10 --browser_path /storage/research/dbmr_yilmaz_gut_evo/Stoma/SV_test/s02.SVs_pipeline_cohort.full.res/test/test.html --csv_output






















