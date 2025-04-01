#!/bin/sh
#!/bin/bash
#SBATCH --mail-user=leven.terziev@dbmr.unibe.ch
#SBATCH --mail-type=end
#SBATCH --job-name="Chip67_68"
#SBATCH --nodes=1
#SBATCH --time=8:00:00
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
# convert_fastaqual_fastq.py -c fastq_to_fastaqual -f Chip68.fastq -o fastaqual_Chip68

# split_libraries.py -m Chip67_map.txt -f fastaqual_Chip67/Chip67.fna -q fastaqual_Chip67/Chip67.qual -o split_output_Chip67 -b variable_length
# split_libraries.py -m Chip68_map.txt -f fastaqual_Chip68/Chip68.fna -q fastaqual_Chip68/Chip68.qual -o split_output_Chip68 -b variable_length

# split_libraries.py -m Chip67_map_5000.txt -f fastaqual_Chip67/Chip67.fna -q fastaqual_Chip67/Chip67.qual -o split_output_Chip67_5000 -b variable_length
# split_libraries.py -m Chip68_map_5000.txt -f fastaqual_Chip68/Chip68.fna -q fastaqual_Chip68/Chip68.qual -o split_output_Chip68_5000 -b variable_length

# cat split_output_Chip67_5000/seqs.fna split_output_Chip68_5000/seqs.fna Stoma_Combined_March2017.fna > Stoma_November_2018.fna


# pick_otus.py -i Stoma_November_2018.fna -o uclust_picked_otus
# pick_rep_set.py -i uclust_picked_otus/Stoma_November_2018_otus.txt -f Stoma_November_2018.fna -o rep_set.fna
assign_taxonomy.py -i rep_set.fna -r 97_otus.fasta -t 97_otu_taxonomy.txt  -o taxonomy_results/
make_otu_table.py -i uclust_picked_otus/Stoma_November_2018_otus.txt -t taxonomy_results/rep_set_tax_assignments.txt -o taxonomy_results/otu_table.biom

filter_taxa_from_otu_table.py -i taxonomy_results/otu_table.biom -o taxonomy_results/otu_table_curated.biom -n Unassigned
filter_otus_from_otu_table.py -i taxonomy_results/otu_table_curated.biom -o taxonomy_results/otu_table_curated_0.001.biom --min_count_fraction 0.001
filter_otus_from_otu_table.py -i taxonomy_results/otu_table_curated.biom -o taxonomy_results/otu_table_curated_0.0001.biom --min_count_fraction 0.0001
filter_otus_from_otu_table.py -i taxonomy_results/otu_table_curated.biom -o taxonomy_results/otu_table_curated_0.00001.biom --min_count_fraction 0.00001

summarize_taxa.py -i taxonomy_results/otu_table_curated_0.001.biom -L 2,3,4,5,6 -o taxonomy_L2-6_0.001/
summarize_taxa.py -i taxonomy_results/otu_table_curated_0.0001.biom -L 2,3,4,5,6 -o taxonomy_L2-6_0.0001/
summarize_taxa.py -i taxonomy_results/otu_table_curated_0.00001.biom -L 2,3,4,5,6 -o taxonomy_L2-6_0.00001/

# align_seqs.py -i rep_set.fna -o alignment/
# filter_alignment.py -i alignment/rep_set_aligned.fasta -o alignment/
# make_phylogeny.py -i alignment/rep_set_aligned_pfiltered.fasta -o rep_set_tree.tre
