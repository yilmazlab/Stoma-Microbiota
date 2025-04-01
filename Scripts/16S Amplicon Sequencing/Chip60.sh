#!/bin/sh
#!/bin/bash
#SBATCH --mail-user=leven.terziev@dkf.unibe.ch
#SBATCH --mail-type=end
#SBATCH --job-name="Chip60"
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
cp Chip60.fastq $TMPDIR/file.fastq
cp Chip_60_MappingFile_3000.txt $TMPDIR/map.txt
mkdir ./Results_$SLURM_JOB_ID/

## Start of QIIME pipeline ##
## Convert fastq file into .fna and .qual files ##
convert_fastaqual_fastq.py -c fastq_to_fastaqual -f $TMPDIR/file.fastq -o $TMPDIR/fastaqual
cp -R $TMPDIR/fastaqual ./Results_$JOB_ID

## Split libraries based on barcodes provided in map.txt ##
split_libraries.py -m $TMPDIR/map.txt -f $TMPDIR/fastaqual/file.fna -q $TMPDIR/fastaqual/file.qual -o $TMPDIR/split_output/ -b variable_length
cp -R $TMPDIR/split_output ./Results_$JOB_ID

## Pick Operational Taxonomic Units and assign taxonomies using Greengenes. The output folder is defined with -o rather than using the default output which would be in the HOME directory ##
#Method for picking OTUs. Valid choices are: sortmerna, mothur, trie, uclust_ref, usearch, usearch_ref, blast, usearch61, usearch61_ref, sumaclust, swarm, prefix_suffix, cdhit, uclust. The mothur method requires an input file of aligned sequences. usearch will enable the usearch quality filtering pipeline. [default: uclust]
pick_otus.py -i $TMPDIR/split_output/seqs.fna -o $TMPDIR/uclust_picked_otus
cp -R $TMPDIR/uclust_picked_otus ./Results_$JOB_ID

pick_rep_set.py -i $TMPDIR/uclust_picked_otus/seqs_otus.txt -f $TMPDIR/split_output/seqs.fna -o $TMPDIR/rep_set.fna
cp $TMPDIR/rep_set.fna ./Results_$JOB_ID/rep_set.fna

## Assign taxonomies using GreenGenes ##
#Taxon assignment method -m, must be one of rdp, blast, rtax, mothur, uclust, sortmerna [default: uclust]
assign_taxonomy.py -i $TMPDIR/rep_set.fna -r /db/SOFTWARE/qiime/1.8.0/gg_12_10_otus/rep_set/97_otus.fasta -t /db/SOFTWARE/qiime/1.8.0/gg_12_10_otus/taxonomy/97_otu_taxonomy.txt -o $TMPDIR/taxonomy_results/
cp -R $TMPDIR/taxonomy_results ./Results_$JOB_ID

## Plot the taxonomy data ##
make_otu_table.py -i $TMPDIR/uclust_picked_otus/seqs_otus.txt -t $TMPDIR/taxonomy_results/rep_set_tax_assignments.txt -o $TMPDIR/taxonomy_results/otu_table.biom
cp $TMPDIR/taxonomy_results/otu_table.biom ./Results_$JOB_ID/otu_table.biom

## exlcude specific samples or samples with less than a given number of reads from the otu_table ##
#filter_samples_from_otu_table.py -i otu_table.biom -o otu_table_new.biom --sample_id_fp map.txt (Map.txt contians the samples to be kept!)
#filter_samples_from_otu_table.py -i otu_table.biom -o otu_table_no_low_coverage_samples.biom -n 500 (at least 500 sequences to be kept)

## sort otu_table.biom based on Description +/- chimeric sequences removed ##
sort_otu_table.py -i $TMPDIR/taxonomy_results/otu_table.biom -o $TMPDIR/taxonomy_results/sorted_otu_table.biom -m $TMPDIR/map.txt -s Description

## Remove samples based on Description ##
#filter_samples_from_otu_table.py -i $TMPDIR/taxonomy_results/sorted_otu_table.biom -o $TMPDIR/taxonomy_results/sorted_otu_table_removed.biom -m $TMPDIR/map.txt -s 'Description:*,!remove'

## Remove Unassigend OTU's from sorted_otu_table.biom ##
filter_taxa_from_otu_table.py -i $TMPDIR/taxonomy_results/sorted_otu_table.biom -o $TMPDIR/taxonomy_results/otu_table_curated.biom -n Unassigned
cp $TMPDIR/taxonomy_results/otu_table_curated.biom ./Results_$JOB_ID/otu_table_curated.biom

## Generate plots inlcusive unassigned reads ##
summarize_taxa.py -i $TMPDIR/taxonomy_results/sorted_otu_table.biom -L 2,3,4,5,6,7 -o $TMPDIR/taxonomy_summaries_incl_unassigned/
plot_taxa_summary.py -i $TMPDIR/taxonomy_summaries_incl_unassigned/sorted_otu_table_L2.txt,$TMPDIR/taxonomy_summaries_incl_unassigned/sorted_otu_table_L3.txt,$TMPDIR/taxonomy_summaries_incl_unassigned/sorted_otu_table_L4.txt,$TMPDIR/taxonomy_summaries_incl_unassigned/sorted_otu_table_L5.txt,$TMPDIR/taxonomy_summaries_incl_unassigned/sorted_otu_table_L6.txt,$TMPDIR/taxonomy_summaries_incl_unassigned/sorted_otu_table_L7.txt -o $TMPDIR/taxonomy_plots_2-7_incl_unassigned/ -c bar,pie,area
cp -R $TMPDIR/taxonomy_plots_2-7_incl_unassigned ./Results_$JOB_ID

## Generate plots with unassigned reads removed ##
summarize_taxa.py -i $TMPDIR/taxonomy_results/otu_table_curated.biom -L 2,3,4,5,6,7 -o $TMPDIR/taxonomy_summaries/
plot_taxa_summary.py -i $TMPDIR/taxonomy_summaries/otu_table_curated_L2.txt,$TMPDIR/taxonomy_summaries/otu_table_curated_L3.txt,$TMPDIR/taxonomy_summaries/otu_table_curated_L4.txt,$TMPDIR/taxonomy_summaries/otu_table_curated_L5.txt,$TMPDIR/taxonomy_summaries/otu_table_curated_L6.txt,$TMPDIR/taxonomy_summaries/otu_table_curated_L7.txt -o $TMPDIR/taxonomy_plots_2-7/ -c bar,pie,area
cp -R $TMPDIR/taxonomy_plots_2-7 ./Results_$JOB_ID

## End of qiime taxonomy pipeline ##

## Start of Alpha Diversity ##
alpha_diversity.py -i $TMPDIR/taxonomy_results/otu_table_curated.biom -o $TMPDIR/alpha_diversity.txt -m shannon,simpson,observed_species,chao1
cp $TMPDIR/alpha_diversity.txt ./Results_$JOB_ID/alpha_diversity.txt

## Beta Diversity Measures ##
## Start of PCoA analysis ##
align_seqs.py -i $TMPDIR/rep_set.fna -o $TMPDIR/alignment/

filter_alignment.py -i $TMPDIR/alignment/rep_set_aligned.fasta -o $TMPDIR/alignment/

cp -R $TMPDIR/alignment ./Results_$JOB_ID

make_phylogeny.py -i $TMPDIR/alignment/rep_set_aligned_pfiltered.fasta -o $TMPDIR/rep_set_tree.tre

cp $TMPDIR/rep_set_tree.tre ./Results_$JOB_ID/rep_set_tree.tre

## Beta diversity: inlcude parameter -e for numer of sequences to be sampled to normalize number of sequences? ##
beta_diversity_through_plots.py -i $TMPDIR/taxonomy_results/otu_table_curated.biom -m $TMPDIR/map.txt -o $TMPDIR/PCoA_folder -t $TMPDIR/rep_set_tree.tre

cp -R $TMPDIR/PCoA_folder ./Results_$JOB_ID

## Make weighted 2D plots ##
make_2d_plots.py -i $TMPDIR/PCoA_folder/weighted_unifrac_pc.txt -m $TMPDIR/map.txt -b 'Disease' -o $TMPDIR/2d_plots_weighted/
cp -R $TMPDIR/2d_plots_weighted ./Results_$JOB_ID

compare_categories.py --method anosim -i $TMPDIR/PCoA_folder/weighted_unifrac_dm.txt -m $TMPDIR/map.txt -c Description -o $TMPDIR/anosim_out_weighted -n 9999
cp -R $TMPDIR/anosim_out_weighted ./Results_$JOB_ID

## Make unweighted 2D plots ##
make_2d_plots.py -i $TMPDIR/PCoA_folder/unweighted_unifrac_pc.txt -m $TMPDIR/map.txt -b 'Disease' -o $TMPDIR/2d_plots_unweighted/
cp -R $TMPDIR/2d_plots_unweighted ./Results_$JOB_ID

compare_categories.py --method anosim -i $TMPDIR/PCoA_folder/unweighted_unifrac_dm.txt -m $TMPDIR/map.txt -c Description -o $TMPDIR/anosim_out_unweighted -n 9999
cp -R $TMPDIR/anosim_out_unweighted ./Results_$JOB_ID


## End of script ##
