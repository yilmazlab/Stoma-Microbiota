# This script for generating the heatmap that shows metabolomic changes in ileal stoma samples from recovered CRC patients before and after they had eaten a standardized breakfast having fasted overnight. Unsupervised clustering of the annotated metabolome at different times after the meal and the inferred source of metabolites from the HMDB and food database is annotated adjacent to the y-axis. 185 metabolites are depicted on the heatmap based on their dynamics of increasing or decreasing from (a) after feeding at 4 out of 6 patients within hours. 
# This script doesn't only include the figure that we generated for the manuscript. It also helps to generate the supervised and non-supervised plot based on column and row information.

library(ape)
library(xlsx)
library(plyr)
library(dplyr)
library(gdata)
library(ggplot2)
library(ape)
library(extrafont)
library(scales)
library(MASS)
library(PMCMR)
library(psych)
library(pheatmap)
library(RColorBrewer)
library(phyloseq)
library(vegan)
library(ggthemes)
library(gridExtra)
library(ggrepel)
library(phangorn)
library(picante)
library(reshape2)
library(reshape)
library(reshape)
library(dada2)
library(ggpubr)
library(bioDist)
library(vegan)
library(mvtnorm)
library(reshape)
library(ggplot2)
library(igraph)
library(phangorn)
library(picante)
library(reshape2)
library(reshape)
library(fdrtool)
library(phyloseq)
library(qiime2R)
library(data.table)
library(MicrobeR)
library(plotly)
library(hablar)
library(ggpubr)
library(gridExtra)

setwd("~/Dropbox/Ongoing Analysis/StomaMerged_November_2018/Metabolomics/RawData/Z-score/HeatmapForSFig")

annotation_table = read.table("annotation.txt")
annotation_row = read.table("annotation_row.txt")
Taxa <- read.table("CRC.CD.metabolites.txt", row.names=NULL)
row.names(Taxa) <- make.unique(as.character(Taxa$row.names))
Taxa$row.names <- NULL


ann_colors = list(
  Location = c(Ileostoma="#B5BBE3", Colostoma_Right="#aaf2ff", Colostoma_Left="#3399FF"),
  Class = c(Benzene_and_substituted_derivatives = "#7fc97f", Benzenoids ="#beaed4", Carboxylic_acids_and_derivatives="#fdc086", Cinnamic_acids_and_derivatives="#ffff99", Organic_acids_and_derivatives="#386cb0", Organic_oxygen_compounds = "#f0027f", Organoheterocyclic_compounds ="#bf5b17", Organooxygen_compounds ="#666666"),
  Disease = c(CRC = "green", CD = "red"))

rwbcols <- c( "#9c0e09","#bdaf1e", "#04bd2f")


feedingmetabolites <- pheatmap(log2(Taxa+1),
                               fontsize_row = 6, fontsize_col = 5, 
                               cluster_cols = TRUE, cluster_rows =TRUE,
                               color =  brewer.pal(100, "Spectral"),
                               cellwidth = 8, cellheight = 6, # cutree_rows = 4,  
                               annotation_col=annotation_table, annotation_row=annotation_row, annotation_colors = ann_colors)
feedingmetabolites

feedingmetabolites <- pheatmap(log2(Taxa+1),
                               fontsize_row = 6, fontsize_col = 5, 
                               cluster_cols = TRUE, cluster_rows =FALSE,
                               color =  brewer.pal(100, "Spectral"),
                               cellwidth = 8, cellheight = 6, # cutree_rows = 4,  
                               annotation_col=annotation_table, annotation_row=annotation_row, annotation_colors = ann_colors)
feedingmetabolites

feedingmetabolites <- pheatmap(log2(Taxa+1),
                               fontsize_row = 6, fontsize_col = 5, 
                               cluster_cols = FALSE, cluster_rows =FALSE,
                               color =  brewer.pal(100, "Spectral"),
                               cellwidth = 8, cellheight = 6, # cutree_rows = 4,  
                               annotation_col=annotation_table, annotation_row=annotation_row, annotation_colors = ann_colors)
feedingmetabolites

feedingmetabolites <- pheatmap(log2(Taxa+1),
                               fontsize_row = 6, fontsize_col = 5, 
                               cluster_cols = FALSE, cluster_rows = TRUE,
                               color =  brewer.pal(100, "Spectral"),
                               cellwidth = 8, cellheight = 6, # cutree_rows = 4,  
                               annotation_col=annotation_table, annotation_row=annotation_row, annotation_colors = ann_colors)
feedingmetabolites

ggsave(feedingmetabolites, file="/Users/bahti/Dropbox/Ongoing\ Analysis/StomaMerged_November_2018/Metabolomics/bandit\ subset\ 2021/heatmap_increase_decrrease/feedingmetabolite.increasing.decreasing.V2.pdf", width=20, height=30, useDingbats=FALSE,limitsize = FALSE)


annotation_table = read.table("CRC/annotation.txt")
annotation_row = read.table("CRC/annotation_row.txt")
Taxa <- read.table("CRC/CRC.metabolites.txt", row.names=NULL)
row.names(Taxa) <- make.unique(as.character(Taxa$row.names))
Taxa$row.names <- NULL

ann_colors = list(
  Class = c(Benzene_and_substituted_derivatives = "#7fc97f", Benzenoids ="#beaed4", Carboxylic_acids_and_derivatives="#fdc086", Cinnamic_acids_and_derivatives="#ffff99", Organic_acids_and_derivatives="#386cb0", Organic_oxygen_compounds = "#f0027f", Organoheterocyclic_compounds ="#bf5b17", Organooxygen_compounds ="#666666"),
  Location = c(Ileostoma="#a1d76a", Colostoma_Right="#f7f7f7", Colostoma_Left="#e9a3c9"))

feedingmetabolites <- pheatmap(log2(Taxa+1),
                               fontsize_row = 6, fontsize_col = 5, 
                               cluster_cols =TRUE, cluster_rows =TRUE,
                               color =  brewer.pal(100, "Spectral"),
                               cellwidth = 8, cellheight = 6, # cutree_rows = 4,  
                               annotation_col=annotation_table, annotation_row=annotation_row, annotation_colors = ann_colors)
feedingmetabolites
ggsave(feedingmetabolites, file="CRC/LocationalChanges.pdf", width=20, height=30, useDingbats=FALSE,limitsize = FALSE)

annotation_table = read.table("CD/annotation.txt")
annotation_row = read.table("CD/annotation_row.txt")
Taxa <- read.table("CD/CD.metabolites.txt", row.names=NULL)
row.names(Taxa) <- make.unique(as.character(Taxa$row.names))
Taxa$row.names <- NULL

ann_colors = list(
  Class = c(Benzene_and_substituted_derivatives = "#7fc97f", Benzenoids ="#beaed4", Carboxylic_acids_and_derivatives="#fdc086", Cinnamic_acids_and_derivatives="#ffff99", Organic_acids_and_derivatives="#386cb0", Organic_oxygen_compounds = "#f0027f", Organoheterocyclic_compounds ="#bf5b17", Organooxygen_compounds ="#666666"),
  Location = c(Ileostoma="#a1d76a", Colostoma_Right="#f7f7f7", Colostoma_Left="#e9a3c9"))

feedingmetabolites <- pheatmap(log2(Taxa+1),
                               fontsize_row = 6, fontsize_col = 5, 
                               cluster_cols = TRUE, cluster_rows =TRUE,
                               color =  brewer.pal(100, "Spectral"),
                               cellwidth = 8, cellheight = 6, # cutree_rows = 4,  
                               annotation_col=annotation_table, annotation_row=annotation_row, annotation_colors = ann_colors)
feedingmetabolites
ggsave(feedingmetabolites, file="CD/LocationalChanges.pdf", width=20, height=30, useDingbats=FALSE,limitsize = FALSE)



