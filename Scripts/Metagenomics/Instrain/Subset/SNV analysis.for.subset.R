library(ape)
library(plyr)
library(dplyr)
library(gdata)
library(ggplot2)
library(ape)
library(extrafont)
library(scales)
library(MASS)
library(ggpmisc)
library(PMCMR)
library(psych)
library(pheatmap)
library(RColorBrewer)
library(phyloseq)
library(varhandle)
library(vegan)
library(ggthemes)
library(gridExtra)
library(ggrepel)
library(phangorn)
library(picante)
library(reshape2)
library(reshape)
library(reshape)
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
library(icesTAF)
library(data.table)
library(ggsankey)
library(xlsx)
library(openxlsx)
library(tidyr)


setwd("~/Dropbox/Ongoing Analysis/StomaMerged_November_2018/Variant Analysis/Subsetting for E.coli and Klebsiella")


# For E.coli SNPs --------------------------------------------------------------

SF2<-read.table(file = 'E.coli_sub/E.coli_SF2.IS/output/E.coli_SF2.IS_SNVs.tsv', sep = '\t', header = TRUE)
SF2$Sample <- "SF2"

SF3<-read.table(file = 'E.coli_sub/E.coli_SF3.IS/output/E.coli_SF3.IS_SNVs.tsv', sep = '\t', header = TRUE)
SF3$Sample <- "SF3"

SF12<-read.table(file = 'E.coli_sub/E.coli_SF12.IS/output/E.coli_SF12.IS_SNVs.tsv', sep = '\t', header = TRUE)
SF12$Sample <- "SF12"

SF14<-read.table(file = 'E.coli_sub/E.coli_SF14.IS/output/E.coli_SF14.IS_SNVs.tsv', sep = '\t', header = TRUE)
SF14$Sample <- "SF14"

SF17<-read.table(file = 'E.coli_sub/E.coli_SF17.IS/output/E.coli_SF17.IS_SNVs.tsv', sep = '\t', header = TRUE)
SF17$Sample <- "SF17"

SF20<-read.table(file = 'E.coli_sub/E.coli_SF20.IS/output/E.coli_SF20.IS_SNVs.tsv', sep = '\t', header = TRUE)
SF20$Sample <- "SF20"

SF26<-read.table(file = 'E.coli_sub/E.coli_SF26.IS/output/E.coli_SF26.IS_SNVs.tsv', sep = '\t', header = TRUE)
SF26$Sample <- "SF26"

SF21<-read.table(file = 'E.coli_sub/E.coli_SF21.IS/output/E.coli_SF21.IS_SNVs.tsv', sep = '\t', header = TRUE)
SF21$Sample <- "SF21"

SF24<-read.table(file = 'E.coli_sub/E.coli_SF24.IS/output/E.coli_SF24.IS_SNVs.tsv', sep = '\t', header = TRUE)
SF24$Sample <- "SF24"

SF25<-read.table(file = 'E.coli_sub/E.coli_SF25.IS/output/E.coli_SF25.IS_SNVs.tsv', sep = '\t', header = TRUE)
SF25$Sample <- "SF25"

# Somehow you need to call dplyr library again here if you have any problem.
GUT_GENOME144544_all<-rbind(SF2,SF3,SF12,SF14,SF17,SF20,SF21,SF24,SF25,SF26)

# Subset for E.coli genome
GUT_GENOME144544_all_subset <- GUT_GENOME144544_all  %>%  
  filter(scaffold %in% c("GUT_GENOME144544_1", "GUT_GENOME144544_2"))

# Filtering the data for allele count more than one and variant frequency higher than 20%
GUT_GENOME144544_all_subset <- GUT_GENOME144544_all_subset  %>%
  filter(position_coverage>10, allele_count>1,var_freq>0.2)

GUT_GENOME144544_subset <- GUT_GENOME144544_all_subset %>% select(Sample,position,var_freq)

GUT_GENOME144544_wide <- pivot_wider(GUT_GENOME144544_subset, id_cols = position, names_from = Sample, values_from = var_freq)   
df <- apply(GUT_GENOME144544_wide,2,as.character)

write.csv(GUT_GENOME144544_all_subset_2,"GUT_GENOME144544_all_subset.csv")

# Data further subsetted using excel document. I didn't use in this part the R. 
# I wanteded to look at the data in details and I did filtering  in excel.

library(pheatmap)
library(RColorBrewer)
library(gdata)
library(tibble)

annotation_table= read.table("E.coli_sub/Heatmap/annotation.2.txt") # Contains the patient ID and sampling time
Taxa<-read.table("E.coli_sub/Heatmap/E.coli.dNdS.AF.synonymous.nonsynonymous.combined.txt", header = TRUE) # This file contains combinedd filtered and non-filtered original nonsynonymous and synonomous allelic frequency data of the manuscript. 
annotation_row= read.table("/Users/bahti/Dropbox/Ongoing Analysis/StomaMerged_November_2018/Variant Analysis/E.coli/E.coli/New/dNdS variant plot/annotation_row.2.synonymous.nonsynonymous.txt") # It contains the SNPEff information.



ann_colors = list(
  TimePoint = c(Zero ="#fff5f0", Two="#fcbba1", Four="#fc9272", Six="#ef3b2c", Eight="#cb181d", Ten = "#a50f15", Fortyeight="#67000d"),
  Patient = c(P1="#B5BBE3", P2="#aaf2ff", P3="#3399FF",P4="#efb94c", P6="#cc6666", P7="#f4fcfc"))

rwbcols <- c( "#781d1a","#E07B91", "#4A6FE3", "white")


E.coli.syn.nonsyn.combined.plot <- pheatmap(Taxa, 
                             fontsize_row = 5, fontsize_col = 14, 
                             cluster_cols = TRUE, cluster_rows = TRUE,
                             # color = colorRampPalette(rev(rwbcols))(100),
                             cellwidth = 40, cellheight = 5, #cutree_rows = 6,  
                             annotation_col=annotation_table, annotation_row=annotation_row, annotation_colors = ann_colors, border_color =  "grey60")

ggsave(E.coli.syn.nonsyn.combined.plot, file="E.coli_sub/Heatmap/E.coli.syn.nonsyn.combined.plot.for.subsetvsoriginaldata.pdf", width=25, height=50, useDingbats=FALSE,limitsize = FALSE)


# For Klebsiella pneumaniae --------------------------------------------------------------


SF2<-read.table(file = 'Klebsiella_sub/Klebsiella_SF2.IS/output/Klebsiella_SF2.IS_SNVs.tsv', sep = '\t', header = TRUE)
SF2$Sample <- "SF2"

SF3<-read.table(file = 'Klebsiella_sub/Klebsiella_SF3.IS/output/Klebsiella_SF3.IS_SNVs.tsv', sep = '\t', header = TRUE)
SF3$Sample <- "SF3"

SF9<-read.table(file = 'Klebsiella_sub/Klebsiella_SF9.IS/output/Klebsiella_SF9.IS_SNVs.tsv', sep = '\t', header = TRUE)
SF9$Sample <- "SF9"

SF10<-read.table(file = 'Klebsiella_sub/Klebsiella_SF10.IS/output/Klebsiella_SF10.IS_SNVs.tsv', sep = '\t', header = TRUE)
SF10$Sample <- "SF10"

SF20<-read.table(file = 'Klebsiella_sub/Klebsiella_SF20.IS/output/Klebsiella_SF20.IS_SNVs.tsv', sep = '\t', header = TRUE)
SF20$Sample <- "SF20"

SF21<-read.table(file = 'Klebsiella_sub/Klebsiella_SF21.IS/output/Klebsiella_SF21.IS_SNVs.tsv', sep = '\t', header = TRUE)
SF21$Sample <- "SF21"

SF23<-read.table(file = 'Klebsiella_sub/Klebsiella_SF23.IS/output/Klebsiella_SF23.IS_SNVs.tsv', sep = '\t', header = TRUE)
SF23$Sample <- "SF23"

SF25<-read.table(file = 'Klebsiella_sub/Klebsiella_SF25.IS/output/Klebsiella_SF25.IS_SNVs.tsv', sep = '\t', header = TRUE)
SF25$Sample <- "SF25"

SF26<-read.table(file = 'Klebsiella_sub/Klebsiella_SF26.IS/output/Klebsiella_SF26.IS_SNVs.tsv', sep = '\t', header = TRUE)
SF26$Sample <- "SF26"


GUT_GENOME147598_all<-rbind(SF2,SF3,SF9,SF10,SF20,SF21,SF23,SF25,SF26)

GUT_GENOME147598_all_subset <- GUT_GENOME147598_all  %>%  
  filter(scaffold %in% c("GUT_GENOME147598_1","GUT_GENOME147598_2","GUT_GENOME147598_3","GUT_GENOME147598_4"))

GUT_GENOME147598_all_subset_2 <- GUT_GENOME147598_all_subset %>%
  filter(position_coverage>10,allele_count>1,var_freq>0.2) # If I keep the original criteria we lost around 100 data point. Problem is the allele count. We have some with value of 1.

GUT_GENOME147598_all_subset_2 <- GUT_GENOME147598_all_subset %>% select(Sample,position,var_freq)

GUT_GENOME147598_wide <- pivot_wider(GUT_GENOME147598_all_subset_2, id_cols = position, names_from = Sample, values_from = var_freq)   
df <- apply(GUT_GENOME147598_wide,2,as.character)


write.csv(df,"Klebsiella_sub/Analysis/GUT_GENOME144544_wide_all_subset.csv")

# Data further subsetted using excel document. I didn't use in this part the R. 
# I wanteded to look at the data in details and I did filtering  in excel.


library(pheatmap)
library(RColorBrewer)
library(gdata)
library(tibble)

annotation_table= read.table("Klebsiella_sub/Heatmap/annotation.2.txt") # Contains the patient ID and sampling time
Taxa<-read.table("Klebsiella_sub/Heatmap/Klebsiella.dNdS.AF.synonymous.nonsynonymous.combined.txt", header = TRUE) # This file contains combinedd filtered and non-filtered original nonsynonymous and synonomous allelic frequency data of the manuscript. 
annotation_row= read.table("Klebsiella_sub/Heatmap/annotation_row.2.synonymous.nonsynonmous.txt") # It contains the SNPEff information.



ann_colors = list(
  TimePoint = c(Zero ="#fff5f0", Two="#fcbba1", Four="#fc9272", Six="#ef3b2c", Eight="#cb181d", Ten = "#a50f15", Fortyeight="#67000d"),
  Patient = c(P1="#B5BBE3", P2="#aaf2ff", P3="#3399FF",P4="#efb94c", P6="#cc6666", P7="#f4fcfc"))

rwbcols <- c( "#781d1a","#E07B91", "#4A6FE3", "white")


Klebsiella.syn.nonsyn.combined.plot <- pheatmap(Taxa, 
                                            fontsize_row = 5, fontsize_col = 14, 
                                            cluster_cols = TRUE, cluster_rows = TRUE,
                                            # color = colorRampPalette(rev(rwbcols))(100),
                                            cellwidth = 40, cellheight = 5, #cutree_rows = 6,  
                                            annotation_col=annotation_table, annotation_row=annotation_row, annotation_colors = ann_colors, border_color =  "grey60")

ggsave(Klebsiella.syn.nonsyn.combined.plot, file="Klebsiella_sub/Heatmap/Klebsiella.syn.nonsyn.combined.plot.for.subsetvsoriginaldata.pdf", width=25, height=50, useDingbats=FALSE,limitsize = FALSE)


pheatmap(Taxa, 
         fontsize_row = 10, fontsize_col = 10, 
         cluster_cols = FALSE, cluster_rows = TRUE,
         color = colorRampPalette(rev(rwbcols))(100),
         cellwidth = 60, cellheight = 10, #cutree_cols = 2, # cutree_rows = 3, 
         annotation_col=annotation_table, annotation_colors = ann_colors, border_color =  "grey60", useRaster = T)




