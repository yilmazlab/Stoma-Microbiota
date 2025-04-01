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


setwd("~/Dropbox/Ongoing Analysis/StomaMerged_November_2018/Variant Analysis/Instrain Output")


# For E.coli pNpS_variants --------------------------------------------------------------

# Loading the Instrain output files
GUT_GENOME144544 <- read.xls("GUT_GENOME144544.xlsx", sheet = 1, header = TRUE)

SF2<-read.table(file = 'SF2.IS/output/SF2.IS_gene_info.tsv', sep = '\t', header = TRUE)
SF2$Sample <- "SF2"

SF3<-read.table(file = 'SF3.IS/output/SF3.IS_gene_info.tsv', sep = '\t', header = TRUE)
SF3$Sample <- "SF3"

SF12<-read.table(file = 'SF12.IS/output/SF12.IS_gene_info.tsv', sep = '\t', header = TRUE)
SF12$Sample <- "SF12"

SF14<-read.table(file = 'SF14.IS/output/SF14.IS_gene_info.tsv', sep = '\t', header = TRUE)
SF14$Sample <- "SF14"

SF17<-read.table(file = 'SF17.IS/output/SF17.IS_gene_info.tsv', sep = '\t', header = TRUE)
SF17$Sample <- "SF17"

SF20<-read.table(file = 'SF20.IS/output/SF20.IS_gene_info.tsv', sep = '\t', header = TRUE)
SF20$Sample <- "SF20"

SF26<-read.table(file = 'SF26.IS/output/SF26.IS_gene_info.tsv', sep = '\t', header = TRUE)
SF26$Sample <- "SF26"

SF21<-read.table(file = 'SF21.IS/output/SF21.IS_gene_info.tsv', sep = '\t', header = TRUE)
SF21$Sample <- "SF21"

SF24<-read.table(file = 'SF24.IS/output/SF24.IS_gene_info.tsv', sep = '\t', header = TRUE)
SF24$Sample <- "SF24"

SF25<-read.table(file = 'SF25.IS/output/SF25.IS_gene_info.tsv', sep = '\t', header = TRUE)
SF25$Sample <- "SF25"

# Combining them into one file
GUT_GENOME144544_all<-rbind(SF2,SF3,SF12,SF14,SF17,SF20,SF26,SF21,SF24,SF25)

# Subsetting the data for coverage >10 and breadth > 0.8
GUT_GENOME144544_all_subset <- GUT_GENOME144544_all  %>%  
                        filter(scaffold %in% c("GUT_GENOME144544_1", "GUT_GENOME144544_2")) %>%
                        filter(coverage>10,breadth>0.8,pNpS_variants>1)

unique_gene <- GUT_GENOME144544_all_subset %>% distinct(gene, .keep_all = TRUE)


SF2_149gene<-SF2[match(unique_gene$gene, SF2$gene), ]
SF3_149gene<-SF3[match(unique_gene$gene, SF3$gene), ]
SF12_149gene<-SF12[match(unique_gene$gene, SF12$gene), ]
SF17_149gene<-SF17[match(unique_gene$gene, SF17$gene), ]
SF14_149gene<-SF14[match(unique_gene$gene, SF14$gene), ]
SF20_149gene<-SF20[match(unique_gene$gene, SF20$gene), ]
SF21_149gene<-SF21[match(unique_gene$gene, SF21$gene), ]
SF24_149gene<-SF24[match(unique_gene$gene, SF24$gene), ]
SF25_149gene<-SF25[match(unique_gene$gene, SF25$gene), ]
SF26_149gene<-SF26[match(unique_gene$gene, SF26$gene), ]

# Combining them into one file
GUT_GENOME144544_149gene<-rbind(SF2_149gene,SF3_149gene,SF12_149gene,SF14_149gene,SF17_149gene,SF20_149gene,SF26_149gene,SF21_149gene,SF24_149gene,SF25_149gene)
write.csv(GUT_GENOME144544_149gene, "GUT_GENOME144544_149gene.csv")

count_gene<-GUT_GENOME144544_149gene %>% 
  count(gene)
count_gene <- count_gene  %>%  
  filter(n==10)

GUT_GENOME144544_73gene_final <- GUT_GENOME144544_149gene  %>%  
  filter(gene %in% count_gene$gene)
  
write.csv(GUT_GENOME144544_73gene_final, "GUT_GENOME144544_73gene_final.csv")

GUT_GENOME144544_73gene_final_list <- GUT_GENOME144544_73gene_final %>% select(gene,pNpS_variants,Sample)

GUT_GENOME144544_73gene_final_list.wide <- pivot_wider(GUT_GENOME144544_73gene_final_list, names_from = Sample, values_from = pNpS_variants)   
write.csv(GUT_GENOME144544_73gene_final_list.wide,"GUT_GENOME144544_73gene_final_list.wide.csv")

# We then organized in the excel the gene name annotated with SNPeff
# Then we used this file to generate the heatmap.

library(pheatmap)
library(RColorBrewer)
library(gdata)
library(tibble)
setwd("~/Dropbox/Ongoing Analysis/StomaMerged_November_2018/Variant Analysis/E.coli") # Location where you can find the text files for the generation of heatmap
annotation_table= read.table("annotation.txt")
Taxa<-read.table("data_GUT_GENOME144544.txt", sep = "\t", header = TRUE)
Taxa<- Taxa %>% remove_rownames %>% column_to_rownames(var="X")
# Taxa[is.na(Taxa)] <- as.double("NA")

ann_colors = list(
  TimePoint = c(Zero ="#fff5f0", Two="#fcbba1", Four="#fc9272", Six="#ef3b2c", Eight="#cb181d", Ten = "#a50f15", Fortyeight="#67000d"),
  Patient = c(P1="#B5BBE3", P2="#aaf2ff", P3="#3399FF",P4="#efb94c", P6="#cc6666", P7="#f4fcfc"))


rwbcols <- c( "#781d1a","#D33F6A","#E07B91", "#4A6FE3", "white")

pheatmap(Taxa, 
         fontsize_row = 10, fontsize_col = 10, cluster_cols = FALSE, cluster_rows = TRUE,
         color = colorRampPalette(rev(rwbcols))(100),
         cellwidth =30, cellheight = 11, # cutree_cols = 2, # cutree_rows = 3, 
         annotation_col=annotation_table, annotation_colors = ann_colors, border_color = "grey60", na_col = "grey90")

heatmap.2(Taxa, Rowv = F, Colv = F, trace = "none", na.color = "Green")

# For Klebsiella pneumaniae pNpS_variants --------------------------------------------------------------
setwd("~/Dropbox/Ongoing Analysis/StomaMerged_November_2018/Variant Analysis/Instrain Output") # Set library back where the Instrain Output files were

 # Loading the Instrain Output files
GUT_GENOME147598 <- read.xls("GUT_GENOME147598_2.xlsx", sheet = 1, header = TRUE)
SF2<-read.table(file = 'SF2.IS/output/SF2.IS_gene_info.tsv', sep = '\t', header = TRUE)
SF2$Sample <- "SF2"
SF3<-read.table(file = 'SF3.IS/output/SF3.IS_gene_info.tsv', sep = '\t', header = TRUE)
SF3$Sample <- "SF3"
SF9<-read.table(file = 'SF9.IS/output/SF9.IS_gene_info.tsv', sep = '\t', header = TRUE)
SF9$Sample <- "SF9"
SF10<-read.table(file = 'SF10.IS/output/SF10.IS_gene_info.tsv', sep = '\t', header = TRUE)
SF10$Sample <- "SF10"
SF20<-read.table(file = 'SF20.IS/output/SF20.IS_gene_info.tsv', sep = '\t', header = TRUE)
SF20$Sample <- "SF20"
SF21<-read.table(file = 'SF21.IS/output/SF21.IS_gene_info.tsv', sep = '\t', header = TRUE)
SF21$Sample <- "SF21"
SF23<-read.table(file = 'SF23.IS/output/SF23.IS_gene_info.tsv', sep = '\t', header = TRUE)
SF23$Sample <- "SF23"
SF25<-read.table(file = 'SF25.IS/output/SF25.IS_gene_info.tsv', sep = '\t', header = TRUE)
SF25$Sample <- "SF25"
SF26<-read.table(file = 'SF26.IS/output/SF26.IS_gene_info.tsv', sep = '\t', header = TRUE)
SF26$Sample <- "SF26"

# Combining the files into one file
GUT_GENOME147598_all<-rbind(SF2,SF3,SF9,SF10,SF20,SF21,SF23,SF26,SF21,SF25)

# Subsetting the data for coverage >10 and breadth > 0.8
GUT_GENOME147598_all_subset <- GUT_GENOME147598_all  %>%  
  filter(scaffold %in% c("GUT_GENOME147598_1", "GUT_GENOME147598_2", "GUT_GENOME147598_3", "GUT_GENOME147598_4")) %>%
  filter(coverage>10,breadth>0.8,pNpS_variants>1.0)

unique_gene <- GUT_GENOME147598_all_subset %>% distinct(gene, .keep_all = TRUE)


SF2_92gene<-SF2[match(unique_gene$gene, SF2$gene), ]
SF3_92gene<-SF3[match(unique_gene$gene, SF3$gene), ]
SF9_92gene<-SF9[match(unique_gene$gene, SF9$gene), ]
SF10_92gene<-SF10[match(unique_gene$gene, SF10$gene), ]
SF20_92gene<-SF20[match(unique_gene$gene, SF20$gene), ]
SF21_92gene<-SF21[match(unique_gene$gene, SF21$gene), ]
SF23_92gene<-SF23[match(unique_gene$gene, SF23$gene), ]
SF25_92gene<-SF25[match(unique_gene$gene, SF25$gene), ]
SF26_92gene<-SF26[match(unique_gene$gene, SF26$gene), ]

GUT_GENOME147598_92gene<-rbind(SF2_92gene,SF3_92gene,SF9_92gene,SF10_92gene,SF20_92gene,SF21_92gene,SF23_92gene,SF26_92gene,SF25_92gene)
write.csv(GUT_GENOME147598_92gene, "GUT_GENOME147598_92gene.csv")

count_gene<-GUT_GENOME147598_92gene %>% 
  count(gene)
count_gene <- count_gene  %>%  
  filter(n==9)

GUT_GENOME147598_83gene_final <- GUT_GENOME147598_92gene  %>%  
  filter(gene %in% count_gene$gene)

write.csv(GUT_GENOME147598_83gene_final, "GUT_GENOME147598_83gene_final.csv")

GUT_GENOME147598_83gene_final_list <- GUT_GENOME147598_83gene_final %>% select(gene,pNpS_variants,Sample) # dplyr again here

library(tidyr)
GUT_GENOME147598_83gene_final_list.wide <- pivot_wider(GUT_GENOME147598_83gene_final_list, names_from = Sample, values_from = pNpS_variants)   
write.csv(GUT_GENOME147598_83gene_final_list.wide,"GUT_GENOME147598_83gene_final_list.wide.csv")

# We then organized in the excel the gene name annotated with SNPeff
# Then we used this file to generate the heatmap.

library(pheatmap)
library(RColorBrewer)
library(gdata)

annotation_table= read.table("/Users/bahti/Dropbox/Ongoing\ Analysis/StomaMerged_November_2018/Variant\ Analysis/Klebsiella\ pneumoniae/annotation.txt")
Taxa<-read.table("/Users/bahti/Dropbox/Ongoing\ Analysis/StomaMerged_November_2018/Variant\ Analysis/Klebsiella\ pneumoniae/data_GUT_GENOME147598.txt", sep = "\t", header = TRUE)
Taxa<- Taxa %>% remove_rownames %>% column_to_rownames(var="X")

ann_colors = list(
  TimePoint = c(Zero ="#fff5f0", Two="#fcbba1", Four="#fc9272", Six="#ef3b2c", Eight="#cb181d", Ten = "#a50f15", Fortyeight="#67000d"),
  Patient = c(P1="#B5BBE3", P2="#aaf2ff", P3="#3399FF",P4="#efb94c", P6="#cc6666", P7="#f4fcfc"))


rwbcols <- c( "#781d1a","#E07B91", "#4A6FE3", "white")

pheatmap(Taxa, 
         fontsize_row = 10, fontsize_col = 10, cluster_cols = FALSE, cluster_rows = TRUE,
         color = colorRampPalette(rev(rwbcols))(100),
         cellwidth =30, cellheight = 11, #cutree_cols = 2, # cutree_rows = 3, 
         annotation_col=annotation_table, annotation_colors = ann_colors, border_color = "grey60")



# For E.coli SNPs --------------------------------------------------------------
setwd("~/Dropbox/Ongoing Analysis/StomaMerged_November_2018/Variant Analysis/Instrain Output") # Set library back where the Instrain Output files were
SF2<-read.table(file = 'SF2.IS/output/SF2.IS_SNVs.tsv', sep = '\t', header = TRUE)
SF2$Sample <- "SF2"

SF3<-read.table(file = 'SF3.IS/output/SF3.IS_SNVs.tsv', sep = '\t', header = TRUE)
SF3$Sample <- "SF3"

SF12<-read.table(file = 'SF12.IS/output/SF12.IS_SNVs.tsv', sep = '\t', header = TRUE)
SF12$Sample <- "SF12"

SF14<-read.table(file = 'SF14.IS/output/SF14.IS_SNVs.tsv', sep = '\t', header = TRUE)
SF14$Sample <- "SF14"

SF17<-read.table(file = 'SF17.IS/output/SF17.IS_SNVs.tsv', sep = '\t', header = TRUE)
SF17$Sample <- "SF17"

SF20<-read.table(file = 'SF20.IS/output/SF20.IS_SNVs.tsv', sep = '\t', header = TRUE)
SF20$Sample <- "SF20"

SF26<-read.table(file = 'SF26.IS/output/SF26.IS_SNVs.tsv', sep = '\t', header = TRUE)
SF26$Sample <- "SF26"

SF21<-read.table(file = 'SF21.IS/output/SF21.IS_SNVs.tsv', sep = '\t', header = TRUE)
SF21$Sample <- "SF21"

SF24<-read.table(file = 'SF24.IS/output/SF24.IS_SNVs.tsv', sep = '\t', header = TRUE)
SF24$Sample <- "SF24"

SF25<-read.table(file = 'SF25.IS/output/SF25.IS_SNVs.tsv', sep = '\t', header = TRUE)
SF25$Sample <- "SF25"

GUT_GENOME144544_all<-rbind(SF2,SF3,SF12,SF14,SF17,SF20,SF21,SF24,SF25,SF26)

GUT_GENOME144544_all_subset <- GUT_GENOME144544_all  %>%  
  filter(scaffold %in% c("GUT_GENOME144544_1", "GUT_GENOME144544_2"))

GUT_GENOME144544_all_subset <- GUT_GENOME144544_all_subset  %>%
  filter(position_coverage>10,allele_count>1,var_freq>0.2) %>%
  filter(mutation_type %in% c("N") )

GUT_GENOME144544_subset <- GUT_GENOME144544_all_subset %>% select(Sample,position,var_freq)

GUT_GENOME144544_wide <- pivot_wider(GUT_GENOME144544_subset, id_cols = position, names_from = Sample, values_from = var_freq)   
df <- apply(GUT_GENOME144544_wide,2,as.character)
write.csv(df,"GUT_GENOME144544_wide.csv")

# We then filtered and annotated the data in excel sheet.
# We included the SNPs that have value higher than zero in 7 samples out of 9 samples in this part of the analysis. We also filtered out the samples SNPS that have a value of 1 in more than 7 samples.
# Then we run the heatmap.

library(pheatmap)
library(RColorBrewer)
library(gdata)
library(tibble)

setwd("~/Dropbox/Ongoing Analysis/StomaMerged_November_2018/Variant Analysis/E.coli/E.coli/New/dNdS variant plot")

GUT_GENOME144544_vcf <- read.vcfR(file="SnpEff/GUT_GENOME144544_all_SNVs_snpeff.vcf", verbose = FALSE) # call vcf file
GUT_GENOME144544_vcf_tidy <- vcfR2tidy(GUT_GENOME144544_vcf)
GUT_GENOME144544_vcf_fix <- GUT_GENOME144544_vcf_tidy[["fix"]]

dNdS.E.coli<-read.table("~/Dropbox/Ongoing Analysis/StomaMerged_November_2018/Variant Analysis/E.coli/E.coli/New/dNdS variant plot/E.coli.dNdS.AF.txt", header = TRUE) # call sample file

dNdS.E.coli_vcf<-GUT_GENOME144544_vcf_fix %>% # subset based on patient data information
  filter(POS %in% dNdS.E.coli$position)

dNdS.E.coli<-read.table("~/Dropbox/Ongoing Analysis/StomaMerged_November_2018/Variant Analysis/E.coli/E.coli/New/dNdS variant plot/E.coli.dNdS.AF.txt", header = TRUE) # call sample file

dNdS.E.coli_vcf<-GUT_GENOME144544_vcf_fix %>% # subset based on patient data information
  filter(POS %in% dNdS.E.coli$position)

# We then prepared the annotation and SNVs allele frequency tables (called Taxa). We included non-synonmous and synonymous SNVs.

ann_colors = list(
  TimePoint = c(Zero ="#fff5f0", Two="#fcbba1", Four="#fc9272", Six="#ef3b2c", Eight="#cb181d", Ten = "#a50f15", Fortyeight="#67000d"),
  Patient = c(P1="#B5BBE3", P2="#aaf2ff", P3="#3399FF",P4="#efb94c", P6="#cc6666", P7="#f4fcfc"),
  SnpEff_Effect =c(INTRAGENIC="#ffffd9", DOWNSTREAM="#edf8b1" , UPSTREAM="#c7e9b4", STOP_GAINED="#7fcdbb", STOP_LOST="#41b6c4" , SYNONYMOUS_STOP="#1d91c0", SYNONYMOUS_CODING="#225ea8", NON_SYNONYMOUS_START="#253494", NON_SYNONYMOUS_CODING="#081d58"),
  SnpEff_Impact=c(LOW="#abf5be",MODERATE="#63c97d",HIGH="#24853c",MODIFIER="#023b10"),
  SnpEff_Functional_Class=c(NONE="white", MISSENSE="#8c705b",NONSENSE="#5e3e26",SILENT="#f5f2f0"))


rwbcols <- c( "#781d1a","#E07B91", "#4A6FE3", "white")

# Calling the files for the heatmap
annotation_table= read.table("/Users/bahti/Dropbox/Ongoing Analysis/StomaMerged_November_2018/Variant Analysis/E.coli/E.coli/New/dNdS variant plot/annotation.txt")
annotation_row= read.table("/Users/bahti/Dropbox/Ongoing Analysis/StomaMerged_November_2018/Variant Analysis/E.coli/E.coli/New/dNdS variant plot/annotation_row.2.synonymous.nonsynonymous.txt")
Taxa<-read.table("/Users/bahti/Dropbox/Ongoing Analysis/StomaMerged_November_2018/Variant Analysis/E.coli/E.coli/New/dNdS variant plot/E.coli.dNdS.AF.synonymous.nonsynonymous.txt", row.names=NULL)
row.names(Taxa) <- make.unique(as.character(Taxa$row.names))
Taxa$row.names <- NULL

dNdS.E.coli.plot <- pheatmap(Taxa, 
                             fontsize_row = 5, fontsize_col = 10, 
                             cluster_cols = FALSE, cluster_rows = TRUE,
                             #color = colorRampPalette(rev(rwbcols))(100),
                             cellwidth = 80, cellheight = 5, #cutree_rows = 6,  
                             annotation_col=annotation_table, annotation_row=annotation_row, annotation_colors = ann_colors, border_color =  "grey60")


# For Klebsiella pneumaniae --------------------------------------------------------------


SF2<-read.table(file = 'SF2.IS/output/SF2.IS_SNVs.tsv', sep = '\t', header = TRUE)
SF2$Sample <- "SF2"

SF3<-read.table(file = 'SF3.IS/output/SF3.IS_SNVs.tsv', sep = '\t', header = TRUE)
SF3$Sample <- "SF3"

SF9<-read.table(file = 'SF9.IS/output/SF9.IS_SNVs.tsv', sep = '\t', header = TRUE)
SF9$Sample <- "SF9"

SF10<-read.table(file = 'SF10.IS/output/SF10.IS_SNVs.tsv', sep = '\t', header = TRUE)
SF10$Sample <- "SF10"

SF20<-read.table(file = 'SF20.IS/output/SF20.IS_SNVs.tsv', sep = '\t', header = TRUE)
SF20$Sample <- "SF20"

SF21<-read.table(file = 'SF21.IS/output/SF21.IS_SNVs.tsv', sep = '\t', header = TRUE)
SF21$Sample <- "SF21"

SF23<-read.table(file = 'SF23.IS/output/SF23.IS_SNVs.tsv', sep = '\t', header = TRUE)
SF23$Sample <- "SF23"

SF25<-read.table(file = 'SF25.IS/output/SF25.IS_SNVs.tsv', sep = '\t', header = TRUE)
SF25$Sample <- "SF25"

SF26<-read.table(file = 'SF26.IS/output/SF26.IS_SNVs.tsv', sep = '\t', header = TRUE)
SF26$Sample <- "SF26"


GUT_GENOME147598_all<-rbind(SF2,SF3,SF9,SF10,SF20,SF21,SF23,SF25,SF26)

GUT_GENOME147598_all_subset <- GUT_GENOME147598_all  %>%  
  filter(scaffold %in% c("GUT_GENOME147598_1", "GUT_GENOME147598_2", "GUT_GENOME147598_3", "GUT_GENOME147598_4"))

GUT_GENOME147598_all_subset <- GUT_GENOME147598_all_subset %>%
  filter(position_coverage>10,allele_count>1,var_freq>0.2) %>%
  filter(mutation_type %in% c("N") )

GUT_GENOME147598_subset <- GUT_GENOME147598_all_subset %>% select(Sample,position,var_freq)

GUT_GENOME147598_wide <- pivot_wider(GUT_GENOME147598_subset, id_cols = position, names_from = Sample, values_from = var_freq)   
df <- apply(GUT_GENOME147598_wide,2,as.character)
write.csv(df,"GUT_GENOME147598_wide.csv")

# We then prepared the annotation and SNVs allele frequency tables (called Taxa). We included non-synonmous and synonymous SNVs.
setwd("/Users/bahti/Dropbox/Ongoing\ Analysis/StomaMerged_November_2018/Variant\ Analysis/Klebsiella\ pneumoniae/dNdS\ variant\ plot")
list.files()

GENOME147598_vcf <- read.vcfR(file="SnpEff/GUT_GENOME147598_all_SNVs_snpeff.vcf", verbose = FALSE) # call vcf file
GENOME147598_vcf_tidy <- vcfR2tidy(GENOME147598_vcf)
GENOME147598_vcf_fix <- GENOME147598_vcf_tidy[["fix"]]

dNdS.K.pneumoniae<-read.table("/Users/bahti/Dropbox/Ongoing\ Analysis/StomaMerged_November_2018/Variant\ Analysis/Klebsiella\ pneumoniae/dNdS\ variant\ plot/K.pneumoniae.dNdS.AF.txt", header=TRUE) # call sample file

dNdS.K.pneumoniae_vcf<-GENOME147598_vcf_fix %>% # subset based on patient data information
  filter(POS %in% dNdS.K.pneumoniae$position)

dNdS.K.pneumoniae_vcf_simple<-dNdS.K.pneumoniae_vcf[!duplicated(dNdS.K.pneumoniae_vcf$POS), ] # remove the duplicates
write.csv(dNdS.K.pneumoniae_vcf_simple, "/Users/bahti/Dropbox/Ongoing\ Analysis/StomaMerged_November_2018/Variant\ Analysis/Klebsiella\ pneumoniae/dNdS\ variant\ plot/dNdS.K.pneumoniae_vcf_simple_zeros.csv")

# We then filtered and annotated the data in excel sheet.
# We included the SNPs that have value higher than zero in 7 samples out of 9 samples in this part of the analysis. We also filtered out the samples SNPS that have a value of 1 in more than 7 samples.
# Then we run the heatmap.
# We then prepared the annotation and SNVs allele frequency tables (called Taxa). We included non-synonmous and synonymous SNVs.

library(pheatmap)
library(RColorBrewer)
library(gdata)
library(tibble)

# Calling the files for heatmap

annotation_table= read.table("/Users/bahti/Dropbox/Ongoing\ Analysis/StomaMerged_November_2018/Variant\ Analysis/Klebsiella\ pneumoniae/dNdS\ variant\ plot/annotation.txt")
annotation_row= read.table("/Users/bahti/Dropbox/Ongoing\ Analysis/StomaMerged_November_2018/Variant\ Analysis/Klebsiella\ pneumoniae/dNdS\ variant\ plot/annotation_row.2.synonymous.nonsynonmous.klebsiella.txt")
Taxa<-read.table("/Users/bahti/Dropbox/Ongoing\ Analysis/StomaMerged_November_2018/Variant\ Analysis/Klebsiella\ pneumoniae/dNdS\ variant\ plot/K.pneumoniae.dNdS.AF.synonymous.nonsynonmous.txt", row.names=NULL)
row.names(Taxa) <- make.unique(as.character(Taxa$row.names))
Taxa$row.names <- NULL

ann_colors = list(
  TimePoint = c(Zero ="#fff5f0", Two="#fcbba1", Four="#fc9272", Six="#ef3b2c", Eight="#cb181d", Ten = "#a50f15", Fortyeight="#67000d"),
  Patient = c(P1="#B5BBE3", P2="#aaf2ff", P3="#3399FF",P4="#efb94c", P6="#cc6666", P7="#f4fcfc"),
  SnpEff_Effect =c(INTRAGENIC="#ffffd9", DOWNSTREAM="#edf8b1" , UPSTREAM="#c7e9b4", STOP_GAINED="#7fcdbb", STOP_LOST="#41b6c4" , SYNONYMOUS_STOP="#1d91c0", SYNONYMOUS_CODING="#225ea8", NON_SYNONYMOUS_START="#253494", NON_SYNONYMOUS_CODING="#081d58"),
  SnpEff_Impact=c(LOW="#abf5be",MODERATE="#63c97d",HIGH="#24853c",MODIFIER="#023b10"),
  SnpEff_Functional_Class=c(NONE="white", MISSENSE="#8c705b",NONSENSE="#5e3e26",SILENT="#f5f2f0"))



rwbcols <- c( "#781d1a","#E07B91", "#4A6FE3", "white")


dNdS.K.pneumoniae.plot <- pheatmap(Taxa, 
                                   fontsize_row = 5, fontsize_col = 10, 
                                   cluster_cols = FALSE, cluster_rows = TRUE,
                                   #color = colorRampPalette(rev(rwbcols))(100),
                                   cellwidth = 80, cellheight = 6, #cutree_rows = 6,  
                                   annotation_col=annotation_table, annotation_row=annotation_row, annotation_colors = ann_colors, border_color =  "grey60")


# BY-P.copri --------------------------------------------------------------

BY1<-read.table(file = 'BY1.IS/output/BY1.IS_SNVs.tsv', sep = '\t', header = TRUE)
BY1$Sample <- "BY1"

BY2<-read.table(file = 'BY2.IS/output/BY2.IS_SNVs.tsv', sep = '\t', header = TRUE)
BY2$Sample <- "BY2"

BY3<-read.table(file = 'BY3.IS/output/BY3.IS_SNVs.tsv', sep = '\t', header = TRUE)
BY3$Sample <- "BY3"

BY4<-read.table(file = 'BY4.IS/output/BY4.IS_SNVs.tsv', sep = '\t', header = TRUE)
BY4$Sample <- "BY4"

PBC1<-read.table(file = 'PBC1.IS/output/PBC1.IS_SNVs.tsv', sep = '\t', header = TRUE)
PBC1$Sample <- "PBC1"

PBC2<-read.table(file = 'PBC2.IS/output/PBC2.IS_SNVs.tsv', sep = '\t', header = TRUE)
PBC2$Sample <- "PBC2"

PBC3<-read.table(file = 'PBC3.IS/output/PBC3.IS_SNVs.tsv', sep = '\t', header = TRUE)
PBC3$Sample <- "PBC3"

PBC4<-read.table(file = 'PBC4.IS/output/PBC4.IS_SNVs.tsv', sep = '\t', header = TRUE)
PBC4$Sample <- "PBC4"


combined_all<-rbind(BY1,BY2,BY3,BY4,PBC1,PBC2,PBC3,PBC4)
write.csv(combined_all,"combined_all.csv")

saveRDS(combined_all, "combined_all.rds")
combined_all = readRDS("combined_all.rds")

# P.copri subsetting
P.copri  <- combined_all %>% filter(grepl("GUT_GENOME140299", scaffold))
write.csv(P.copri,"P.copri.csv")


P.copri_subset <- P.copri  %>%
  filter(position_coverage>10,allele_count>1,var_freq>0.3) %>%
  filter(mutation_type %in% c("N") )

detach("package:dplyr", unload = TRUE)
library(dplyr)
P.copri_subset <- P.copri_subset %>% select(Sample,position,var_freq)

P.copri_subset_wide <- pivot_wider(P.copri_subset, id_cols = position, names_from = Sample, values_from = var_freq)   
df <- apply(P.copri_subset_wide,2,as.character)
write.csv(df,"P.copri_subset_wide.csv")

P.copri.SNP.Filtered <- read.csv("P.copri_subset_wide.csv")

library(pheatmap)
library(RColorBrewer)
library(gdata)
library(tibble)

annotation_table= read.table("/Users/bahti/Dropbox/Ongoing\ Analysis/StomaMerged_November_2018/Variant\ Analysis/Instrain\ Output/BY_PBC_annotation.txt")
Taxa<-read.table("P.copri_subset_wide.csv", sep = ",", header = TRUE)
Taxa<- Taxa %>% remove_rownames %>% column_to_rownames(var="position")


ann_colors = list(
  TimePoint = c(Zero ="#fff5f0", Two="#fcbba1", Four="#fc9272", Six="#ef3b2c"),
  Patient = c(BY="#B5BBE3", PBC="#3cc23c"))


rwbcols <- c( "#781d1a","#E07B91", "#4A6FE3", "white")

P.copri.heatmap<-pheatmap(Taxa, 
                          # fontsize_row = 5, fontsize_col = 10, 
                          cluster_cols = FALSE, cluster_rows = TRUE,
                          # color = colorRampPalette(rev(rwbcols))(100),
                          cellwidth = 100, cellheight = 0.2, cutree_cols = 2, # cutree_rows = 3, 
                          annotation_col=annotation_table, annotation_colors = ann_colors, border_color =  "grey60", useRaster = T)

ggsave(P.copri.heatmap, file="P.copri.heatmap.pdf", width=30, height=30, useDingbats=FALSE)

pheatmap(Taxa, 
         fontsize_row = 0.3, fontsize_col = 10, 
         cluster_cols = FALSE, cluster_rows = TRUE,
         color = colorRampPalette(rev(rwbcols))(100),
         cellwidth = 60, cellheight = 10, cutree_cols = 2, # cutree_rows = 3, 
         annotation_col=annotation_table, annotation_colors = ann_colors, border_color =  "grey60", useRaster = T)

showLegend = list(
  TimePoint = c(Zero ="#fff5f0", Two="#fcbba1", Four="#fc9272", Six="#ef3b2c"),
  Patient = c(BY="#B5BBE3", PBC="#3cc23c"))



# Make helper function to map metadata category to color
gAnnotationData = read.table("/Users/bahti/Dropbox/Ongoing\ Analysis/StomaMerged_November_2018/Variant\ Analysis/Instrain\ Output/BY_PBC_annotation.txt")

heatmap3(Taxa, showColDendro = T,  showRowDendro = T)





# R.bromii subsetting -----------------------------------------------------

R.bromii  <- combined_all %>% filter(grepl("GUT_GENOME143522", scaffold))
write.csv(R.bromii,"R.bromii.csv")


R.bromii_subset <- R.bromii  %>%
  filter(position_coverage>10,allele_count>1,var_freq>0.3) %>%
  filter(mutation_type %in% c("N") )

detach("package:dplyr", unload = TRUE)
library(dplyr)
R.bromii_subset <- R.bromii_subset %>% select(Sample,position,var_freq)

R.bromii_subset_wide <- pivot_wider(R.bromii_subset, id_cols = position, names_from = Sample, values_from = var_freq)   
df <- apply(R.bromii_subset_wide,2,as.character)
write.csv(df,"R.bromii_subset_wide.csv")

R.bromii.SNP.Filtered <- read.csv("R.bromii_subset_wide.csv")

library(pheatmap)
library(RColorBrewer)
library(gdata)
library(tibble)

annotation_table= read.table("/Users/bahti/Dropbox/Ongoing\ Analysis/StomaMerged_November_2018/Variant\ Analysis/Instrain\ Output/BY_PBC_annotation.txt")
Taxa<-read.table("R.bromii_subset_wide.csv", sep = ",", header = TRUE)
Taxa<- Taxa %>% remove_rownames %>% column_to_rownames(var="position")


ann_colors = list(
  TimePoint = c(Zero ="#fff5f0", Two="#fcbba1", Four="#fc9272", Six="#ef3b2c"),
  Patient = c(BY="#B5BBE3", PBC="#3cc23c"))


rwbcols <- c( "#781d1a","#E07B91", "#4A6FE3", "white")

R.bromii.heatmap<-pheatmap(Taxa, 
                           # fontsize_row = 5, fontsize_col = 10, 
                           cluster_cols = FALSE, cluster_rows = TRUE,
                           # color = colorRampPalette(rev(rwbcols))(100),
                           cellwidth = 100, cellheight = 0.2, cutree_cols = 2, # cutree_rows = 3, 
                           annotation_col=annotation_table, annotation_colors = ann_colors, border_color =  "grey60", useRaster = T)

ggsave(R.bromii.heatmap, file="R.bromii.heatmap.pdf", width=30, height=30, useDingbats=FALSE)



# Make helper function to map metadata category to color
gAnnotationData = read.table("/Users/bahti/Dropbox/Ongoing\ Analysis/StomaMerged_November_2018/Variant\ Analysis/Instrain\ Output/BY_PBC_annotation.txt")

heatmap3(Taxa, showColDendro = T,  showRowDendro = T)


# Faecalibacterium subsetting ---------------------------------------------

Faecalibacterium  <- combined_all %>% filter(grepl("GENOME042029", scaffold))
write.csv(Faecalibacterium,"Faecalibacterium.csv")

BY4_Faecalibacterium  <- BY4 %>% filter(grepl("GENOME042029", scaffold))
write.csv(BY4_Faecalibacterium,"BY4_Faecalibacterium.csv")

Faecalibacterium_subset <- Faecalibacterium  %>%
  filter(position_coverage>10,allele_count>1,var_freq>0.3) %>%
  filter(mutation_type %in% c("N") )

detach("package:dplyr", unload = TRUE)
library(dplyr)
Faecalibacterium_subset <- Faecalibacterium_subset %>% select(Sample,position,var_freq)

Faecalibacterium_subset_wide <- pivot_wider(Faecalibacterium_subset, id_cols = position, names_from = Sample, values_from = var_freq)   
df <- apply(Faecalibacterium_subset_wide,2,as.character)
write.csv(df,"Faecalibacterium_subset_wide.csv")

library(pheatmap)
library(RColorBrewer)
library(gdata)
library(tibble)

annotation_table= read.table("/Users/bahti/Dropbox/Ongoing\ Analysis/StomaMerged_November_2018/Variant\ Analysis/Instrain\ Output/BY_PBC_annotation.txt")
Taxa<-read.table("Faecalibacterium_subset_wide.csv", sep = ",", header = TRUE)
Taxa<- Taxa %>% remove_rownames %>% column_to_rownames(var="position")


ann_colors = list(
  TimePoint = c(Zero ="#fff5f0", Two="#fcbba1", Four="#fc9272", Six="#ef3b2c"),
  Patient = c(BY="#B5BBE3", PBC="#3cc23c"))


rwbcols <- c( "#781d1a","#E07B91", "#4A6FE3", "white")

Faecalibacterium.heatmap<-pheatmap(Taxa, 
                                   # fontsize_row = 5, fontsize_col = 10, 
                                   cluster_cols = FALSE, cluster_rows = TRUE,
                                   # color = colorRampPalette(rev(rwbcols))(100),
                                   cellwidth = 100, cellheight = 0.2, cutree_cols = 2, # cutree_rows = 3, 
                                   annotation_col=annotation_table, annotation_colors = ann_colors, border_color =  "grey60", useRaster = T)

ggsave(Faecalibacterium.heatmap, file="Faecalibacterium.heatmap.pdf", width=30, height=30, useDingbats=FALSE)



# Make helper function to map metadata category to color
gAnnotationData = read.table("/Users/bahti/Dropbox/Ongoing\ Analysis/StomaMerged_November_2018/Variant\ Analysis/Instrain\ Output/BY_PBC_annotation.txt")

heatmap3(Taxa, showColDendro = T,  showRowDendro = T)




# F.saccharivorans subsetting ---------------------------------------------

F.saccharivorans  <- combined_all %>% filter(grepl("GENOME001587", scaffold))
write.csv(F.saccharivorans,"F.saccharivorans.csv")

F.saccharivorans_subset <- F.saccharivorans  %>%
  filter(position_coverage>10,allele_count>1,var_freq>0.3) %>%
  filter(mutation_type %in% c("N") )

detach("package:dplyr", unload = TRUE)
library(dplyr)
F.saccharivorans_subset <- F.saccharivorans_subset %>% select(Sample,position,var_freq)

F.saccharivorans_subset_wide <- pivot_wider(F.saccharivorans_subset, id_cols = position, names_from = Sample, values_from = var_freq)   
df <- apply(F.saccharivorans_subset_wide,2,as.character)
write.csv(df,"F.saccharivorans_subset_wide.csv")

library(pheatmap)
library(RColorBrewer)
library(gdata)
library(tibble)

annotation_table= read.table("/Users/bahti/Dropbox/Ongoing\ Analysis/StomaMerged_November_2018/Variant\ Analysis/Instrain\ Output/BY_PBC_annotation.txt")
Taxa<-read.table("F.saccharivorans_subset_wide.csv", sep = ",", header = TRUE)
Taxa<- Taxa %>% remove_rownames %>% column_to_rownames(var="position")


ann_colors = list(
  TimePoint = c(Zero ="#fff5f0", Two="#fcbba1", Four="#fc9272", Six="#ef3b2c"),
  Patient = c(BY="#B5BBE3", PBC="#3cc23c"))


rwbcols <- c( "#781d1a","#E07B91", "#4A6FE3", "white")

F.saccharivorans.heatmap<-pheatmap(Taxa, 
                                   # fontsize_row = 5, fontsize_col = 10, 
                                   cluster_cols = FALSE, cluster_rows = TRUE,
                                   # color = colorRampPalette(rev(rwbcols))(100),
                                   cellwidth = 100, cellheight = 0.2, cutree_cols = 2, # cutree_rows = 3, 
                                   annotation_col=annotation_table, annotation_colors = ann_colors, border_color =  "grey60", useRaster = T)

ggsave(F.saccharivorans.heatmap, file="F.saccharivorans.heatmap.pdf", width=30, height=30, useDingbats=FALSE)



# Make helper function to map metadata category to color
gAnnotationData = read.table("/Users/bahti/Dropbox/Ongoing\ Analysis/StomaMerged_November_2018/Variant\ Analysis/Instrain\ Output/BY_PBC_annotation.txt")

heatmap3(Taxa, showColDendro = T,  showRowDendro = T)



# B.adolescentis subsetting
B.adolescentis  <- combined_all %>% filter(grepl("GENOME142522", scaffold))
write.csv(B.adolescentis,"B.adolescentis.csv")

B.adolescentis_subset <- B.adolescentis  %>%
  filter(position_coverage>10,allele_count>1,var_freq>0.3) %>%
  filter(mutation_type %in% c("N") )

detach("package:dplyr", unload = TRUE)
library(dplyr)
B.adolescentis_subset <- B.adolescentis_subset %>% select(Sample,position,var_freq)

B.adolescentis_subset_wide <- pivot_wider(B.adolescentis_subset, id_cols = position, names_from = Sample, values_from = var_freq)   
df <- apply(B.adolescentis_subset_wide,2,as.character)
write.csv(df,"B.adolescentis_subset_wide.csv")

library(pheatmap)
library(RColorBrewer)
library(gdata)
library(tibble)

annotation_table= read.table("/Users/bahti/Dropbox/Ongoing\ Analysis/StomaMerged_November_2018/Variant\ Analysis/Instrain\ Output/BY_PBC_annotation.txt")
Taxa<-read.table("/Users/bahti/Dropbox/Ongoing\ Analysis/StomaMerged_November_2018/Variant\ Analysis/Instrain\ Output/B.adolescentis_subset_wide.csv", sep = ",", header = TRUE)
Taxa<- Taxa %>% remove_rownames %>% column_to_rownames(var="position")


ann_colors = list(
  TimePoint = c(Zero ="#fff5f0", Two="#fcbba1", Four="#fc9272", Six="#ef3b2c"),
  Patient = c(BY="#B5BBE3", PBC="#3cc23c"))


rwbcols <- c( "#781d1a","#E07B91", "#4A6FE3", "white")

B.adolescentis.heatmap<-pheatmap(Taxa, 
                                 # fontsize_row = 5, fontsize_col = 10, 
                                 cluster_cols = FALSE, cluster_rows = TRUE,
                                 # color = colorRampPalette(rev(rwbcols))(100),
                                 cellwidth = 80, cellheight = 0.2, cutree_cols = 2, # cutree_rows = 3, 
                                 annotation_col=annotation_table, annotation_colors = ann_colors, border_color =  "grey60", useRaster = T)

ggsave(B.adolescentis.heatmap, file="B.adolescentis.heatmap2.pdf", width=30, height=30, useDingbats=FALSE)



# Make helper function to map metadata category to color
gAnnotationData = read.table("/Users/bahti/Dropbox/Ongoing\ Analysis/StomaMerged_November_2018/Variant\ Analysis/Instrain\ Output/BY_PBC_annotation.txt")

heatmap3(Taxa, showColDendro = T,  showRowDendro = T)



# Blautia subsetting ------------------------------------------------------


Blautia  <- combined_all %>% filter(grepl("GENOME052918", scaffold))
write.csv(Blautia,"Blautia.csv")

Blautia_subset <- Blautia  %>%
  filter(position_coverage>10,allele_count>1,var_freq>0.3) %>%
  filter(mutation_type %in% c("N") )

detach("package:dplyr", unload = TRUE)
library(dplyr)
Blautia_subset <- Blautia_subset %>% select(Sample,position,var_freq)

Blautia_subset_wide <- pivot_wider(Blautia_subset, id_cols = position, names_from = Sample, values_from = var_freq)   
df <- apply(Blautia_subset_wide,2,as.character)
write.csv(df,"Blautia_subset_wide.csv")

library(pheatmap)
library(RColorBrewer)
library(gdata)
library(tibble)

annotation_table= read.table("/Users/bahti/Dropbox/Ongoing\ Analysis/StomaMerged_November_2018/Variant\ Analysis/Instrain\ Output/BY_PBC_annotation.txt")
Taxa<-read.table("/Users/bahti/Dropbox/Ongoing\ Analysis/StomaMerged_November_2018/Variant\ Analysis/Instrain\ Output/Blautia_subset_wide.csv", sep = ",", header = TRUE)
Taxa<- Taxa %>% remove_rownames %>% column_to_rownames(var="position")


ann_colors = list(
  TimePoint = c(Zero ="#fff5f0", Two="#fcbba1", Four="#fc9272", Six="#ef3b2c"),
  Patient = c(BY="#B5BBE3", PBC="#3cc23c"))


rwbcols <- c( "#781d1a","#E07B91", "#4A6FE3", "white")

Blautia.heatmap<-pheatmap(Taxa, 
                          # fontsize_row = 5, fontsize_col = 10, 
                          cluster_cols = FALSE, cluster_rows = TRUE,
                          # color = colorRampPalette(rev(rwbcols))(100),
                          cellwidth = 80, cellheight = 0.2, cutree_cols = 2, # cutree_rows = 3, 
                          annotation_col=annotation_table, annotation_colors = ann_colors, border_color =  "grey60", useRaster = T)

ggsave(Blautia.heatmap, file="Blautia.heatmap.2.pdf", width=30, height=30, useDingbats=FALSE)



# Make helper function to map metadata category to color
gAnnotationData = read.table("/Users/bahti/Dropbox/Ongoing\ Analysis/StomaMerged_November_2018/Variant\ Analysis/Instrain\ Output/BY_PBC_annotation.txt")

heatmap3(Taxa, showColDendro = T,  showRowDendro = T)



