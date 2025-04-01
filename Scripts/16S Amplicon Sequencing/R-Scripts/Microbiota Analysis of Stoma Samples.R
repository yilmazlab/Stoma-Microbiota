# StomaMerged_January_2020

# Phyloseq Analysis
# Loading libraries
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
library(btools)
library(pairwiseAdonis)

setwd("~/Dropbox/Ongoing Analysis/StomaMerged_November_2018/January2020/QIIME")

source("miseqR.R")
# Loading QIIME Files -----------------------------------------------------
My_Theme = theme(
  axis.title.x = element_text(size = 18),
  axis.text.x = element_text(size = 18),
  axis.title.y = element_text(size = 18),
  axis.text.y = element_text(size = 18))

# Loading necessary files for phyloseq #
StomaMerged = readRDS("StomaMerged.rds") # Read rds

StomaMerged
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 4028 taxa and 373 samples ]
# sample_data() Sample Data:       [ 373 samples by 39 sample variables ]
# tax_table()   Taxonomy Table:    [ 4028 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 4028 tips and 3551 internal nodes ]


# General Figures for Demographic Info ------------------------------------
# Checking the sequencing depth
sdt = data.table(as(sample_data(StomaMerged), "data.frame"),
                 TotalReads = sample_sums(StomaMerged), keep.rownames = TRUE)
setnames(sdt, "rn", "SampleID")
pSeqDepth = ggplot(sdt, aes(TotalReads)) + geom_histogram() + ggtitle("Sequencing Depth")
pSeqDepth
write.csv(sdt, "map_with_totalreads.csv")

My_Theme = theme(
  axis.title.x = element_text(size = 18),
  axis.text.x = element_text(size = 18),
  axis.title.y = element_text(size = 18),
  axis.text.y = element_text(size = 18))

# Plotting Sequencing Depth vs. Disease
sdt = data.table(as(sample_data(StomaMerged), "data.frame"), TotalReads = sample_sums(StomaMerged), keep.rownames = TRUE)
setnames(sdt, "rn", "SampleID")
pSeqDepth = ggplot(sdt, aes(TotalReads)) + geom_histogram() + ggtitle("Sequencing Depth")
pSeqDepth

sdt$Diagnosis_f = factor(sdt$Diagnosis_v2, levels=c('CRC','IBD','Other_cured','Other_chronic'))

SequencingDepthvsDisease <- ggplot (sdt, aes(SampleID, TotalReads, color = Diagnosis_v2)) + 
  geom_point(size = 2) +  geom_smooth(method = lm) +
  ggtitle("Sequencing Depth vs. Disease") + scale_y_log10() +
  theme(axis.line = element_line(colour = "black"), axis.title=element_text(face="bold", size="50", color="black"),
        axis.text.x = element_text(colour="grey20",size=1,angle=45,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=25,angle=0,hjust=1,vjust=0,face="plain"),
        legend.text=element_text(face="bold", size="50", color="black"),
        legend.title=element_text(face="bold", size="50", color="black")) + theme_bw() + theme(plot.title = element_text(hjust=0.5)) + 
  facet_grid(Diagnosis_f ~ .)
SequencingDepthvsDisease
ggsave(SequencingDepthvsDisease, file="Sequencing Depth vs. Disease.pdf", width=11.69, height=8.27, useDingbats=FALSE)


# CRC ---------------------------------------------------------------
CRC <- subset_samples(StomaMerged, Diagnosis_v4 == "C-CRC")
setwd("~/Dropbox/Ongoing Analysis/StomaMerged_November_2018/January2020/")
mkdir("CRC") # Make Directory to save generated figure in it
setwd("~/Dropbox/Ongoing Analysis/StomaMerged_November_2018/January2020/CRC")

# Beta Diversity in only CRC Samples
CRC_Ord_PCoA_bray = ordinate(CRC, "NMDS", distance="bray")
Plot_CRC <- plot_ordination(CRC, CRC_Ord_PCoA_bray, color = "Anatomic_Location_2") + geom_point(size=5, alpha=1)
Plot_CRC_Modified <- Plot_CRC + scale_colour_manual(values=c("#CC6666", "#00CC00","#3399FF", "#615350")) + scale_fill_manual(values=c("#CC6666", "#00CC00","#3399FF", "#615350")) + geom_point(size=6, alpha=1)  + stat_ellipse(aes(group=Anatomic_Location_2)) + theme_bw()
Plot_CRC_Modified 
ggsave(Plot_CRC_Modified, file="Beta_diversity_CRC_Bray_AnatomicLocation.pdf", width=11.69, height=8.27, useDingbats=FALSE)


#  Beta Diversity CRC Bag samples
CRC <- subset_samples(StomaMerged, Diagnosis_v2 == "CRC" & Sample_Style == "Bag")
CRC_Ord_PCoA_bray = ordinate(CRC, "PCoA", distance="bray")
Plot_CRC <- plot_ordination(CRC, CRC_Ord_PCoA_bray, color = "Anatomic_Location_2") + geom_point(size=5, alpha=1)
Plot_CRC_Modified <- Plot_CRC + scale_colour_manual(values=c("#CC6666", "#00CC00","#3399FF", "#615350")) + scale_fill_manual(values=c("#CC6666", "#00CC00","#3399FF", "#615350")) + geom_point(size=6, alpha=1)  + stat_ellipse(aes(group=Anatomic_Location)) + theme_bw()
Plot_CRC_Modified 
ggsave(Plot_CRC_Modified, file="Beta_diversity_CRC_Bray_AnatomicLocation_Bag.pdf", width=11.69, height=8.27, useDingbats=FALSE)

# Statistical analysis of this dataset
ent = prune_taxa(!(taxa_names(CRC) %in% "-1"), CRC)
ent = prune_taxa(taxa_sums(ent) > 0.0, ent)
ent <- transform_sample_counts(ent, function(x) x/sum(x))
ent <- subset_samples(ent, !is.na(CRC))
df = data.frame(sample_data(ent))
BC = distance(ent, "bray")
Stat_Beta_Diversity <- adonis(BC ~ Anatomic_Location*Age, data=df)
Stat_Beta_Diversity

#  Beta Diversity CRC Stoma samples
CRC <- subset_samples(StomaMerged, Diagnosis_v2 == "CRC" & Sample_Style == "Stoma")
CRC_Ord_PCoA_bray = ordinate(CRC, "PCoA", distance="bray")
Plot_CRC <- plot_ordination(CRC, CRC_Ord_PCoA_bray, color = "Anatomic_Location") + geom_point(size=5, alpha=1)
Plot_CRC_Modified <- Plot_CRC + scale_colour_manual(values=c("#CC6666", "#00CC00","#3399FF", "#615350")) + scale_fill_manual(values=c("#CC6666", "#00CC00","#3399FF", "#615350")) + geom_point(size=6, alpha=1)  + stat_ellipse(aes(group=Anatomic_Location)) + theme_bw()
Plot_CRC_Modified 
ggsave(Plot_CRC_Modified, file="Beta_diversity_CRC_Stoma.pdf", width=11.69, height=8.27, useDingbats=FALSE)

# Statistical analysis of this dataset
ent = prune_taxa(!(taxa_names(CRC) %in% "-1"), CRC)
ent = prune_taxa(taxa_sums(ent) > 0.0, ent)
ent <- transform_sample_counts(ent, function(x) x/sum(x))
ent <- subset_samples(ent, !is.na(CRC))
df = data.frame(sample_data(ent))
BC = distance(ent, "bray")
Stat_Beta_Diversity <- adonis(BC ~ Anatomic_Location*Age, data=df)
Stat_Beta_Diversity

# Alpha Diversity in CRC Samples
CRC <- subset_samples(StomaMerged, Diagnosis_v2 == "CRC")
my_comparisons <- list( c("Ileostoma", "Colostoma_Left"),
                        c("Ileostoma", "Colostoma_Right"),
                        c("Colostoma_Right", "Colostoma_Left"))
alpha_diversity_StomaMerged_plot <- plot_richness(CRC, color = "Anatomic_Location_2", x = "Anatomic_Location_2", measures = c("Shannon","Simpson")) + geom_boxplot() + stat_compare_means(comparisons = my_comparisons) + theme_bw() 
alpha_diversity_StomaMerged_plot
ggsave(alpha_diversity_StomaMerged_plot, file="Alpha_diversity_CRC_Diagnosis_plot.pdf", width=11.27, height=11.69, useDingbats=FALSE)

#  Alpha Diversity CRC Bag samples
CRC <- subset_samples(StomaMerged, Diagnosis_v2 == "CRC" & Sample_Style == "Bag")
my_comparisons <- list( c("Ileostoma", "Colostoma_Left"),
                        c("Ileostoma", "Colostoma_Right"),
                        c("Colostoma_Right", "Colostoma_Left"))
alpha_diversity_StomaMerged_plot <- plot_richness(CRC, color = "Anatomic_Location_2", x = "Anatomic_Location_2", measures = c("Shannon","Simpson")) + geom_boxplot() + stat_compare_means(comparisons = my_comparisons) + theme_bw() 
alpha_diversity_StomaMerged_plot
ggsave(alpha_diversity_StomaMerged_plot, file="Alpha_diversity_CRC_AnatomicLocation_Bag_samples.pdf", width=11.27, height=11.69,  useDingbats=FALSE)

#  Alpha Diversity CRC Stoma samples
CRC <- subset_samples(StomaMerged, Diagnosis_v2 == "CRC" & Sample_Style == "Stoma")
my_comparisons <- list( c("Ileostoma", "Colostoma_Left"),
                        c("Ileostoma", "Colostoma_Right"),
                        c("Colostoma_Right", "Colostoma_Left"))
alpha_diversity_StomaMerged_plot <- plot_richness(CRC, color = "Anatomic_Location_2", x = "Anatomic_Location_2", measures = c("Shannon","Simpson")) + geom_boxplot() + stat_compare_means(comparisons = my_comparisons) + theme_bw() 
alpha_diversity_StomaMerged_plot
ggsave(alpha_diversity_StomaMerged_plot, file="Alpha_diversity_CRC_AnatomicLocation_Stoma_samples.pdf", width=11.27, height=11.69,  useDingbats=FALSE)


# Taxonomy Analysis Using MaAsLin2 ----------------------------------------
# CD
input_data_bacteria = read.table("StomaSamples_Chip44_68_MappingFile_5000_L6.txt", header=TRUE,row.names="sample")
input_metadata = read.table("CD_metadata.txt", header=TRUE,row.names="sample")

input_metadata_subset <- input_metadata  %>%  
  tibble::rownames_to_column('sample') %>% 
  dplyr::filter(Diagnosis_v3 %in% c("CD"))  %>% 
  tibble::column_to_rownames('sample')
input_data<-select(input_metadata_subset, Age,Gender,BMI_Categorized,Col_Ile,Anatomic_Location,Sample_Style,Smoking_Status,Antibiotics_3months)

fit_data <- Maaslin2(
  input_data_bacteria, input_metadata, 'CD_relativeAbundance', transform = "LOG",normalization="CSS", max_significance=0.2, min_abundance=0.00000001, min_prevalence=0.1,
  fixed_effects = c('Gender','Col_Ile', 'Sample_Style'),
  random_effects = c('Smoking_Status', 'Antibiotics_3months', 'BMI_Categorized','Age'),
  standardize = FALSE)

input_data_bacteria = read.table("bacteria_count.txt", header=TRUE,row.names="sample")
input_metadata = read.table("CD_metadata.txt", header=TRUE,row.names="sample")
fit_data <- Maaslin2(
  input_data_bacteria, input_metadata, 'CD_BacterialCount', transform = "LOG",normalization="CSS", max_significance=0.2, min_abundance=0.00000001, min_prevalence=0.1,
  fixed_effects = c('Gender','Col_Ile', 'Sample_Style'),
  random_effects = c('Smoking_Status', 'Antibiotics_3months', 'BMI_Categorized','Age'),
  standardize = FALSE)

# CRC
input_data_bacteria = read.table("bacteria.txt", header=TRUE,row.names="sample")
input_metadata = read.table("data_CRC.txt", header=TRUE,row.names="sample")

fit_data <- Maaslin2(
  input_data_bacteria, input_metadata, 'CRC_relativeAbundance', transform = "LOG",normalization="CLR", max_significance=0.2, min_abundance=0.00000001, min_prevalence=0.1,
  fixed_effects = c('Gender','Col_Ile', 'Sample_Style'),
  random_effects = c('Smoking_Status', 'Antibiotics_3months', 'BMI_Categorized','Age'),
  standardize = FALSE)


input_data_bacteria = read.table("bacteria_count.txt", header=TRUE,row.names="sample")
input_metadata = read.table("data_CRC.txt", header=TRUE,row.names="sample")

fit_data <- Maaslin2(
  input_data_bacteria, input_metadata, 'CRC_BacterialCount', transform = "LOG",normalization="CLR", max_significance=0.2, min_abundance=0.00000001, min_prevalence=0.1,
  fixed_effects = c('Gender','Col_Ile', 'Sample_Style'),
  random_effects = c('Smoking_Status', 'Antibiotics_3months', 'BMI_Categorized','Age'),
  standardize = FALSE)



# Beta diversity of Bacterial Count samples -------------------------------
# https://fromthebottomoftheheap.net/2012/04/11/customising-vegans-ordination-plots/
# https://fukamilab.github.io/BIO202/06-B-constrained-ordination.html



data.les.pcoa=read.csv("Stoma_BacterialCount_withmetadata_nochronic_nocontrol.csv")
row.names(data.les.pcoa) <- paste(data.les.pcoa[,1],data.les.pcoa[,2],sep="_")
head(data.les.pcoa)

data.les.pcoa <- data.les.pcoa  %>%  
  dplyr::filter(Diagnosis_v2!=("Other_C"))


data2 <- data.les.pcoa[,16:85]
data1 <- data.les.pcoa[,2:15]

mod <- cca(data2 ~ Diagnosis_v2 +  BMI +  Gender + Col_Ile + Sample_Style + Smoking_Status + Antibiotics_3months, data1)
mod <- cca(data2 ~  BMI +  Gender + Col_Ile + Sample_Style + Smoking_Status + Antibiotics_3months, data1)
plot(mod, scaling = 3)

with(data1, levels(Col_Ile))
scl <- 3 ## scaling = 3
colvec <- c("black", "red", "green", "blue")
colvec <- c("green", "blue")

with(data1, points(mod, display = "sites", pch=20, col = colvec[Diagnosis_v2],
                   scaling = scl, bg = colvec[Diagnosis_v2]))

with(data1, points(mod, display = "sites", pch=20, col = colvec[Col_Ile],
                   scaling = scl, bg = colvec[Col_Ile]))

head(with(data1, colvec[Diagnosis_v2]))
head(with(data1, colvec[Col_Ile]))
# text(mod, display = "species", scaling = scl, cex = 0.8, col = "black")
with(data1, legend("bottomright", legend = levels(Diagnosis_v2), bty = "n",
                   col = colvec, pch = 21, pt.bg = colvec))

with(data1, legend("bottomright", legend = levels(Col_Ile), bty = "n",
                   col = colvec, pch = 21, pt.bg = colvec))

with(data1, ordihull(mod, Diagnosis_v2, kind="sd", conf=0.95, lwd=2, col=1:4, draw = "polygon",
                     label=TRUE, border=TRUE, scaling = 3))
with(data1, ordibar(mod, Diagnosis_v2, kind="se", conf=0.95, lwd=2, col=1:4,
                    label=TRUE))

with(data1, ordihull(mod, Col_Ile, kind="sd", conf=0.95, lwd=2, col=1:4, draw = "polygon",
                     label=TRUE, border=TRUE, scaling = 3))
with(data1, ordibar(mod, Col_Ile, kind="se", conf=0.95, lwd=2, col=1:4,
                    label=TRUE))


spe.hel <- decostand(data2, "hellinger")
bc<-vegdist(spe.hel, method="bray", binary=FALSE) 

set.seed(100)
bci.mds<-metaMDS(spe.hel, distance = "bray", k = 2)

# extract x and y coordinates from MDS plot into new dataframe, so you can plot with ggplot 
MDS_xy <- data.frame(bci.mds$points)
bci.mds$stress # 0.09184608

# colour by island
MDS_Plot<-ggplot(MDS_xy, aes(MDS1, MDS2, col=data1$Diagnosis_v2)) + theme_bw() + geom_point(size=5, alpha=1,aes(shape=data1$Anatomic_Location)) + theme_bw() + 
  stat_ellipse(aes(x = MDS1,y=MDS2,lty=data1$Diagnosis_v2,fill=data1$Diagnosis_v2), geom="polygon",level=0.8,alpha=0.2)
MDS_Plot
ggsave(MDS_Plot, file="MDS_Plot for Relative Abundance.pdf", width=11.69, height=8.27,useDingbats=FALSE)

ggplot(MDS_xy, aes(MDS1, MDS2, col=data1$Col_Ile)) + geom_point(size=2) + theme_bw() + geom_point(size=5, alpha=1) + theme_bw() + 
  stat_ellipse(aes(x = MDS1,y=MDS2,lty=data1$Col_Ile,fill=data1$Col_Ile), geom="polygon",level=0.8,alpha=0.2)



# Surgery Samples (Bandit Cohort) -----------------------------------------

setwd("~/Dropbox/Ongoing Analysis/StomaMerged_November_2018/Bandit Combined")
# Phyloseq Analysis
# Loading libraries
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


# Loading necessary files for phyloseq #
# Set the size of X and Y axis
My_Theme = theme(
  axis.title.x = element_text(size = 18),
  axis.text.x = element_text(size = 18),
  axis.title.y = element_text(size = 18),
  axis.text.y = element_text(size = 18))


# Load the qiime 2 files into phyloseq object
Bandit <-qza_to_phyloseq(features="table_Chip1_7.qza", tree="rooted-tree_Chip1_7.qza", taxonomy="taxonomy_Chip1_7.qza", metadata="Chip1_7_Bandit_v2.txt")
map <- sample_data(Bandit) #build sample data object
write.csv(map,"map.csv")   #write table of sample data object
saveRDS(Bandit, "Bandit.rds") # Save rds
Bandit = readRDS("Bandit.rds") # Read rds


# Subseting for the samples used in this part of the study
Bandit<-subset_samples(Bandit, Patient %in% c("P11","P13","P14","P15","P16","P17","P18","P19", "P20", "P24","P26","P27","P29","P32"))

# Overview ----------------------------------------------------------------

sdt = data.table(as(sample_data(Bandit), "data.frame"), TotalReads = sample_sums(Bandit), keep.rownames = TRUE)
write.csv(sdt, "MappingFilewithTotalReads.csv")


Bandit_subset<-subset_samples(Bandit, !BeforeAfter %in% c("1","2","1y"))

Bandit_filtered = prune_samples(sample_sums(Bandit_subset)>=930, Bandit_subset)
Bandit_filtered_rarefy <- rarefy_even_depth(Bandit_filtered)

# Beta Diversity
Bandit_Ord_PCoA_bray = ordinate(Bandit_subset, "PCoA", distance="bray") #Calculate Bray-Curtis dissimilarity
saveRDS(Bandit_Ord_PCoA_bray, "Bandit_Ord_PCoA_bray.rds")
Bandit_Ord_PCoA_unweighted = ordinate(Bandit_subset, "PCoA", distance="bray",weighted=FALSE) #Calculate Bray-Curtis dissimilarity
Bandit_Ord_PCoA_weighted = ordinate(Bandit_subset, "PCoA", distance="unifrac",weighted=TRUE) #Calculate Bray-Curtis dissimilarity

# PcoA Plot
Plot_Bandit <- plot_ordination(Bandit_subset, Bandit_Ord_PCoA_bray, color="BeforeAfter") + geom_point(size=5, alpha=1)  + theme_bw() +  scale_shape_manual(values=rep(c(0:3,5:6), times=5)) + theme_bw()  +  
  stat_ellipse(aes(group=BeforeAfter,fill=BeforeAfter), geom="polygon",alpha=0.3)
Plot_Bandit
ggsave(Plot_Bandit, file="Plot_Bandi Before and After Surgery.pdf", width=11.69, height=8.27,useDingbats=FALSE)

Plot_Bandit <- plot_ordination(Bandit_subset, Bandit_Ord_PCoA_bray, color="TimePoint") + geom_point(size=5, alpha=1)  + theme_bw() +  scale_shape_manual(values=rep(c(0:3,5:6), times=5))  + stat_ellipse(aes(group=TimePoint))
Plot_Bandit
ggsave(Plot_Bandit, file="Plot_Bandit for Time Points.pdf", width=11.69, height=8.27,useDingbats=FALSE)

# Clean up data before adonis
ent = prune_taxa(!(taxa_names(Bandit_subset) %in% "-1"), Bandit_subset) 
ent = prune_taxa(taxa_sums(ent) > 0.0, ent)
ent <- transform_sample_counts(ent, function(x) x/sum(x))
ent <- subset_samples(ent, !is.na(Bandit_subset))
df = data.frame(sample_data(ent)) 
BC = distance(ent, "bray")

# First method
all <- adonis(BC ~ Patient*BeforeAfter + Patient*Location + Location*BeforeAfter, data=df)
all
beta <- betadisper(BC, df$Before_After)
permutest(beta)

# Extraction of Alpha diversity indices for downstream analysis -----------
alpha_diversity_Bandit <-estimate_richness(Bandit, split = TRUE, measures = NULL)
data_cbind <- cbind(sample_data(Bandit), alpha_diversity_Bandit)
saveRDS(data_cbind, "alpha_diversity_Bandit.rds")
# alpha_diversity_Bandit= readRDS("alpha_diversity_Bandit.rds")



# Some Alpha Diversity Plots ----------------------------------------------

alpha_diversity_Bandit_plot <- plot_richness(Bandit, color = "TimePoint", x = "Patient", measures = c("Shannon","Simpson")) + geom_boxplot(aes(fill = TimePoint), alpha = 0.2) + theme_bw() 
alpha_diversity_Bandit_plot
ggsave(alpha_diversity_Bandit_plot, file="alpha_diversity_Bandit_plot.pdf", width=14.69, height=8.27,useDingbats=FALSE)

my_comparisons <- list( c("Before", "After") )
sample_data(Bandit) <- within(sample_data(Bandit), 
                              Before_After <- factor(BeforeAfter, 
                                                     levels=names(sort(table(BeforeAfter), 
                                                                       decreasing=FALSE))))
alpha_diversity_Bandit_plot <- plot_richness(Bandit, color = "BeforeAfter", x = "BeforeAfter", measures = c("Shannon","Simpson")) + geom_boxplot(aes(fill = BeforeAfter), alpha = 0.2) + stat_compare_means(comparisons = my_comparisons) # + scale_colour_manual(values=c("#CC6666", "#3399FF")) + scale_fill_manual(values=c("#CC6666", "#3399FF")) + theme_bw() 
alpha_diversity_Bandit_plot

# This includes all samples. We just wanted to see the whole picture and statistically we compared later the longitudinal samples.
data_cbind_Bandit <- cbind(sample_data(Bandit_subset), alpha_diversity_Bandit)
anova_Shannon <- aov(Shannon ~ BeforeAfter, data_cbind_Bandit)
summary(anova_Shannon)
attach(data_cbind_Bandit)
pairwise.wilcox.test(data_cbind_Bandit$Shannon, data_cbind_Bandit$BeforeAfter, p.adjust.method="fdr")


# Stoma Over Time ---------------------------------------------------------

Stoma <-subset_samples(Bandit, Location %in% c("Stoma") & !TimePoint %in% c("E0"))
# Beta Diversity
Stoma_Ord_PCoA_bray = ordinate(Stoma, "PCoA", distance="bray") #Calculate Bray-Curtis dissimilarity

# Clean up data before adonis
ent = prune_taxa(!(taxa_names(Stoma) %in% "-1"), Stoma) 
ent = prune_taxa(taxa_sums(ent) > 0.0, ent)
ent <- transform_sample_counts(ent, function(x) x/sum(x))
ent <- subset_samples(ent, !is.na(Stoma))
df = data.frame(sample_data(ent)) 
BC = distance(ent, "bray")

# First method
all <- adonis(BC ~ Patient*TimePoint, data=df)
all
beta <- betadisper(BC, df$TimePoint)
permutest(beta)

# PcoA Plot
Plot_Stoma <- plot_ordination(Stoma, Stoma_Ord_PCoA_bray, color="TimePoint") + geom_point(size=5, alpha=1)  + theme_bw() + scale_colour_manual(values=c("lightblue","forestgreen", "red", "purple", "orange", "beige", "gold4", "black")) + 
  scale_fill_manual(values=c("lightblue","forestgreen", "red", "purple", "orange", "beige", "gold4", "black")) +
  scale_shape_manual(values=rep(c(15:25), times=7)) + stat_ellipse(aes(group=TimePoint), geom="polygon",alpha=0.05)
Plot_Stoma
ggsave(Plot_Stoma, file="Plot_Stoma Over Time.pdf", width=11.69, height=8.27,useDingbats=FALSE)

Plot_Stoma <- plot_ordination(Stoma, Stoma_Ord_PCoA_bray, color="Patient", shape="Patient") + geom_point(size=5, alpha=1)  + theme_bw() +
  scale_shape_manual(values=rep(c(15:25), times=7)) + stat_ellipse(aes(group=Patient), geom="polygon",alpha=0.05)
Plot_Stoma
ggsave(Plot_Stoma, file="Plot_Stoma Over per Patient.pdf", width=11.69, height=8.27,useDingbats=FALSE)

alpha_diversity_Stoma_plot <- plot_richness(Stoma, color = "TimePoint", x = "TimePoint", measures = c("Shannon")) + geom_boxplot(aes(fill = TimePoint), alpha = 0.2) + theme_bw()  + scale_colour_manual(values=c("lightblue","forestgreen", "red", "purple", "orange", "beige", "gold4", "black")) + 
  scale_fill_manual(values=c("lightblue","forestgreen", "red", "purple", "orange", "beige", "gold4", "black"))
alpha_diversity_Stoma_plot
ggsave(alpha_diversity_Stoma_plot, file="alpha_diversity_Stoma_plot.pdf", width=14.69, height=8.27,useDingbats=FALSE)


# Rectum Content Over Time ------------------------------------------------

RectumContent <-subset_samples(Bandit, Location %in% c("3c"))
# Beta Diversity
RectumContent_Ord_PCoA_bray = ordinate(RectumContent, "PCoA", distance="bray") #Calculate Bray-Curtis dissimilarity


# Clean up data before adonis
ent = prune_taxa(!(taxa_names(RectumContent) %in% "-1"), RectumContent) 
ent = prune_taxa(taxa_sums(ent) > 0.0, ent)
ent <- transform_sample_counts(ent, function(x) x/sum(x))
ent <- subset_samples(ent, !is.na(RectumContent))
df = data.frame(sample_data(ent)) 
BC = distance(ent, "bray")

# First method
all <- adonis(BC ~ Patient*TimePoint, data=df)
all
beta <- betadisper(BC, df$TimePoint)
permutest(beta)

# Post Hoc comparison
pairwise.adonis(BC, sample_data(RectumContent)$TimePoint)

adonis.site <- phyloseq_to_adonis(
  physeq =RectumContent, 
  dist = "bray", 
  formula = "TimePoint"
)


TP <-subset_samples(RectumContent, TimePoint %in% c("E0","E1"))
# Clean up data before adonis
ent = prune_taxa(!(taxa_names(TP) %in% "-1"), TP) 
ent = prune_taxa(taxa_sums(ent) > 0.0, ent)
ent <- transform_sample_counts(ent, function(x) x/sum(x))
ent <- subset_samples(ent, !is.na(TP))
df = data.frame(sample_data(ent)) 
BC = distance(ent, "bray")

# First method
all <- adonis(BC ~ Patient*TimePoint, data=df)
all
beta <- betadisper(BC, df$TimePoint)
permutest(beta)


# PcoA Plot
Plot_RectumContent <- plot_ordination(RectumContent, RectumContent_Ord_PCoA_bray, color="TimePoint", shape="TimePoint") + geom_point(size=5, alpha=1)  + theme_bw() + scale_colour_manual(values=c("lightblue","forestgreen", "red", "purple", "orange", "pink", "gold4", "black")) + 
  scale_fill_manual(values=c("lightblue","forestgreen", "red", "purple", "orange", "pink", "gold4", "black")) +
  scale_shape_manual(values=rep(c(15:25), times=7)) + stat_ellipse(aes(group=TimePoint), geom="polygon",alpha=0.05)
Plot_RectumContent
ggsave(Plot_RectumContent, file="Plot_RectumContent Over Time.pdf", width=11.69, height=8.27,useDingbats=FALSE)

Plot_RectumContent <- plot_ordination(RectumContent, RectumContent_Ord_PCoA_bray, color="Patient", shape="Patient") + geom_point(size=5, alpha=1)  + theme_bw() +
  scale_shape_manual(values=rep(c(15:25), times=7)) + stat_ellipse(aes(group=Patient))
Plot_RectumContent
ggsave(Plot_RectumContent, file="Plot_RectumContent Over  per Patient.pdf", width=11.69, height=8.27,useDingbats=FALSE)


alpha_diversity_RectumContent_plot <- plot_richness(RectumContent, color = "TimePoint", x = "TimePoint", measures = c("Shannon")) + geom_boxplot(aes(fill = TimePoint), alpha = 0.6) + theme_bw()  + scale_colour_manual(values=c("lightblue","forestgreen", "red", "purple", "orange", "pink", "gold4", "black")) + 
  scale_fill_manual(values=c("lightblue","forestgreen", "red", "purple", "orange", "pink", "gold4", "black"))
alpha_diversity_RectumContent_plot
ggsave(alpha_diversity_RectumContent_plot, file="alpha_diversity_RectumContent_plot.pdf", width=14.69, height=8.27,useDingbats=FALSE)


# Alpha Diversity Analysis with Linear Mixed Model for longitudinal samples
library(tidyverse)
library(lme4)
library(lmerTest)  # Note, we might have to install this first! 

pn_long <- read_csv("AlphaDiversityRectumContent.csv")
mem_rt_4 <- lmer(Diversity ~ TimePoint + BMI + Antibiotic + Gender + (1|PatientID) + (1|SmokingStatus), data = pn_long)
coef(summary(mem_rt_4))
head(coef(mem_rt_4)$PatientID)
summary(mem_rt_4)

pn_long <- read_csv("AlphaDiversityIleum.csv")
mem_rt_4 <- lmer(Diversity ~ TimePoint + BMI + Antibiotic + Gender + (1|PatientID) + (1|SmokingStatus), data = pn_long)
coef(summary(mem_rt_4))
head(coef(mem_rt_4)$PatientID)
summary(mem_rt_4)
                               
                               
### CRC vs CD comparison
# Calling metadata and bacteria file
input_data_bacteria = read.table("bacteria.txt", header=TRUE,row.names="sample")
input_metadata = read.table("CD_CRC_metadata.txt", header=TRUE,row.names="sample")

# Calling subseting the metadata for Ileostoma analysis
input_metadata_subset <- input_metadata  %>%  
  tibble::rownames_to_column('sample') %>% 
  dplyr::filter(Col_Ile %in% c("Ileostoma"))  %>% 
  tibble::column_to_rownames('sample')

input_data<-select(input_metadata_subset,Diagnosis_v3,Age,Gender,BMI_Categorized,Smoking_Status,Antibiotics_3months)

fit_data <- Maaslin2(
  input_data_bacteria, input_data, 'CD_CRC_Ileostoma',transform = "LOG",normalization="CLR", min_abundance=0.0000001, min_prevalence=0.1,max_significance=0.2,
  fixed_effects = c('Diagnosis_v3', 'Gender', 'Sample_Style','BMI','Age'),
  random_effects = c('Smoking_Status', 'Antibiotics_3months'),
  standardize = FALSE)



# Calling subseting the metadata for Colostoma analysis
input_metadata_subset <- input_metadata  %>%  
  tibble::rownames_to_column('sample') %>% 
  dplyr::filter(Col_Ile %in% c("Colostoma"))  %>% 
  tibble::column_to_rownames('sample')

input_data<-select(input_metadata_subset,Diagnosis_v3,Age,Gender,BMI_Categorized,Smoking_Status,Antibiotics_3months)

fit_data <- Maaslin2(
  input_data_bacteria, input_data, 'CD_CRC_Colostoma',transform = "LOG",normalization="CLR", min_abundance=0.000001, min_prevalence=0.1,max_significance=0.2,
  fixed_effects = c('Diagnosis_v3', 'Gender', 'Sample_Style','BMI','Age'),
  random_effects = c('Smoking_Status', 'Antibiotics_3months'),
  standardize = FALSE)

