setwd("~/Desktop")
library(phyloseq)
source("~/Desktop/metaphlanToPhyloseq.R")
indat = read.delim("plant_metagenome.txt",sep = "")
inmetadat = read.delim("feeding_metadata.txt",sep = "")

My_Theme = theme(
  axis.title.x = element_text(size = 18),
  axis.text.x = element_text(size = 18),
  axis.title.y = element_text(size = 18),
  axis.text.y = element_text(size = 18))

Plant= metaphlanToPhyloseq(indat, inmetadat, simplenames = TRUE, 
                              roundtointeger = FALSE)
rank_names(Plant)

par(mar = c(30, 4, 0, 0) + 0.1) # make more room on bottom margin
barplot(sort(taxa_sums(Plant), TRUE)[1:30]/nsamples(Plant), las=2)

Plant.int = Plant
otu_table(Plant.int) = round(otu_table(Plant)*1e5)
alpha_meas = c("Shannon", "Simpson")
p <- plot_richness(Plant.int, "Patient",  measures=alpha_meas, color="Patient") + geom_boxplot(alpha = 0) + theme_bw() 
p 
ggsave(p, file="alpha_diversity_for all patient.pdf", width=16.69, height=8.27,useDingbats=FALSE)

alpha_diversity_Plant <-estimate_richness(Plant.int, split = TRUE, measures = NULL)
data_cbind <- cbind(sample_data(Plant.int), alpha_diversity_Plant)
write.csv(data_cbind, "alphadiversity_Plant.csv")

ord= ordinate(Plant, method="PCoA", distance="gower")
plot<-plot_ordination(Plant, ord, color="Patient", shape="TimePoint") + geom_point(size=5, alpha=1) + theme_bw() + scale_shape_manual(values=rep(c(15:19), times=6)) + My_Theme + geom_text_repel(size=5,  aes(label = PatientID), max.overlaps = Inf) # + stat_ellipse(aes(group=Patient))
plot  
ggsave(plot, file="Beta Diversity of Plant for Metagenomic Data.pdf", width=14, height=9,useDingbats=FALSE)


rank_names(Plant)
subset_taxa(Plant, !is.na(Strain))
(Plant.sp_strain = subset_taxa(Plant, !is.na(Species)))
taxonomy.level = apply(tax_table(Plant), 1, function(x) sum(!is.na(x)))
Plant.phy = prune_taxa(taxonomy.level==2, Plant)
taxa_names(Plant.phy)

Plant.species = prune_taxa(taxonomy.level==6, Plant)
taxa_names(Plant.species)

plot_heatmap(Plant.species, method="PCoA", distance="bray", sample.label="PatientID", sample.order = "PatientID", low="#000033", high="#CCFF66", trans = log_trans(2))
heatmap(otu_table(Plant.species), Colv = NA)

plot(hclust(phyloseq::distance(Plant, method="bray")), 
     main="Bray-Curtis Dissimilarity", xlab="", sub = "")

ent = prune_taxa(!(taxa_names(Plant) %in% "-1"), Plant)
ent = prune_taxa(taxa_sums(ent) > 0.0, ent)
ent <- transform_sample_counts(ent, function(x) x/sum(x))
ent <- subset_samples(ent, !is.na(Plant))
df = data.frame(sample_data(ent)) 
BC = distance(ent, "bray")
Stat_Beta_Diversity <- adonis(BC ~ Patient*TimePoint, data=df)
Stat_Beta_Diversity
beta <- betadisper(BC, df$Disease)
permutest(beta)

# Maaslin -----------------------------------------------------------------
library(Maaslin2)
# Analyzing at family level
input_data_family = indat
input_metadata = inmetadat
input_metadata_subset<- input_metadata %>%
  dplyr::filter(TimePoint %in% c("0h","2h","4h","6h", "8h", "10h"))

input_data<-select(input_metadata_subset,TimePoint_num)

fit_data <- Maaslin2(
  input_data_family, input_data, 'PlantTaxonomy', max_significance=1, min_abundance=0.001, min_prevalence=0.5,
  standardize = FALSE)






