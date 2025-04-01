# Pavian
options(shiny.maxRequestSize=2000*2048^2)
pavian::runApp(port=5000)

setwd("~/Dropbox/Ongoing Analysis/StomaMetagenomicMergedMerged_November_2018/Variant Analysis")

setwd("~/Dropbox/Ongoing Analysis/StomaMerged_November_2018/Variant Analysis/kreports/Braken")

# Loading libraries -------------------------------------------------------
library(xlsx)
library(plyr)
library(dplyr)
library(gdata)
library(ggplot2)
library(extrafont)
library(scales)
library(MASS) 
library(PMCMR)
library(psych)
require(stats)
library(RColorBrewer)
library(ggrepel)
library(phyloseq)
library(vegan)
library(varhandle)
library(ggplot2)
library(ape)
library(phangorn)
library(picante)
library(ggthemes)
library("dada2")
library(gridExtra)
library(reshape2)
library(reshape)
library(caret)
library(phyloseqGraphTest)
library(magrittr)
library(grid)
library(randomForest)
library(knitr)
library(data.table)
library(microbiomeSeq) 
source("http://peterhaschke.com/Code/multiplot.R")
loadfonts(device = "pdf", quiet = FALSE) #font_import()
## Color codes for disease =  #CC6666", "#00CC00","#3399FF"
require(gdata)
library(tidyverse)
library(RColorBrewer)
library(taxonomizr)
library(vegan)
library(ggtree)
library(treeio)
library(ggridges)
library(kableExtra)
library(ggrepel)
library(patchwork)
library(ggpubr)

My_Theme = theme(
  axis.title.x = element_text(size = 18),
  axis.text.x = element_text(size = 18),
  axis.title.y = element_text(size = 18),
  axis.text.y = element_text(size = 18))


# Demographic -------------------------------------------------------------

StomaMetagenomicPercentage = read.xls ("StomaMetagenomicSamples.xlsx", sheet = 1, header = TRUE)
StomaMetagenomicPercentage=read.csv("StomaMetagenomicSamples_Percentage.csv")
StomaMetagenomicPercentage$Viral_reads<- format(StomaMetagenomicPercentage$Viral_reads, digit=10)

NumberOfReads <- ggplot(data=StomaMetagenomicPercentage, aes(x=PatientID, y=Number_of_raw_reads, fill=PatientID)) + # geom_boxplot(notch = FALSE) +  
  ggtitle("Number_of_raw_reads per sample") +
  scale_y_continuous(name = "Reads per sample", trans = 'log10') +
  scale_x_discrete(name = "PatientID") + geom_jitter(colour = "black", size = 8, width = 0.2, height = 0.1) + 
  labs(fill = "PatientID") + 
  theme(axis.title=element_text(face="plain", size="30", color="black"), 
        axis.text.x = element_text(colour="grey20",size=20,angle=0,hjust=.5,vjust=.5,face="plain"), 
        axis.text.y = element_text(colour="grey20",size=40,angle=0,hjust=1,vjust=0,face="plain"), 
        legend.text=element_text(face="plain", size="30", color="black"),
        legend.title=element_text(face="plain", size="30", color="black")) + theme_bw() + My_Theme
NumberOfReads
ggsave(NumberOfReads, file="Number of Reads Per Sample.pdf", width=40, height=15, useDingbats=FALSE)

Classified_reads <- ggplot(data=StomaMetagenomicPercentage, aes(x=PatientID, y=Classified_reads, fill=PatientID)) + #geom_boxplot()  +
  ggtitle("% Classified reads per sample") +
  scale_y_continuous(name = "Classified reads per sample", trans = 'log10') +
  scale_x_discrete(name = "PatientID") + geom_point(size=8, alpha=1, aes(colour = factor(Classified_reads))) +
  labs(fill = "PatientID") + 
  theme(axis.title=element_text(face="plain", size="30", color="black"), 
        axis.text.x = element_text(colour="grey20",size=20,angle=0,hjust=.5,vjust=.5,face="plain"), 
        axis.text.y = element_text(colour="grey20",size=40,angle=0,hjust=1,vjust=0,face="plain"), 
        legend.text=element_text(face="plain", size="30", color="black"),
        legend.title=element_text(face="plain", size="30", color="black")) + theme_bw() 
Classified_reads
ggsave(Classified_reads, file="Classified reads percentage of Stoma Patients.pdf", width=25, height=15, useDingbats=FALSE)

Microbial_reads <- ggplot(data=StomaMetagenomicPercentage, aes(x=PatientID, y=Microbial_reads, fill=PatientID)) + # geom_boxplot()  +
  ggtitle("Number of Microbial Reads per Sample") +
  scale_y_continuous(name = "Microbial reads per sample") +
  scale_x_discrete(name = "PatientID") + geom_point(size=8, alpha=1,aes(colour = factor(Microbial_reads))) +
  labs(fill = "PatientID") + 
  theme(axis.title=element_text(face="plain", size="30", color="black"), 
        axis.text.x = element_text(colour="grey20",size=40,angle=45,hjust=.5,vjust=.5,face="plain"), 
        axis.text.y = element_text(colour="grey20",size=40,angle=0,hjust=1,vjust=0,face="plain"), 
        legend.text=element_text(face="plain", size="30", color="black"),
        legend.title=element_text(face="plain", size="30", color="black")) + theme_bw() 
Microbial_reads
ggsave(Microbial_reads, file="Microbial reads percentage of Stoma Patients.pdf", width=25, height=15, useDingbats=FALSE)

Bacterial_reads <- ggplot(data=StomaMetagenomicPercentage, aes(x=PatientID, y=Bacterial_reads, fill=PatientID)) + # geom_boxplot()  +
  ggtitle("Number of Bacterial Reads per Sample") +
  scale_y_continuous(name = "Bacterial reads per sample") +
  scale_x_discrete(name = "PatientID") + geom_point(size=8, alpha=1,aes(colour = factor(Bacterial_reads))) +
  labs(fill = "PatientID") + 
  theme(axis.title=element_text(face="plain", size="30", color="black"), 
        axis.text.x = element_text(colour="grey20",size=20,angle=0,hjust=.5,vjust=.5,face="plain"), 
        axis.text.y = element_text(colour="grey20",size=40,angle=0,hjust=1,vjust=0,face="plain"), 
        legend.text=element_text(face="plain", size="30", color="black"),
        legend.title=element_text(face="plain", size="30", color="black")) + theme_bw() 
Bacterial_reads
ggsave(Bacterial_reads, file="Bacterial reads percentage of Stoma Patients.pdf", width=25, height=15, useDingbats=FALSE)


Reads_m <- melt(Reads, id.vars = "PatientID", measure.vars = c("Chordate_reads", "Unclassified_reads", "Bacterial_reads", "Viral_reads", "Fungal_reads", "Protozoan_reads"))
Reads_m <- melt(Reads2, id.vars = "PatientID")
theme_set(theme_classic())



Viral_reads <- ggplot(data=StomaMetagenomicPercentage, aes(x=PatientID, y=Viral_reads, fill=PatientID)) + # geom_boxplot()  +
  ggtitle("Number of Viral reads per sample") +
  scale_y_continuous(name = "Viral reads per sample") +
  scale_x_discrete(name = "PatientID") + geom_point(size=8, alpha=1,aes(colour = factor(Viral_reads))) +
  labs(fill = "PatientID") + 
  theme(axis.title=element_text(face="plain", size="30", color="black"), 
        axis.text.x = element_text(colour="grey20",size=20,angle=0,hjust=.5,vjust=.5,face="plain"), 
        axis.text.y = element_text(colour="grey20",size=40,angle=0,hjust=1,vjust=0,face="plain"), 
        legend.text=element_text(face="plain", size="30", color="black"),
        legend.title=element_text(face="plain", size="30", color="black")) + theme_bw() 
Viral_reads
ggsave(Viral_reads, file="Viral reads percentage of Stoma Patients.pdf", width=25, height=15, useDingbats=FALSE)


Protozoan_reads <- ggplot(data=StomaMetagenomicPercentage, aes(x=PatientID, y=Protozoan_reads, fill=PatientID)) + # geom_boxplot()  +
  ggtitle( "Number of Protozoan reads per sample") +
  scale_y_continuous(name = "Protozoan reads per sample") +
  scale_x_discrete(name = "PatientID") + geom_point(size=8, alpha=1,aes(colour = factor(Protozoan_reads))) +
  labs(fill = "PatientID") + 
  theme(axis.title=element_text(face="plain", size="30", color="black"), 
        axis.text.x = element_text(colour="grey20",size=20,angle=0,hjust=.5,vjust=.5,face="plain"), 
        axis.text.y = element_text(colour="grey20",size=40,angle=0,hjust=1,vjust=0,face="plain"), 
        legend.text=element_text(face="plain", size="30", color="black"),
        legend.title=element_text(face="plain", size="30", color="black")) + theme_bw() 
Protozoan_reads
ggsave(Protozoan_reads, file="Viral reads percentage of Stoma Patients.pdf", width=25, height=15, useDingbats=FALSE)



# Call the files to combine metadata samples ------------------------------

read_kreport <- function(report_path, dataset_name) {
  df <- read_tsv(report_path,
                 col_names = c("percentage_rooted", "n_reads_rooted", "n_reads", "rank", "taxid", "sci_name")) %>%
    add_column(dataset = dataset_name)
  return(df)
}

read_kreports <- function(named_paths) {
  df <- tibble()
  for (path in named_paths) {
    name <- names(which(named_paths == path))
    t <- read_kreport(path, name)
    df <- bind_rows(df, t)
  }
  return(df)
}

metadata <- read.table("metadata.tsv", sep = "\t", row.names = 1, header = TRUE)
samples <- c("SF1",	"SF2",	"SF3",	"SF4",	"SF5",	"SF6",	"SF7",	"SF8",	"SF9",	"SF10",	"SF11",	"SF12",	"SF13",	"SF14",	"SF15",	"SF16",	"SF17",	"SF18",	"SF19",	"SF20",	"SF21",	"SF22",	"SF23",	"SF24",	"SF25",	"SF26",	"SF27",	"SF34",	"SF35",	"SF36",	"SF37")
bracken_paths <- dir("~/Dropbox/Ongoing Analysis/StomaMerged_November_2018/Variant Analysis/kreports/Braken/", pattern = "bracken_species.kreport.tsv", full.names = TRUE) %>% set_names(samples)

# Joining the metadata and samples
bracken_data <- read_kreports(bracken_paths) %>%
  inner_join(metadata, by = "dataset")
saveRDS(bracken_data, file = "results/bracken_data.Rds") # save the bracken_data for later to call
bracken_data = readRDS("results/bracken_data.rds") # read the bracken_data


# Basic Plots of Statistic of The Data ------------------------------
## The ranks are a one letter code: (U)nclassified, (D)omain, (K)ingdom, (P)hylum, (C)lass, (O)rder, (F)amily, (G)enus, or (S)pecies. All other ranks are '-'

### Just to extract the number of Domain and Unclassified in the list
bracken_data %>%
  filter(rank == "D") %>%
  select(n_reads_rooted, sci_name, dataset) %>%
  group_by(sci_name) %>%
  summarise(total = sum(n_reads_rooted)) %>%
  kable() %>%
  kable_styling(bootstrap_options = "striped", full_width = F)

domain<-bracken_data %>% filter(rank == "D")
write.csv(domain, "Domain Read Count.csv")
# Plot for Domain
bracken_data %>% filter(rank == "D") %>%
  ggplot(.) +
  geom_boxplot(aes(sci_name, n_reads_rooted), alpha = 1,
               color = "black") + 
  geom_jitter(aes(sci_name, n_reads_rooted, color = Patient), size=4, stroke=1,
              width = 0.3) +
  scale_color_brewer(palette = 'Paired', name = 'Patients') +
  ggtitle("Distribution of Domains across datasets") + theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Domain") +
  ylab("Counts") + My_Theme

# Plot for Phylum
bracken_data %>% filter(rank == "P") %>%
  ggplot(.) +
  geom_boxplot(aes(sci_name, n_reads_rooted), alpha = 1,
               color = "black") +
  geom_jitter(aes(sci_name, n_reads_rooted, color = Patient), size=4, stroke=1,
              width = 0.3) +
  scale_color_brewer(palette = 'Paired', name = 'Patients') +
  ggtitle("Distribution of Phylum across datasets") +
  xlab("Phylum") +
  ylab("Counts")



# Analyzing Bacterial Fraction ------------------------------
## There are clearly quite a lot of Eukaryotic sequences. We'll leave them aside for now and focus on the bacteria and viruses.

is_taxonomy <- function(query_tax, target_tax, database) {
  # check if query tax has target tax as a parent
  query_id <- getId(query_tax, database)
  query_lineage <- getTaxonomy(query_id, database)
  if (target_tax %in% query_lineage) {
    paste(query_lineage)
    return(TRUE)
  }
  else {
    return(FALSE)
  }
}
is_taxonomy_vec <- Vectorize(is_taxonomy, vectorize.args = "query_tax")


### To generate taxonomy.sqlite file
# getNamesAndNodes(outDir = "db",
#     url = "ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz",
#    fileNames = c("names.dmp", "nodes.dmp"))
read.nodes.sql("~/Dropbox/Ongoing Analysis/Wild-Type Mice vs Wild Mice Analysis/Metagenomics/Metagenomics/Bracken/data/taxdump/nodes.dmp", overwrite=TRUE,
               sqlFile = "~/Dropbox/Ongoing Analysis/Wild-Type Mice vs Wild Mice Analysis/Metagenomics/Metagenomics/Bracken/data/taxonomy.sqlite")
read.names.sql("~/Dropbox/Ongoing Analysis/Wild-Type Mice vs Wild Mice Analysis/Metagenomics/Metagenomics/Bracken/data/taxdump/names.dmp", overwrite=TRUE,
               sqlFile = "~/Dropbox/Ongoing Analysis/Wild-Type Mice vs Wild Mice Analysis/Metagenomics/Metagenomics/Bracken/data/taxonomy.sqlite")



bracken_phyla <- bracken_data %>%
  filter(rank == "P")
is_bact <- is_taxonomy_vec(bracken_phyla$sci_name, "Bacteria", "~/Dropbox/Ongoing Analysis/Wild-Type Mice vs Wild Mice Analysis/Metagenomics/Metagenomics/Bracken/data/taxonomy.sqlite")
saveRDS(is_bact, file = "results/is_bact.Rds") # save the is_bact for later to call
is_bact = readRDS("results/is_bact.rds") # read the is_bact


bracken_phyla <- bracken_phyla[is_bact, ]

most_abundant_bracken_phyla <- bracken_phyla %>%
  group_by(sci_name) %>%
  mutate(count = sum(n_reads_rooted)) %>%
  ungroup() %>%
  distinct(sci_name, .keep_all = TRUE) %>%
  select(sci_name, count) %>%
  top_n(10)

top_bracken_phyla <- bracken_phyla %>%
  filter(sci_name %in% most_abundant_bracken_phyla$sci_name)

background <- bracken_phyla %>%
  anti_join(top_bracken_phyla) %>%
  group_by(dataset) %>%
  mutate(count = sum(n_reads_rooted)) %>%
  select(-n_reads_rooted, -sci_name, -percentage_rooted) %>%
  dplyr::rename(n_reads_rooted = count) %>%
  add_column(sci_name = "Other") %>%
  distinct(dataset, .keep_all = TRUE) %>%
  ungroup()

top_phyla <- top_bracken_phyla %>%
  bind_rows(background)

top_phyla %>%
  arrange(desc(n_reads_rooted)) %>%
  kable() %>%
  # kableExtra::group_rows(index = table(top_phyla$dataset)) %>%
  kable_styling(bootstrap_options = "striped", full_width = F)

# cab2d6 (lightpurple) = Bacteroidetes or (6a3d9a)
# steelblue = Proteobacteria
# ff7f00 (orange) = Firmicutes
# mediumseagreen = Cyanobacteria
palette <- c("#6a3d9a", "mediumseagreen", "#ff7f00", "red", "steelblue")
# names(palette) <- levels(tpf)

phyla <- ggplot(top_phyla) +
  geom_bar(aes(dataset, n_reads_rooted, fill = sci_name),
           stat = "identity", position = "fill", alpha = 1) +
  # theme_minimal() +
  theme(text = element_text(family = "serif", size = 7)) +
  # scale_fill_manual(values = palette, name = "Phylum") +
  facet_grid(~ Patient, scales = "free") +
  xlab("Samples") +
  ylab("Percentage of bacterial reads") +
  theme(axis.line = element_line(colour = "black"), axis.title=element_text(face="bold", size="50", color="black"),
         axis.text.x = element_text(colour="grey20",size=1,angle=45,hjust=.5,vjust=.5,face="plain"),
         axis.text.y = element_text(colour="grey20",size=25,angle=0,hjust=1,vjust=0,face="plain"),
         legend.text=element_text(face="bold", size="50", color="black"),
         legend.title=element_text(face="bold", size="50", color="black")) + 
  ggtitle("Main Phyla Changes in Patients") + theme_bw() +  theme(plot.title = element_text(hjust=0.5))
phyla
ggsave(phyla, file="Main Phyla Changes in Patients.pdf", width=15.69, height=8.27, useDingbats=FALSE)


# Chanced top_n to 20
most_abundant_bracken_phyla <- bracken_phyla %>%
  group_by(sci_name) %>%
  mutate(count = sum(n_reads_rooted)) %>%
  ungroup() %>%
  distinct(sci_name, .keep_all = TRUE) %>%
  select(sci_name, count) %>%
  top_n(10)

most_abundant_bracken_phyla

top_bracken_phyla <- bracken_phyla %>%
  filter(sci_name %in% most_abundant_bracken_phyla$sci_name)

background <- bracken_phyla %>%
  anti_join(top_bracken_phyla) %>%
  group_by(dataset) %>%
  mutate(count = sum(n_reads_rooted)) %>%
  select(-n_reads_rooted, -sci_name, -percentage_rooted) %>%
  dplyr::rename(n_reads_rooted = count) %>%
  add_column(sci_name = "Other") %>%
  distinct(dataset, .keep_all = TRUE) %>%
  ungroup()

top_phyla <- top_bracken_phyla %>%
  bind_rows(background)

top_phyla %>%
  arrange(desc(n_reads_rooted)) %>%
  kable() %>%
  # kableExtra::group_rows(index = table(top_phyla$dataset)) %>%
  kable_styling(bootstrap_options = "striped", full_width = F)


phyla <- ggplot(top_phyla) +
  geom_bar(aes(dataset, n_reads_rooted, fill = sci_name),
           stat = "identity", position = "fill", alpha = 1) +
  # theme_minimal() +
  theme(text = element_text(family = "serif", size = 7)) +
  # scale_fill_manual(values = palette, name = "Phylum") +
  facet_grid(~ Patient, scales = "free") +
  xlab("Samples") +
  ylab("Percentage of bacterial reads") +
  theme(axis.line = element_line(colour = "black"), axis.title=element_text(face="bold", size="50", color="black"),
        axis.text.x = element_text(colour="grey20",size=1,angle=45,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=25,angle=0,hjust=1,vjust=0,face="plain"),
        legend.text=element_text(face="bold", size="50", color="black"),
        legend.title=element_text(face="bold", size="50", color="black")) + 
  ggtitle("Main Phyla Changes in Patients") + theme_bw() +  theme(plot.title = element_text(hjust=0.5))
phyla
ggsave(phyla, file="Main Phyla Changes in Patients.pdf", width=15.69, height=8.27, useDingbats=FALSE)


# Proteobacteria fraction -------------------------------------------------

## First we create subsets oif the data containing the bacterial genera and species read counts, while retaining phylum information.
## These datasets will be useful for investigating subsets of the bacterial population (aka zooming in)

get_phylum_genus <-  function(taxid) {
  phylum <- getTaxonomy(taxid, "~/Dropbox/Ongoing Analysis/StomaMerged_November_2018/Variant Analysis/kreports/Braken/taxonomy.sqlite") %>%
    as_tibble() %>%
    select(phylum, genus) %>%
    add_column(taxid)
  return(phylum)
}

phylum <- bracken_data %>%
  filter(rank == "P")

genera <- bracken_data %>%
  filter(rank == "G")

families <- bracken_data %>%
  filter(rank == "F")

species <- bracken_data %>%
  filter(rank == "S")

genera_to_phyla <- get_phylum_genus(genera$taxid) %>%
  distinct(taxid, .keep_all = TRUE)
families_to_phyla <- get_phylum_genus(families$taxid) %>%
  distinct(taxid, .keep_all = TRUE)
species_to_phyla <- get_phylum_genus(species$taxid) %>%
  distinct(taxid, .keep_all = TRUE)
phylum_to_phyla <- get_phylum_genus(phylum$taxid) %>%
  distinct(taxid, .keep_all = TRUE)

bact_genera <- genera %>%
  left_join(genera_to_phyla, by = "taxid") %>%
  filter(phylum %in% bracken_phyla$sci_name)
bact_families <- families %>%
  left_join(families_to_phyla, by = "taxid") %>%
  filter(phylum %in% bracken_phyla$sci_name)
bact_species <- species %>%
  left_join(species_to_phyla, by = "taxid") %>%
  filter(phylum %in% bracken_phyla$sci_name)
bact_phylum <- phylum %>%
  left_join(phylum_to_phyla, by = "taxid") %>%
  filter(phylum %in% bracken_phyla$sci_name)


saveRDS(bracken_data, file = "~/Desktop/Bracken/results/bracken.Rds")
saveRDS(bact_genera, file = "~/Desktop/Bracken/results/genera.Rds")
saveRDS(bact_families, file = "~/Desktop/Bracken/results/families.Rds")
saveRDS(bact_species, file = "~/Desktop/Bracken/results/species.Rds")
saveRDS(bact_phylum, file = "~/Desktop/Bracken/results/phylum.Rds")

# Now we can visualise the genera and species distribution of the datasets
most_abundant_species <- bact_species %>%
  group_by(sci_name) %>%
  mutate(count = sum(percentage_rooted)) %>%
  ungroup() %>%
  distinct(sci_name, .keep_all = TRUE) %>%
  select(sci_name, count) %>%
  top_n(20)

top_s <- bact_species %>%
  filter(sci_name %in% most_abundant_species$sci_name) %>%
  group_by(taxid) %>%
  mutate(total = sum(percentage_rooted)) %>%
  ungroup()

# define a palette to map species name to phyla
tp <- top_s %>% arrange(total) %>% select(phylum) %>% filter(row_number() %% 6 == 0)
tpf <- as_factor(tp$phylum)
n_colors <- length(levels(tpf)) # How many colors you need
colors <- scales::brewer_pal('qual')
palette <- c("mediumpurple", "darkgreen", "red", "orange")
names(palette) <- levels(tpf)

ggplot(top_s) +
  geom_bar(aes(reorder(sci_name, total, max), percentage_rooted, fill = TimePoint),
           stat = "identity", position = "stack", alpha = 1) +
  theme_minimal() +
  theme(text = element_text(family = "serif")) +
  facet_grid(~ Patient, scales = "free") +
  coord_flip() +
  ylim(0, 1) +
  ylab("Bacterial fraction (in %)") +
  xlab("Species name") +
  theme(axis.text.y = element_text(colour = palette[tpf], size = 7))
ggsave("~/Desktop/Bracken/plots/top_50_sp.pdf", device = "pdf", dpi = 600, width = 20, height = 30, units = "cm")

# Now we can visualise the phylum and phylum distribution of the datasets
most_abundant_phylum <- bact_phylum %>%
  group_by(sci_name) %>%
  mutate(count = sum(percentage_rooted)) %>%
  ungroup() %>%
  distinct(sci_name, .keep_all = TRUE) %>%
  select(sci_name, count) %>%
  top_n(10)

top_s <- bact_phylum %>%
  filter(sci_name %in% most_abundant_phylum$sci_name) %>%
  group_by(taxid) %>%
  mutate(total = sum(percentage_rooted)) %>%
  ungroup()

# define a palette to map phylum name to phyla
tp <- top_s %>% arrange(total) %>% select(phylum) %>% filter(row_number() %% 6 == 0)
tpf <- as_factor(tp$phylum)
n_colors <- length(levels(tpf)) # How many colors you need
colors <- scales::brewer_pal('qual')
palette <- c("mediumpurple", "darkgreen", "red", "orange")
names(palette) <- levels(tpf)

ggplot(top_s) +
  geom_bar(aes(reorder(sci_name, total, max), percentage_rooted, fill = WildMiceSPF),
           stat = "identity", position = "stack", alpha = 1) +
  theme_minimal() +
  theme(text = element_text(family = "serif")) +
  facet_grid(~ WildMiceSPF, scales = "free") +
  coord_flip() +
  ylim(0, 1) +
  ylab("Bacterial fraction (in %)") +
  xlab("Phyla name") +
  theme(axis.text.y = element_text(colour = palette[tpf], size = 7))
ggsave("~/Desktop/Bracken/plots/top_50_sp.pdf", device = "pdf", dpi = 600, width = 20, height = 30, units = "cm")



# Now we can visualise the genera and genera distribution of the datasets
most_abundant_genera <- bact_genera %>%
  group_by(sci_name) %>%
  mutate(count = sum(percentage_rooted)) %>%
  ungroup() %>%
  distinct(sci_name, .keep_all = TRUE) %>%
  select(sci_name, count) %>%
  top_n(50)

top_s <- bact_genera %>%
  filter(sci_name %in% most_abundant_genera$sci_name) %>%
  group_by(taxid) %>%
  mutate(total = sum(percentage_rooted)) %>%
  ungroup()

# define a palette to map genera name to phyla
tp <- top_s %>% arrange(total) %>% select(phylum) %>% filter(row_number() %% 6 == 0)
tpf <- as_factor(tp$phylum)
n_colors <- length(levels(tpf)) # How many colors you need
colors <- scales::brewer_pal('qual')
palette <- c("mediumpurple","orange", "steelblue", "darkgreen", "red" )
names(palette) <- levels(tpf)

ggplot(top_s) +
  geom_bar(aes(reorder(sci_name, total, max), percentage_rooted, fill = Patient),
           stat = "identity", position = "stack", alpha = 1) +
  theme_minimal() +
  theme(text = element_text(family = "serif")) +
  facet_grid(~ TimePoint, scales = "free") +
  coord_flip() +
  ylim(0, 1) +
  ylab("Bacterial fraction (in %)") +
  xlab("Genera name") +
  theme(axis.text.y = element_text(colour = palette[tpf], size = 7))
ggsave("~/Desktop/Bracken/plots/top_50_genera_per_group.pdf", device = "pdf", dpi = 600, width = 50, height = 30, units = "cm")

write.csv(top_s, "~/Desktop/Bracken/results/top_s.csv")







