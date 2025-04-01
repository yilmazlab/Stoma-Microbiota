# https://fromthebottomoftheheap.net/2012/04/11/customising-vegans-ordination-plots/
# https://fukamilab.github.io/BIO202/06-B-constrained-ordination.html

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
library(xlsx)
library(openxlsx)


setwd("~/Dropbox/Ongoing Analysis/StomaMerged_November_2018/Metabolomics/RawData")
data.les.pcoa <- read.xls("quantile_normalized_data_with_chronic.xls", sheet = 1, header = TRUE)
data.les.pcoa <- read.xls("quantile_normalized_zscored_data_with_chronic.xls", sheet = 1, header = TRUE) #C-CRC
data.les.pcoa <- read.xls("quantile_normalized_zscored_data_with_chronic.xls", sheet = 2, header = TRUE) #C-CRC
data.les.pcoa <- read.xls("quantile_normalized_zscored_data_with_chronic.xls", sheet = 3, header = TRUE) #COMPLEX
data.les.pcoa <- read.xls("quantile_normalized_zscored_data_with_chronic.xls", sheet = 4, header = TRUE) #IBD
data.les.pcoa <- read.csv(file = 'BanditFeedingTime_Zscore.csv')
row.names(data.les.pcoa) <- paste(data.les.pcoa[,1],data.les.pcoa[,2],sep="_")
head(data.les.pcoa)

data2 <- data.les.pcoa[,15:256]
data1 <- data.les.pcoa[,2:14]

pheatmap(data2, color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),cluster_rows = FALSE, cluster_cols = TRUE, cellwidth = 10, cellheight = 2)


mod <- cca(data2 ~ PatientID, data1)
plot(mod, scaling = 2)

with(data1, levels(AnatomicLocation))
scl <- 3 ## scaling = 3
colvec <- c("red", "green", "blue")


with(data1, points(mod, display = "sites", pch=20, col = colvec[Diagnosis_v2],
                   scaling = scl, bg = colvec[Diagnosis_v2]))

with(data1, points(mod, display = "sites", pch=20, col = colvec[AnatomicLocation],
                   scaling = scl, bg = colvec[AnatomicLocation]))

head(with(data1, colvec[Diagnosis_v2]))
head(with(data1, colvec[AnatomicLocation]))
# text(mod, display = "species", scaling = scl, cex = 0.8, col = "black")
with(data1, legend("bottomright", legend = levels(Diagnosis_v2), bty = "n",
                   col = colvec, pch = 21, pt.bg = colvec))

with(data1, legend("bottomright", legend = levels(AnatomicLocation), bty = "n",
                   col = colvec, pch = 21, pt.bg = colvec))

with(data1, ordihull(mod, Diagnosis_v2, kind="sd", conf=0.95, lwd=2, col=1:4, draw = "polygon",
                     label=TRUE, border=TRUE, scaling = 3))
with(data1, ordibar(mod, Diagnosis_v2, kind="se", conf=0.95, lwd=2, col=1:4,
                    label=TRUE))

with(data1, ordihull(mod, AnatomicLocation, kind="sd", conf=0.95, lwd=2, col=1:4, draw = "polygon",
                     label=TRUE, border=TRUE, scaling = 3))
with(data1, ordibar(mod, AnatomicLocation, kind="se", conf=0.95, lwd=2, col=1:4,
                    label=TRUE))


spe.hel <- decostand(data2, "hellinger")
bc<-vegdist(spe.hel, method="bray", binary=FALSE) 

PC<-prcomp(data2, center = TRUE, scale. = FALSE)
summary(PC)
plot(PC$x[,1], PC$x[,2], cex=2, col=factor(data1$PatientID), xlab="PC1", ylab="PC2", pch=16, main="feedingTime", las=1) 
Stat_Beta_Diversity <- adonis(PC_x ~ PatientID+TimePoint, data=data1)
Stat_Beta_Diversity

PC_x <- data.frame(PC$x)

MDS_Plot<-ggplot(PC_x, aes(PC1, PC2)) + geom_point(size=6, alpha=1,aes(color=data1$Diagnosis_v2)) + theme_bw() + scale_fill_brewer(type="qual", palette="Set1") + scale_colour_brewer(type="qual", palette="Set1")  + stat_ellipse(aes(color=data1$Diagnosis_v2, fill=data1$Diagnosis_v2), geom="polygon",level=0.9,alpha=0.1) # + geom_polygon(aes(fill=AnatomicLocation)) 
MDS_Plot 
ggsave(MDS_Plot, file="All MDS_Plot for Metabolomic Bandit Dataset Feeding Time.pdf", width=11.69, height=8.27,useDingbats=FALSE) 
  
  

set.seed(100)
bci.mds<-metaMDS(PC, distance = "bray", k = 3)

# extract x and y coordinates from MDS plot into new dataframe, so you can plot with ggplot 
MDS_xy <- data.frame(bci.mds$points)
bci.mds$stress # 0.09184608

# colour by island
MDS_Plot<-ggplot(MDS_xy, aes(MDS1, MDS2, col=data1$Diagnosis_v2)) + theme_bw() + geom_point(size=5, alpha=1,aes(shape=data1$AnatomicLocation)) + theme_bw() + 
  stat_ellipse(aes(x = MDS1,y=MDS2,lty=data1$Diagnosis_v2,fill=data1$Diagnosis_v2), geom="polygon",level=0.95,alpha=0.1)
MDS_Plot
ggsave(MDS_Plot, file="MDS_Plot for Relative Abundance.pdf", width=11.69, height=8.27,useDingbats=FALSE)

ggplot(MDS_xy, aes(MDS1, MDS2, col=data1$AnatomicLocation)) + geom_point(size=2) + theme_bw() + geom_point(size=5, alpha=1) + theme_bw() + 
  stat_ellipse(aes(x = MDS1,y=MDS2,lty=data1$AnatomicLocation,fill=data1$AnatomicLocation), geom="polygon",level=0.8,alpha=0.2)

simpleRDA <- cca(data2 ~  Diagnosis_v2 + BMI_Cat +  AnatomicLocation  + Gender, data=data1)
summary(simpleRDA)
screeplot(simpleRDA) #bstick not available for constrained ordinations

# canonical coefficients
coef(simpleRDA)

# unadjusted R^2 retreived from the rda result
R2 <- RsquareAdj(simpleRDA)$r.squared
R2

# adjusted R^2
R2adj <- RsquareAdj(simpleRDA)$adj.r.squared
R2adj 

# Triplot: three different entities in the plot: sites, response variables and explanatory variables (arrowheads are on the explanatory variables)
# Scaling 1
plot(simpleRDA, scaling=3, main="Triplot RDA matrix ~ env - scaling 1 - wa scores")

# arrows for species are missing, so lets add them without heads so they look different than the explanatory variables
spe.sc <- scores(simpleRDA, choices=0:5, scaling=1, display="sp")
arrows(0,0,spe.sc[,1], spe.sc[,2], length=0, lty=1, col='red')


library(devtools)
install_github('fawda123/ggord')
library(ggord)
ggord(simpleRDA, data1$AnatomicLocation, ellipse = FALSE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + stat_ellipse(aes(fill=data1$AnatomicLocation), geom="polygon",level=0.8,alpha=0.2)

vif.cca(simpleRDA)

# variance inflation factors in the RDA
anova.cca(simpleRDA, step=1000)
Stat_Beta_Diversity <- adonis(data2 ~ Diagnosis_v2 + BMI_Cat +  AnatomicLocation  + Gender, data=data1)
Stat_Beta_Diversity


data2.scaled <- scale(data2, # cols 13:16 are factors
                      center=TRUE, # center data by subtracting column means 
                      scale=TRUE)

str(data2.scaled)

env.pca <- rda(data2, scale=TRUE)    # rda function makes PCAs
# scale means divide by species unit variance

env.pca
ev <- env.pca$CA$eig
ev>mean(ev)
par(mfrow=c(2,1))
barplot(ev, main="eigenvalues", col="bisque", las=2)
abline(h=mean(ev), col="red")
legend("topright", "Average eignenvalue", lwd=1, col=2)
barplot(t(cbind(100*ev/sum(ev), bsm$p[n:1])), beside=TRUE, main="% variance", col=c("bisque", 2), las=2)

par(mfrow=c(1,2))
biplot(env.pca, scaling=1, main="PCA scaling 1") # scaling 1 for distances between objects
biplot(env.pca, main="PCA scaling 2") # correlation biplot, for seeing correlation between response variables (species) see angles. 


library(ade4)
pr.env2 <- prcomp(data2.scaled)
prin.env2 <- princomp(data2.scaled)
dudi.env2 <-dudi.pca(data2.scaled, center=F, scale=F, scannf=FALSE, nf=4) # ade4 library
pr.env2
prin.env2
dudi.env2

meaneig <-mean(dudi.env2$eig) # .98

library(factoextra)
eig0 <- fviz_eig(pr.env2, geom="bar", width=.4, title="prcomp")+theme(panel.grid=element_blank())
eig1 <- fviz_eig(prin.env2, geom="bar", width=.4, title="princomp")+theme(panel.grid=element_blank())
eig2 <- fviz_eig(dudi.env2, geom="bar", width=.4, title="dudi.pca")+theme(panel.grid=element_blank())

library(gridExtra)
grid.arrange(eig0,eig1, eig2, nrow=1)

p1 <- fviz_pca_ind(pr.env2, geom="point",
                   habillage= data1$Diagnosis_v2,
                   axes= 2:3)

p2 <- fviz_pca_ind(pr.env2, geom="point",
                   habillage= data1$Diagnosis_v2,
                   axes= 3:4)

p3 <- fviz_pca_ind(pr.env2, geom="point",
                   habillage= data1$Diagnosis_v2,
                   axes= c(1,3))

p4 <- fviz_pca_ind(pr.env2, geom="point",
                   habillage= data1$Diagnosis_v2,
                   axes= c(1,4))

p5 <- fviz_pca_ind(pr.env2, geom="point",
                   habillage= data1$Diagnosis_v2,
                   axes= c(2,4))

grid.arrange(p1,p2,p3,p4,p5)



ord <- metaMDS(data2)
ggord(ord, data1$Diagnosis_v2)

ord <- capscale(data2 ~ Diagnosis_v2, data1, dist = "bray")



# Subset Diseases ---------------------------------------------------------

# https://fromthebottomoftheheap.net/2012/04/11/customising-vegans-ordination-plots/
# https://fukamilab.github.io/BIO202/06-B-constrained-ordination.html


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
library(xlsx)
library(openxlsx)


setwd("~/Dropbox/Ongoing Analysis/StomaMerged_November_2018/Metabolomics/RawData")
data.les.pcoa <- read.xls("quantile_normalized_data_with_chronic.xls", sheet = 2, header = TRUE) #C-CRC
data.les.pcoa <- read.xls("quantile_normalized_data_with_chronic.xls", sheet = 3, header = TRUE) #COMPLEX
data.les.pcoa <- read.xls("quantile_normalized_data_with_chronic.xls", sheet = 4, header = TRUE) #IBD

row.names(data.les.pcoa) <- paste(data.les.pcoa[,1],data.les.pcoa[,2],sep="_")
head(data.les.pcoa)

data2 <- data.les.pcoa[,15:256]
data1 <- data.les.pcoa[,2:14]

mod <- cca(data2 ~ Ile_Col, data1)
plot(mod, scaling = 2)

with(data1, levels(Ile_Col))
scl <- 3 ## scaling = 3
colvec <- c("red", "green", "blue")


with(data1, points(mod, display = "sites", pch=20, col = colvec[Ile_Col],
                   scaling = scl, bg = colvec[Ile_Col]))

with(data1, points(mod, display = "sites", pch=20, col = colvec[Ile_Col],
                   scaling = scl, bg = colvec[Ile_Col]))

head(with(data1, colvec[Ile_Col]))
head(with(data1, colvec[Ile_Col]))
# text(mod, display = "species", scaling = scl, cex = 0.8, col = "black")
with(data1, legend("bottomright", legend = levels(Ile_Col), bty = "n",
                   col = colvec, pch = 21, pt.bg = colvec))

with(data1, legend("bottomright", legend = levels(Ile_Col), bty = "n",
                   col = colvec, pch = 21, pt.bg = colvec))

with(data1, ordihull(mod, Ile_Col, kind="sd", conf=0.95, lwd=2, col=1:4, draw = "polygon",
                     label=TRUE, border=TRUE, scaling = 3))
with(data1, ordibar(mod, Ile_Col, kind="se", conf=0.95, lwd=2, col=1:4,
                    label=TRUE))

with(data1, ordihull(mod, Ile_Col, kind="sd", conf=0.95, lwd=2, col=1:4, draw = "polygon",
                     label=TRUE, border=TRUE, scaling = 3))
with(data1, ordibar(mod, Ile_Col, kind="se", conf=0.95, lwd=2, col=1:4,
                    label=TRUE))


spe.hel <- decostand(data2, "hellinger")
bc<-vegdist(spe.hel, method="bray", binary=FALSE) 

set.seed(100)
bci.mds<-metaMDS(spe.hel, distance = "bray", k = 3)

# extract x and y coordinates from MDS plot into new dataframe, so you can plot with ggplot 
MDS_xy <- data.frame(bci.mds$points)
bci.mds$stress # 0.09184608

# colour by island
MDS_Plot<-ggplot(MDS_xy, aes(MDS1, MDS2, col=data1$Ile_Col)) + theme_bw() + geom_point(size=5, alpha=1,aes(shape=data1$Ile_Col)) + theme_bw() + 
  stat_ellipse(aes(x = MDS1,y=MDS2,lty=data1$Ile_Col,fill=data1$Ile_Col), geom="polygon",level=0.95,alpha=0.1) + scale_colour_manual(values=c("#8d8666", "#768cbd")) + scale_fill_manual(values=c("#8d8666", "#768cbd"))
MDS_Plot
ggsave(MDS_Plot, file="MDS_Plot for CRC.pdf", width=11.69, height=8.27,useDingbats=FALSE)

MDS_Plot<-ggplot(MDS_xy, aes(MDS1, MDS2, col=data1$Ile_Col)) + geom_point(size=2) + theme_bw() + geom_point(size=5, alpha=1) + theme_bw() + stat_ellipse(aes(x = MDS1,y=MDS2,lty=data1$Ile_Col,fill=data1$Ile_Col), geom="polygon",level=0.8,alpha=0.2) + scale_colour_manual(values=c("#8d8666", "#768cbd")) + scale_fill_manual(values=c("#8d8666", "#768cbd"))
MDS_Plot
ggsave(MDS_Plot, file="MDS_Plot for IBD.pdf", width=11.69, height=8.27,useDingbats=FALSE)

simpleRDA <- cca(data2 ~  Ile_Col + BMI_Cat + Gender + SampleType, data=data1)
summary(simpleRDA)
screeplot(simpleRDA) #bstick not available for constrained ordinations

# canonical coefficients
coef(simpleRDA)

# unadjusted R^2 retreived from the rda result
R2 <- RsquareAdj(simpleRDA)$r.squared
R2

# adjusted R^2
R2adj <- RsquareAdj(simpleRDA)$adj.r.squared
R2adj 

# Triplot: three different entities in the plot: sites, response variables and explanatory variables (arrowheads are on the explanatory variables)
# Scaling 1
plot(simpleRDA, scaling=2, main="Triplot RDA matrix ~ env - scaling 1 - wa scores")

# arrows for species are missing, so lets add them without heads so they look different than the explanatory variables
spe.sc <- scores(simpleRDA, choices=0:5, scaling=1, display="sp")
arrows(0,0,spe.sc[,1], spe.sc[,2], length=0, lty=1, col='red')


library(devtools)
install_github('fawda123/ggord')
library(ggord)
ggord(simpleRDA, data1$Ile_Col, ellipse = FALSE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + stat_ellipse(aes(fill=data1$Ile_Col), geom="polygon",level=0.8,alpha=0.2)

vif.cca(simpleRDA)

# variance inflation factors in the RDA
anova.cca(simpleRDA, step=1000)
Stat_Beta_Diversity <- adonis(data2 ~ Ile_Col + BMI_Cat +  Ile_Col  + Gender, data=data1)
Stat_Beta_Diversity


data2.scaled <- scale(data2, # cols 13:16 are factors
                      center=TRUE, # center data by subtracting column means 
                      scale=TRUE)

str(data2.scaled)

env.pca <- rda(data2, scale=TRUE)    # rda function makes PCAs
# scale means divide by species unit variance

env.pca
ev <- env.pca$CA$eig
ev>mean(ev)
par(mfrow=c(2,1))
barplot(ev, main="eigenvalues", col="bisque", las=2)
abline(h=mean(ev), col="red")
legend("topright", "Average eignenvalue", lwd=1, col=2)
barplot(t(cbind(100*ev/sum(ev), bsm$p[n:1])), beside=TRUE, main="% variance", col=c("bisque", 2), las=2)

par(mfrow=c(1,2))
biplot(env.pca, scaling=1, main="PCA scaling 1") # scaling 1 for distances between objects
biplot(env.pca, main="PCA scaling 2") # correlation biplot, for seeing correlation between response variables (species) see angles. 


library(ade4)
pr.env2 <- prcomp(data2.scaled)
prin.env2 <- princomp(data2.scaled)
dudi.env2 <-dudi.pca(data2.scaled, center=F, scale=F, scannf=FALSE, nf=4) # ade4 library
pr.env2
prin.env2
dudi.env2

meaneig <-mean(dudi.env2$eig) # .98

library(factoextra)
eig0 <- fviz_eig(pr.env2, geom="bar", width=.4, title="prcomp")+theme(panel.grid=element_blank())
eig1 <- fviz_eig(prin.env2, geom="bar", width=.4, title="princomp")+theme(panel.grid=element_blank())
eig2 <- fviz_eig(dudi.env2, geom="bar", width=.4, title="dudi.pca")+theme(panel.grid=element_blank())

library(gridExtra)
grid.arrange(eig0,eig1, eig2, nrow=1)

p1 <- fviz_pca_ind(pr.env2, geom="point",
                   habillage= data1$Ile_Col,
                   axes= 2:3)

p2 <- fviz_pca_ind(pr.env2, geom="point",
                   habillage= data1$Ile_Col,
                   axes= 3:4)

p3 <- fviz_pca_ind(pr.env2, geom="point",
                   habillage= data1$Ile_Col,
                   axes= c(1,3))

p4 <- fviz_pca_ind(pr.env2, geom="point",
                   habillage= data1$Ile_Col,
                   axes= c(1,4))

p5 <- fviz_pca_ind(pr.env2, geom="point",
                   habillage= data1$Ile_Col,
                   axes= c(2,4))

grid.arrange(p1,p2,p3,p4,p5)



ord <- metaMDS(data2)
ggord(ord, data1$Ile_Col)

ord <- capscale(data2 ~ Ile_Col, data1, dist = "bray")

-

