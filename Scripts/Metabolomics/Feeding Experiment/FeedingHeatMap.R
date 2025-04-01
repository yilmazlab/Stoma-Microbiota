annotation_table= read.table("/Users/bahti/Dropbox/Ongoing\ Analysis/StomaMerged_November_2018/Metabolomics/bandit\ subset\ 2021/heatmap_increase_decrrease/annotation.txt")
annotation_row= read.table("/Users/bahti/Dropbox/Ongoing\ Analysis/StomaMerged_November_2018/Metabolomics/bandit\ subset\ 2021/heatmap_increase_decrrease/annotation_row.txt")
Taxa<-read.table("/Users/bahti/Dropbox/Ongoing\ Analysis/StomaMerged_November_2018/Metabolomics/bandit\ subset\ 2021/heatmap_increase_decrrease/feedingmetabolite.increasing.decreasing.txt", row.names=NULL)
row.names(Taxa) <- make.unique(as.character(Taxa$row.names))
Taxa$row.names <- NULL

ann_colors = list(
  TimePoint = c(Zero ="#fff5f0", Two="#fcbba1", Four="#fc9272", Six="#ef3b2c", Eight="#cb181d", Ten = "#a50f15", Fortyeight="#67000d"),
  Patient = c(P1="#B5BBE3", P2="#aaf2ff", P3="#3399FF",P4="#efb94c", P6="#cc6666", P7="#f4fcfc"),
  Overtime =c(Increase="red", Decrease ="green"),
  Source = c(Bacteria="brown", Drug ="royalblue",Plant ="mediumspringgreen", Protein="pink", General="yellow", Endogenous="tan", Food="firebrick2", Unknown="white"))


rwbcols <- c( "#9c0e09","#bdaf1e", "#04bd2f")


feedingmetabolites <- pheatmap(Taxa,
                                 fontsize_row = 7, fontsize_col = 10, 
                                 cluster_cols = TRUE, cluster_rows = TRUE,
                                 color = turbo(15),
                                 cellwidth = 25, cellheight = 8,# cutree_rows = 4,  
                                 annotation_col=annotation_table, annotation_row=annotation_row, annotation_colors = ann_colors, border_color =  "grey60")
feedingmetabolites
ggsave(feedingmetabolites, file="/Users/bahti/Dropbox/Ongoing\ Analysis/StomaMerged_November_2018/Metabolomics/bandit\ subset\ 2021/heatmap_increase_decrrease/feedingmetabolite.increasing.decreasing.V2.pdf", width=20, height=30, useDingbats=FALSE,limitsize = FALSE)


