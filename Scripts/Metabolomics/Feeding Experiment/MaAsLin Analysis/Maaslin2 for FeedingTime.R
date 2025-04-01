setwd("~/Desktop") # Add the data to Desktop
library(Maaslin2) # Call the pipeline


input_metadata = read.table("metadata_t0_t2.txt", header=TRUE,row.names="sample")
input_metadata = read.table("metadata_t0_t4.txt", header=TRUE,row.names="sample")
input_data = read.table("metabolites.txt", header=TRUE,row.names="sample")


# Overall analysis
fit_data <- Maaslin2(
  input_data, input_metadata, 'Metabolomic_t0_t4', analysis_method="LM", transform = "LOG", normalization ="NONE", max_significance=0.2, min_abundance=0.001, min_prevalence=0.3,
  fixed_effects = c('Age', 'TimePoint'),
  random_effects = c('Patient'),
  standardize = FALSE)
