# This file is used for model construction. 

# set-up 
if (!dir.exists("results")) {
  dir.create("results", recursive = TRUE)
}

wgcna_path <- "data/wgcna_results.rds"
if (!file.exists(processed_path)) {
  system("Rscript 01_wgcna.R")
}

#data_dir <- "data" 
#if (!dir.exists(data_dir)) { 
#  dir.create(data_dir, recursive=TRUE)
#}

# load data
data <- readRDS(wgcna_path)
clinical_factors <- data$pheno
eigengenes <- data$eigengenes

# create dataframe
model_data <- data.frame(pathologic_response = clinical_factors$pathologic_response,eigengenes)

# adjust models to explore what works for training set
m1 <- glm(pathologic_response ~ MEblack + MEyellow + MEpurple + MEbrown, data=model_data, family=binomial)
summary(m1)
