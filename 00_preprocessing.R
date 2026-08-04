# This file is used to retrieve and process the data for the WGCNA. 
# The following script is adopted from hadiazarabad's repository 'Cancer-Drug-Response' from file 'prepare_gse25066.R' 

# load libraries 
library(GEOquery)
library(limma)
library(dplyr)

# environment set-up - make data directory and output file path
data_dir <- "data" 
if (!dir.exists(data_dir)) { 
  dir.create(data_dir, recursive=TRUE)
}

output_file <- file.path(data_dir, "processed_data.rds")

# retrieve data 
gse <- getGEO("GSE25066", GSEMatrix=TRUE)
gse <- gse[[1]]

# get data tables 
expr <- exprs(gse)
feature <- fData(gse)
pheno <- pData(gse)

if (!"Gene Symbol" %in% colnames(feature)) {
  stop("Gene Symbol column not found in feature data")
}

# cleaning gene symbols 
gene_symbols <- feature$'Gene Symbol'
gene_symbols <- ifelse(is.na(gene_symbols) | gene_symbols == "", feature$ID, gene_symbols)
gene_symbols <- sapply(strsplit(gene_symbols, "///"), `[`, 1)

# condense microarray data object to deal with duplicates and empty rows 