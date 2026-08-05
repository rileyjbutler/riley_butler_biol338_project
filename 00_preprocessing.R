# This file is used to retrieve and process the data for the WGCNA. 
# The following script is adopted from hadiazarabad's repository 'Cancer-Drug-Response' from file 'prepare_gse25066.R' 
# The input for WGCNA should have samples as rows and columns as genes with the cells showing the normalised expression values. 

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
gse <- getGEO("GSE25055", GSEMatrix=TRUE)
gse <- gse[[1]]

# get data tables 
expr <- exprs(gse)
feature <- fData(gse)
pheno <- pData(gse)

# expr contains the expression level of the probes for 22283 probes across the 310 training samples 
# feature contains 16 variables for the 22283 features/probes in the data set, including gene symbols and functions
# pheno contains 80 variables for the 310 samples, including patient information, biomarker status and chemotherapy response 


if (!"Gene Symbol" %in% colnames(feature)) {
  stop("Gene Symbol column not found in feature data")
}

# cleaning up gene symbols 
gene_symbols <- feature$'Gene Symbol'
gene_symbols <- ifelse(is.na(gene_symbols) | gene_symbols == "", feature$ID, gene_symbols)
gene_symbols <- sapply(strsplit(gene_symbols, "///"), `[`, 1)

# condense microarray data object to deal with duplicates and empty rows 
expr <- avereps(expr, ID=gene_symbols)
expr <- expr[rownames(expr) != "", ]

# choose and clean the response column containing pCR/RD status 
response_col <- "pathologic_response_pcr_rd:ch1"
if (!response_col %in% colnames(pheno)){
stop("The response column was not found in the phenotype data")  
}

response <- pheno[[response_col]]
response <- trimws(as.character(response))
keep_samples <- !is.na(response) & response %in% c("pCR", "RD")

# choose samples in the datatables that have response data (pCR or RD). This removes 4 samples. 
expr <- expr[, keep_samples] 
response <- response[keep_samples]
pheno <- pheno[keep_samples, , drop=FALSE]

# find gene variances and keep 3000 of the genes with highest variance
variances <- apply(expr, 1, var)
variances <- variances[!is.na(variances)]

n_top_genes <- min(3000, length(variances))
top_genes <- names(sort(variances, decreasing=TRUE))[seq_len(n_top_genes)]
expr_top <- expr[top_genes, ]

processed <- list(
  expression_raw = expr_top,
  response = response,
  phenotype = pheno,
  genes = top_genes
)

hist(expr, main = "Gene Expression Data", xlab = "Expression")
hist(expr_top, main = "Gene Expression Data", xlab = "Expression")

plotDensities(
  expr_top,
  legend = FALSE,
  main = "Expression density across samples"
)

saveRDS(processed, output_file)
