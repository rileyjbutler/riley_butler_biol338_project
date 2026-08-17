# This file is used to retrieve and process the data for the WGCNA. 
# The following script is adopted from hadiazarabad's repository 'Cancer-Drug-Response' from file 'prepare_gse25066.R' 
# The input for WGCNA should have samples as rows and columns as genes with the cells showing the normalised expression values. 

# load libraries 
library(GEOquery)
library(limma)
library(dplyr)
library("WGCNA")

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

# find gene variances and keep 10000 of the genes with highest variance
variances <- apply(expr, 1, var)
variances <- variances[!is.na(variances)]

n_top_genes <- min(10000, length(variances))
top_genes <- names(sort(variances, decreasing=TRUE))[seq_len(n_top_genes)]
expr_top <- expr[top_genes, ]

# filtering for only useful variables 
alltraits <- pheno[, -c(30:80)]
alltraits <- alltraits[, -c(1, 3, 4:11, 19, 23, 28)]


# renaming and cleaning variables - move to other file??? 

# age
alltraits <- alltraits |> mutate(characteristics_ch1.2 = as.numeric(sub("age_years:\\s*", "", characteristics_ch1.2))) |> rename(age = characteristics_ch1.2)

# ER status
alltraits <- alltraits |> mutate(characteristics_ch1.3 = case_when(characteristics_ch1.3 == "er_status_ihc: P" ~ 1,
                                                                   characteristics_ch1.3 == "er_status_ihc: N" ~ 0,                                                              characteristics_ch1.3 %in% c("er_status_ihc: I", "er_status_ihc: NA") ~ NA_real_)) |> rename(ER_status = characteristics_ch1.3)

# PR status
alltraits <- alltraits |> mutate(characteristics_ch1.4 = case_when(characteristics_ch1.4 == "pr_status_ihc: P" ~ 1,
                                                                   characteristics_ch1.4 == "pr_status_ihc: N" ~ 0,
                                                                   characteristics_ch1.4 %in% c("pr_status_ihc: I", "pr_status_ihc: NA") ~ NA_real_)) |> rename(PR_status = characteristics_ch1.4)
# HER2 status
alltraits <- alltraits |> mutate(characteristics_ch1.5 = case_when(characteristics_ch1.5 == "her2_status: P" ~ 1,
                                                                   characteristics_ch1.5 == "her2_status: N" ~ 0,
                                                                   characteristics_ch1.5 %in% c("her2_status: I", "her2_status: NA") ~ NA_real_)) |> rename(HER2_status = characteristics_ch1.5)
# er status
alltraits <- alltraits |> mutate(characteristics_ch1.6 = case_when(characteristics_ch1.6 == "er_status_ihc_esr1_for indeterminate: P" ~ 1,
                                                                   characteristics_ch1.6 == "er_status_ihc_esr1_for indeterminate: N" ~ 0)) |> rename(indeterminate_ER_status = characteristics_ch1.6)
# clinical tumor stage 
alltraits <- alltraits |> mutate(characteristics_ch1.7 = case_when(characteristics_ch1.7 == "clinical_t_stage: T0" ~ 0,
                                                                   characteristics_ch1.7 == "clinical_t_stage: T1" ~ 1, 
                                                                   characteristics_ch1.7 == "clinical_t_stage: T2" ~ 2,
                                                                   characteristics_ch1.7 == "clinical_t_stage: T3" ~ 3, 
                                                                   characteristics_ch1.7 == "clinical_t_stage: T4" ~ 4)) |> rename(tumor_stage = characteristics_ch1.7)
# nodal status
alltraits <- alltraits |> mutate(characteristics_ch1.8 = case_when(characteristics_ch1.8 == "clinical_nodal_status: N0" ~ 0,
                                                                   characteristics_ch1.8 == "clinical_nodal_status: N1" ~ 1, 
                                                                   characteristics_ch1.8 == "clinical_nodal_status: N2" ~ 2, 
                                                                   characteristics_ch1.8 == "clinical_nodal_status: N3" ~ 3)) |> rename(nodal_status = characteristics_ch1.8)

# grade
alltraits <- alltraits |> mutate(characteristics_ch1.10 = case_when(characteristics_ch1.10 == "grade: 1" ~ 1,
                                                                    characteristics_ch1.10 == "grade: 2" ~ 2, 
                                                                    characteristics_ch1.10 == "grade: 3" ~ 3, 
                                                                    characteristics_ch1.10 == "grade: 4=Indeterminate" ~ 4,
                                                                    characteristics_ch1.10 == "grade: NA" ~ NA_real_)) |> rename(grade = characteristics_ch1.10)
# RD/pCR 
alltraits <- alltraits |> mutate(characteristics_ch1.11 = case_when(characteristics_ch1.11 == "pathologic_response_pcr_rd: RD" ~ 1,
                                                                    characteristics_ch1.11 == "pathologic_response_pcr_rd: pCR" ~ 0)) |> rename(pathologic_response = characteristics_ch1.11)

# rcb class
alltraits <- alltraits |> mutate(characteristics_ch1.12 = case_when(characteristics_ch1.12 == "pathologic_response_rcb_class: RCB-0/I" ~ 0,
                                                                    characteristics_ch1.12 == "pathologic_response_rcb_class: RCB-II" ~ 2, 
                                                                    characteristics_ch1.12 == "pathologic_response_rcb_class: RCB-III" ~ 3, 
                                                                    characteristics_ch1.12 == "pathologic_response_rcb_class: NA" ~ NA_real_)) |> rename(pathologic_response_rcb_class = characteristics_ch1.12)
# drfs 
alltraits <- alltraits |> mutate(characteristics_ch1.14 = as.numeric(sub("drfs_even_time_years:\\s*", "", characteristics_ch1.14))) |> rename(drfs = characteristics_ch1.14)

# esr1_status 
alltraits <- alltraits |> mutate(characteristics_ch1.15 = case_when(characteristics_ch1.15 == "esr1_status: P" ~ 1,
                                                                    characteristics_ch1.15 == "esr1_status: N" ~ 0)) |> rename(esr1_status = characteristics_ch1.15)

# erbb2_staus/HER2_status
alltraits <- alltraits |> mutate(characteristics_ch1.16 = case_when(characteristics_ch1.16 == "erbb2_status: P" ~ 1,
                                                                    characteristics_ch1.16 == "erbb2_status: N" ~ 0)) |> rename(erbb2_status = characteristics_ch1.16)
# set class
alltraits <- alltraits |> mutate(characteristics_ch1.17 = case_when(characteristics_ch1.17 == "set_class: SET-High" ~ 2,
                                                                    characteristics_ch1.17 == "set_class: SET-Int" ~ 1, 
                                                                    characteristics_ch1.17 == "set_class: SET-Low" ~ 0)) |> rename(set_class = characteristics_ch1.17)
# ggi class
alltraits <- alltraits |> mutate(characteristics_ch1.19 = case_when(characteristics_ch1.19 == "ggi_class: Low" ~ 0,
                                                                    characteristics_ch1.19 == "ggi_class: High" ~ 1)) |> rename(ggi_class = characteristics_ch1.19)

alltraits <- drop_na(alltraits)


# construct R object with the processed data
processed <- list(
  expression_raw=t(expr_top),
  response=response,
  phenotype=alltraits,
  genes=top_genes,
  metadata=pheno
)

saveRDS(processed, output_file)

### Data Quality Control Checks ###

# checking distribution of data
hist(expr, main = "Gene Expression Data", xlab = "Expression") # before removing low-variance 
hist(expr_top, main = "Gene Expression Data", xlab = "Expression") # after removing low-variance

# densities across samples
plotDensities(
  expr_top,
  legend = FALSE,
  main = "Expression density across samples"
)

# checking for missing entries, entries with weights below a threshold and zero-variance genes
gsg <- goodSamplesGenes(expr_top, verbose=3) 
good_samples <- gsg$allOK  
good_samples # TRUE if all samples are 'good'

# checking that the samples are in the same order for expr_top and pheno
same_samples <- identical(rownames(t(expr_top)), rownames(pheno))
same_samples # TRUE if the samples are in the same order for expression table and phenotype table
