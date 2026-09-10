# load libraries 
library(tidyverse)

# script set-up
data_dir <- "data" 
if (!dir.exists(data_dir)) { 
  dir.create(data_dir, recursive=TRUE)
}

wgcna_path <- "data/wgcna_results.rds"

if (!file.exists(wgcna_path)) {
  system("Rscript 01_wgcna.R")
}

testing_path <- "data/validation_processed_data.rds"
output_file <- file.path(data_dir, "testing_MEs.rds")
processed_path <- "data/processed_data.rds"
output_file <- file.path(data_dir, "test_MEs")

# variable set-up for the loop 
processed_data <- readRDS(processed_path)
training_data <- readRDS(wgcna_path)
testing_data <- readRDS(testing_path)
test_expression <- testing_data$expression_raw

ModuleColors <- training_data$ModuleColors
module_names <- setdiff(unique(ModuleColors), "grey")
module_names <- module_names[
  paste0("ME", module_names) %in% colnames(MEs)
]

test_MEs <- matrix(NA, nrow=nrow(test_expression), ncol=length(module_names))
rownames(test_MEs) <- rownames(test_expression)
colnames(test_MEs) <- paste0("ME", module_names)
expression <- processed_data$expression_raw
MEs <- training_data$eigengenes
projection_summary <- data.frame()
  
# for each module in the list of modules
for (module_color in module_names) {
  # finding the genes in the give module (i.e. per colour)
  module_genes <- colnames(expression)[ModuleColors == module_color]
  
  # find genes that are not part of test set
  missing_genes <- setdiff(module_genes, colnames(test_expression))
  
  # find genes that are in the module and the test set 
  shared_genes <- intersect(module_genes, colnames(test_expression))
  
  # get those genes from the expression data (training)
  module_training <- expression[, shared_genes, drop = FALSE]
  
  # carry out pca for the genes in the training set (calculates the mean, SD, weight of every gene, pc1 score of training set)
  module_pca <- prcomp(module_training, center = TRUE, scale. = TRUE)
  
  # corresponding to WGCNA eignengene
  ME_name <- paste0("ME", module_color)
  
  # check that the correlation sign is correct (should be 1) by comparing to original WGCNA eigengene
  sign_correlation <- cor(
    module_pca$x[, 1], MEs[, ME_name])
  
  # change the direction so that correlation = 1
  direction <- ifelse(sign_correlation < 0, -1, 1)
  module_train_pc1 <- module_pca$x[, 1] * direction
  
  # 
  ME_scale_model <- lm(MEs[, ME_name] ~ module_train_pc1)
  reconstructed_train_ME <- predict(ME_scale_model)
  reconstruction_correlation <- cor(reconstructed_train_ME, MEs[, ME_name])
  
  # get gene subset from test set
  test <- test_expression[, shared_genes, drop = FALSE]
  
  # standardize test set using training mean and SD
  module_test_scaled <- scale(test, center = module_pca$center, scale =  module_pca$scale)
  
  # projection to give 1 value for each test sample
  test_pc1 <- c(module_test_scaled %*% module_pca$rotation[, 1]) * direction
  test_ME <- coef(ME_scale_model)[1] + coef(ME_scale_model)[2] * test_pc1
  
  test_MEs[, ME_name] <- test_ME
  
  # quality control check for each module 
  projection_summary <- rbind(
    projection_summary,
    data.frame(
      module = module_color,
      training_genes = length(module_genes),
      shared_genes = length(shared_genes),
      missing_genes = length(missing_genes),
      retained_percent =
        100 * length(shared_genes) / length(module_genes),
      correlation_with_original_ME =
        abs(sign_correlation)
    )
  )

}

saveRDS(test_MEs, output_file)

# making a gene to module table for PantherDB 
gene_module_table <- data.frame(gene = colnames(expression), module = ModuleColors) |> group_by(module) |> summarise(genes = paste(gene, collapse = ", "))

# checking the df 
head(gene_module_table)

blue <- gene_module_table |> select(module)
