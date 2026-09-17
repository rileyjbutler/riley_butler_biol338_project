# This file is used to explore different model types for predicting pCR. 
# The chosen model is the one with the best balance of specificity, sensitivity and AUC as well as the model that best aligns with the data contextually.

# script set-up 
wgcna_path <- "data/wgcna_results.rds"
data_dir <- "data" 
output_file <- file.path(data_dir, "models.rds")

if (!file.exists(wgcna_path)) {
  system("Rscript src/01_wgcna.R")
  print("Running WGCNA")
} else { 
  print("WGCNA Results Found")
}

# load libraries 
library(tidyverse)
library(performance)
library(glmnet)
library(caret)
library(pROC)
library(car)

# load data
data <- readRDS(wgcna_path)

clinical_factors <- data$pheno
eigengenes <- data$eigengenes

# ensure clinical factors and eigengene samples are in the same order
stopifnot(identical(clinical_factors$geo_accession, rownames(eigengenes)))

# dropping predictors that are too related to the outcome (pathologic response)
clinical_factors <- clinical_factors |> select(-geo_accession, -drfs, -grade, -erbb2_status, -indeterminate_ER_status, -HER2_status)

# create dataframe
model_data <- data.frame(clinical_factors,eigengenes)

# changing outcome to be factor
model_data$pathologic_response <- factor(model_data$pathologic_response, levels = c(1, 0), labels = c("pCR", "RD"))

set.seed(123)

folds <- caret::createMultiFolds(model_data$pathologic_response, k=5, times=5) 

# formatting for the model to be trained against 
control <- caret::trainControl(method = "repeatedcv", # uses repeated cross-validation with GSE25055 training set (model_data)
                               number = 5, # number of folds 
                               repeats = 5, # number of repeated folds
                               index = folds,
                               classProbs = TRUE, # calculated class probabilities for each re-sample 
                               summaryFunction = twoClassSummary, # to compute performance metrics 
                               savePredictions = "final") # indicator of how much of the hold-out predictions for each resample should be saved 


# ELASTIC NET LOGISTIC REGRESSION # 
# using caret to choose parameters (alpha and lambda) and fit elastic net logistic regression model 

# specifying range of alpha and lambda values to try. For reproducibility. 
grid <- expand.grid(alpha = seq(0, 1, by = 0.1), lambda = 10^seq(-4, 0, length.out = 50))

# constructing the elastic net logistic regression model

# running the elastic model 
set.seed(123)
elastic_model <- caret::train(pathologic_response ~ .,
                      data = model_data,
                      method = "glmnet", # using elastic-net logistic regression 
                      preProcess = c("center", "scale"),
                      tuneGrid = grid,
                      metric = "ROC", # use ROC as performance measure 
                      trControl = control) # use information from the control 

elastic_net <- elastic_model$finalModel

# lambda and alpha for best tune
elastic_net$tuneValue

# LOGISTIC REGRESSION # 
# first I will construct 3 general logistic regression models - 1 with only eigegenes, 1 with clinical factors and 1 combined model. 
# I run some model diagnostics by checking that none of the predictors are co-linear (VIF <=10)

# make predictors for ME only model
ME_names <- grep("^ME", names(model_data), value = TRUE)
MEpredictors <- model_data |> select(pathologic_response, all_of(ME_names))

# ME model
ME_logistic_model <- caret::train(pathologic_response ~ ., data = MEpredictors, method = "glm", family = binomial, metric = "ROC", trControl = control)
ME_lm <- ME_logistic_model$finalModel

# full  model
full_logistic <- caret::train(pathologic_response ~ ., data = model_data, method = "glm", family = binomial, metric = "ROC", trControl = control)
full_lm <- full_logistic$finalModel

# clinical factor model
cf_logistic <- caret::train(pathologic_response ~ ER_status + tumor_stage  + age + PR_status + nodal_status, data = model_data, method = "glm", family = binomial, metric = "ROC", trControl = control)
cf_lm <- full_logistic$finalModel

# summary statistics 
summary(ME_logistic_model)
summary(full_logistic)
summary(cf_logistic)

# checking collinearity, high collinearity  at VIF >= 10
full_vif <- vif(full_lm)
full_vif

ME_vif <- vif(ME_lm)
ME_vif 

cf_vif <- vif(cf_lm)
cf_vif

# SUPPORT VECTOR MACHINE (SVM) # 
# SVM finds the optimal hyperplan to separate data into different classes (pCR or RD). 
# The RBF Kernal SVM uses non-linear relationships to map the data into infinite dimensional space 
# The tuning cost (C) defines the influence of a single training example 

# running SVM model with all predictors
set.seed(123)
SVM_model <- caret::train(pathologic_response ~ .,
                                    data = model_data,
                                    method = "svmRadial", # using SVM
                                    family=binomial,
                                    metric = "ROC", # use ROC as performance measure 
                                    trControl = control)

# RANDOM FOREST # 
# Random forest is an ensemble techique that builds multiple decision trees and merges their outputs to improve model accuracy and stability 

# running random forest with all predictors 
set.seed(123)
random_forest_model <- caret::train(pathologic_response ~ .,
                          data = model_data,
                          method = "ranger", # using random forest 
                          metric = "ROC", # use ROC as performance measure 
                          trControl = control,
                          importance = "permutation")

# to observe the importance of each predictor - how much the model's predictive performance depends on that variable.
RF_importance <- varImp(random_forest_model, scale=FALSE)
plot(RF_importance)

EN_importance <- varImp(elastic_model, scale=FALSE)
plot(EN_importance)

SVM_importance <- varImp(SVM_model, scale=FALSE)
plot(SVM_importance)

# comparing models using AUC, Specificity and Sensitivity 
results <- resamples(list(glmnet = elastic_model, MElogistic = ME_logistic_model, CFlogistic = cf_logistic, fullLogistic = full_logistic, SVM_model=SVM_model, randomForest = random_forest_model))
summary(results)
box_whisker <- bwplot(results)

saveRDS(box_whisker, "results/box_wisker.rds")

# construct R object with the models
models <- list(
  elastic_model=elastic_model,
  ME_logistic_model=ME_logistic_model, 
  CF_logistic_model=cf_logistic, 
  full_logistic_model=full_logistic,
  SVM_model=SVM_model, 
  random_forest_model=random_forest_model,
  model_data = model_data
)

saveRDS(models, output_file)
