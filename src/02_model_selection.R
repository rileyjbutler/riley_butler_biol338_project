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

# load data
data <- readRDS(wgcna_path)

clinical_factors <- data$pheno
eigengenes <- data$eigengenes

# ensure clinical factors and eigengene samples are in the same order
stopifnot(identical(clinical_factors$geo_accession, rownames(eigengenes)))

# dropping predictors that are too related to the outcome (pathologic response)
clinical_factors <- clinical_factors |> select(-geo_accession, -drfs, -grade, -erbb2_status, -indeterminate_ER_status)

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

# LOGISTIC REGRESSION # 
# first I will construct 3 general logistic regression models - 1 with only eigegenes, 1 with clinical factors and 1 combined model. 
# I run some model diagnostics by checking that none of the predictors are co-linear (VIF <=10)

# full model
logistic1 <- glm(pathologic_response ~ MEblue + MEtan + MEgreenyellow + MEblack + MEgreen + MEpink + MEsalmon + MEmagenta + MEred + MEpurple + MEturquoise + ER_status + tumor_stage  + age + PR_status + nodal_status, data=model_data, family=binomial)

# ME model
logistic2 <- glm(pathologic_response ~ MEblue + MEtan + MEgreenyellow + MEblack + MEgreen + MEpink + MEsalmon + MEmagenta + MEred + MEpurple + MEturquoise, data=model_data, family=binomial)

# clinical factor model
logistic3 <- glm(pathologic_response ~ ER_status + tumor_stage  + age + PR_status + nodal_status, data=model_data, family=binomial)

# summary statistics 
l1 <- summary(logistic1)
l1
l2 <- summary(logistic2)
l2
l3 <- summary(logistic3)
l3

# check collinearity shows no high collinearity (VIF >= 10)
check_collinearity(logistic1)
check_collinearity(logistic2)
check_collinearity(logistic3)

# using caret for model comparison 
MEpredictors <- model_data |> select(-age, -ER_status, -PR_status, -ggi_class, -HER2_status, -tumor_stage, -nodal_status, -esr1_status, -set_class)

# construct logistic regression model with MEs
set.seed(123)
ME_logistic_model <- caret::train(pathologic_response ~ .,
                              data = MEpredictors,
                              method = "glm", # using logistic regression 
                              family=binomial,
                              metric = "ROC", # use ROC as performance measure 
                              trControl = control)

# construct logistic regression model with clinical factors
CF_predictors <- model_data |> select(-MEblue, -MEtan, -MEgreenyellow, -MEblack, -MEgreen, -MEpink, -MEsalmon, -MEmagenta, -MEred, -MEpurple, -MEturquoise)

# running logistic regression model 
set.seed(123)
CF_logistic_model <- caret::train(pathologic_response ~ .,
                                  data = CF_predictors,
                                  method = "glm", # using logistic regression 
                                  family=binomial,
                                  metric = "ROC", # use ROC as performance measure 
                                  trControl = control)

# construct full logistic regression model 

# running full logistic regression model 
set.seed(123)
full_logistic_model <- caret::train(pathologic_response ~ .,
                                  data = model_data,
                                  method = "glm", # using logistic regression 
                                  family=binomial,
                                  metric = "ROC", # use ROC as performance measure
                                  trControl = control)
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
importance <- varImp(random_forest_model, scale=FALSE)
plot(importance)

# comparing models using AUC, Specificity and Sensitivity 
results <- resamples(list(glmnet = elastic_model, MElogistic = ME_logistic_model, CFlogistic = CF_logistic_model, fullLogistic = full_logistic_model, SVM_model=SVM_model, randomForest = random_forest_model))
summary(results)
bwplot(results)

# construct R object with the models
models <- list(
  elastic_model=elastic_model,
  ME_logistic_model=ME_logistic_model, 
  CF_logistic_model=CF_logistic_model, 
  full_logistic_model=full_logistic_model,
  SVM_model=SVM_model, 
  random_forest_model=random_forest_model,
  model_data = model_data
)

saveRDS(models, output_file)
