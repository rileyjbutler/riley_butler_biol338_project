# set-up 
if (!dir.exists("results")) {
  dir.create("results", recursive = TRUE)
}

wgcna_path <- "data/wgcna_results.rds"
if (!file.exists(wgcna_path)) {
  system("Rscript 01_wgcna.R")
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

# dropping predictors that are too related to the outcome (pathologic response)
clinical_factors <- clinical_factors |> select(-geo_accession, -pathologic_response_rcb_class, -drfs)
eigengenes <- data$eigengenes

# create dataframe
model_data <- data.frame(clinical_factors,eigengenes)

# swapping variables - maybe do this early on? 
model_data$pathologic_response <- ifelse(model_data$pathologic_response == 0, 1, 0)

# x and y set-up for models
x <- model_data |> select(-pathologic_response)
x <- as.matrix(x)
y <- model_data$pathologic_response

# using caret to choose parameters (alpha and lambda) and fit elastic net logistic regression model 

control <- trainControl(method = "repeatedcv", # uses repeated cross-validation with GSE25055 training set (model_data)
                        number = 5, # number of folds 
                        repeats = 5, # number of repeated folds
                        verboseIter = TRUE, # print training log 
                        classProbs = TRUE, # calculated class probabilities for each re-sample 
                        summaryFunction = twoClassSummary, # to compute performance metrics 
                        savePredictions = "final") # indicator of how much of the hold-out predictions for each resample should be saved 

# changing outcome to be factor (pCR or RD)
model_data$pathologic_response <- factor(model_data$pathologic_response, levels = c(1, 0), labels = c("pCR", "RD"))

# specifying range of alpha and lambda values to try. For reproducibility. 
grid <- expand.grid(alpha = seq(0, 1, by = 0.1), lambda = 10^seq(-4, 0, length.out = 50))

# constructing the elastic net logistic regression model
elastic_model <- train(pathologic_response ~ .,
                      data = model_data,
                      method = "glmnet", # using elastic-net logistic regression 
                      preProcess = c("center", "scale"),
                      tuneGrid = grid,
                      metric = "ROC", # use ROC as performance measure 
                      trControl = control) # use information from the control 

elastic_model # shows the average resampled estimates of performance for different lambda and alpha values 

# shows model performance (ROC-AUC) for different alpha and lambda values
ggplot(elastic_model)

# shows the coefficients used for the final model. Note that some correlated coefficients in the model are shrunk. 
coef(elastic_model$finalModel, s = elastic_model$bestTune$lambda)

# predictions for training data, using the elastic model. Show the probability of being pCR or RD. 
pred_probs <- predict(elastic_model, newdata = model_data, type="prob")
pred_probs

# NOTE: because none of the probabilities exceed 0.5, all samples are classified as RD. 
# To improve sensitivity, the threshold is changed using Youden 

# averaging predictions over k-folds
patient_predictions <- elastic_model$pred |> group_by(rowIndex) |> summarise(obs = first(obs), pCR_prob = mean(pCR), .groups = "drop")

# takes actual outcome for each patient and the models predicted probability, calculates how sensitivity and specificity changes across different cutoffs
roc_obj <- roc(response = patient_predictions$obs, predictor = patient_predictions$pCR_prob, levels = c("RD", "pCR"), direction = "<")

# plotting the ROC curve
plot(roc_obj)

# finding the best threshold by balancing sensitivity and specificity by using Youden-selected threshold. 
# which maximizes: sensitivity + specificity - 1
best_threshold <- coords(roc_obj, x = "best", best.method = "youden", ret = c("threshold", "sensitivity", "specificity"))

best_threshold

# assign optimal cutoff value 
cutoff <- as.numeric(best_threshold["threshold"])

patient_predictions$pred <- factor(ifelse(patient_predictions$pCR_prob >= cutoff, "pCR", "RD"))
patient_predictions$obs <- factor(patient_predictions$obs,levels = c("pCR", "RD"))

# make a confusion matrix with the new cutoff 
confusionMatrix(patient_predictions$pred, patient_predictions$obs, positive = "pCR")
