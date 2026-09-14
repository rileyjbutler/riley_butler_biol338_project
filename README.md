# Riley Butler BIOL338 Project 
## Using WGCNA to Predict The Reponse of Taxane-Anthracycline Chemotherapy for HER2-negative Breast Cancer 

### Summary  
This repository uses Weighted Gene Co-expression Network Analysis (WGCNA) and clinical factors to construct a model that attempts to predict patient chemotherapy response. The goal of this project is to improve patient-treatment matching for HER2-negative cancer patients as well as explore how genetic factors can influence patients' chemotherapy response. The training data is processed and run through WGCNA; the results are used to run different model types, and the model with the best performance is used on the test set. The chosen model is run through performance and diagnostics tests. Finally, the test set is processed, and the eigengenes from the training set are projected onto the test set. Then, predictions are made using the test set and the constructed model. 

### Dataset 
The repository uses GEOquery to access two datasets: GSE25055 and GSE25065. The GSE25055 dataset is used as a training set for WGCNA and model construction, whilst GSE25065 is used as a validation set. The training set contains data from 310 patients with (primarily) HER2-negative invasive breast cancer. The validation set contains data from 180 patients. The dataset includes microarray-based gene expression tests from tumour biopsies before surgery and chemotherapy treatment. The clinical metadata includes many factors that might contribute to chemotherapy response and cancer development. The chemotherapy response is measured as RD (residual disease) or pCR (pathologic complete response). 

### Usage 

To set-up the R environment, run the following code in the RStudio console 
```
install.packages("renv")
renv::restore()
```

Open the src folder and run the files in order (00 through 04). To get the processed validation data, change validation to TRUE at the top of the 00_preprocessing file. 

'src' folder

00_preprocessing - cleans and filters the data in preparation for WGCNA. This includes removing unnecessary variables, removing genes with low variation, renaming variables and visualising results.  
01_wgcna - runs a network analysis on the processed data using a chosen soft power threshold of 3. Conducts PCA on the blue module. 
02_model_selection - constructs different models using the caret package to assess which model performs best on training data.
03_ME_projection - projects the training module eigenvalues onto the processed validation data.
04_performance - runs performance and diagnostic checks on the elastic net logistic regression model as well as predicts pCR/RD probabilities for the validation set.

'results' folder contains .png files of the figures constructed throughout the pipeline

'data' folder contains R objects of the data throughout the pipeline 

'genes.txt' contains the gene names of the genes in the blue module for enrichment analysis. These were pasted into PantherDB and annotated. 

### Contact 
If you need help using this repository or run into a bug, please contact Riley Butler at rileyjbutler@gmail.com 
