# This file is used to carry out a Weighted Gene Co-expression Networks Analysis using the processed data from 00_preprocessing

# load libraries 
library("WGCNA")
library(tidyverse)

# set-up 
if (!dir.exists("results")) {
  dir.create("results", recursive = TRUE)
}
processed_path <- "data/processed_data.rds"
if (!file.exists(processed_path)) {
  system("Rscript 00_preprocessing.R")
}

# retrieving processed data 
processed <- readRDS(processed_path)
expression <- processed$expression_raw
response <- factor(processed$response, levels = c("RD", "pCR"))
genes <- processed$genes
pheno <- processed$pheno

# identifying outlier samples 
sampleTree <- hclust(dist(expression), method="average")

par(cex=0.6);
par(mar=c(0,4,2,0))

plot(sampleTree, main="sample clustering to detect outliers", sub="", xlab="", cex.lab=1.5, cex.axis=1.5, cex.main=2)
# there appears to be no clear outlier samples

# choosing soft power threshold 
spt <- pickSoftThreshold(expression)

# plotting soft power thresholds
par(mar=c(1,1,1,1))
plot(spt$fitIndices[,1],spt$fitIndices[,2],
xlab="Soft Threshold", ylab="Scale Free Topology Model Fit", type="n", 
main=paste("Scale independece"))
text(spt$fitIndices[,1], spt$fitIndices[,2], col="red")
abline(h=0.80, col="red")

# plotting mean connectivity
par(mar=c(5,5,3,1))
plot(spt$fitIndices[,1], spt$fitIndices[,5],
     xlab="Soft Threshold", ylab="Mean Connectivity", type="n", xaxt="n", 
     main=paste("Mean connectivity"))
axis(side=1, at=spt$fitIndices[,1], labels=spt$fitIndices[,1])
text(spt$fitIndices[,1], spt$fitIndices[,5], labels=spt$fitIndices[,1], col="red")

# choose soft power = 3
SoftPower <- 3 

# construct adjacency matrix
adjacency <- adjacency(expression, power=SoftPower)

# calculate TOM dissimilarity
TOM <- TOMsimilarity(adjacency)
TOM_dis <- 1-TOM

# plot gene tree
geneTree <- hclust(as.dist(TOM_dis), method='average')

sizeGrWindow(12,9)
plot(geneTree, xlab='', sub='', main="Gene clustering on TOM-based dissimilarity", 
     labels=FALSE, hang=0.04)

# find modules and assign module colours
Modules <- cutreeDynamic(dendro=geneTree, distM=TOM_dis, deepSplit=2, pamRespectsDendro=FALSE, minClusterSize=30)
table(Modules)

ModuleColors <- labels2colors(Modules)
table(ModuleColors)

# plot dendrogram with module colours
plotDendroAndColors(geneTree, ModuleColors, "Module", 
                    dendroLabels=FALSE, hang=0.03, 
                    addGuide=TRUE, guideHange=0.05, 
                    main="Gene dendrogram and module colors")

# module eigengene identification 
MElist <- moduleEigengenes(expression, colors=ModuleColors)
MEs <- MElist$eigengenes 
head(MEs)

# module merging 
ME_dis <- 1-cor(MElist$eigengenes, use="complete") # calculate eigengene dissimilarity 
METree <- hclust(as.dist(ME_dis), method="average") # clustering eigengenes 
par(mar=c(0,4,2,0))
par(cex=0.6)
plot(METree)
abline(h=.25, col="red")

# module-trait matching
# encode pCR and RD as binary 

###################################################################
traits <-data.frame(pCR = ifelse(processed$response == "pCR", 
                        1, 
                        ifelse(processed$response == "RD", 0, NA)))

rownames(traits) <- rownames(pheno)

# define number of genes and samples 
nGenes <- ncol(expression) # 10000
nSamples <- nrow(expression) # 306

module_trait_corr <- cor(MEs, traits, use='p')
module_trait_pvalue <- corPvalueStudent(module_trait_corr, nSamples)

################## USING ALL VARIABLES ###############

# filtering for only useful variables 
alltraits <- pheno
alltraits <- pheno[, -c(30:80)]
alltraits <- alltraits[, -c(1, 3, 4:11)]

# renaming and cleaning variables 

# age
alltraits <- alltraits |> mutate(characteristics_ch1.2 = as.numeric(sub("age_years:\\s*", "", characteristics_ch1.2))) |> rename(age = characteristics_ch1.2)

# ER status
alltraits <- alltraits |> mutate(characteristics_ch1.3 = case_when(characteristics_ch1.3 == "er_status_ihc: P" ~ 1,
                                                                   characteristics_ch1.3 == "er_status_ihc: N" ~ 0,
                                                                   characteristics_ch1.3 %in% c("er_status_ihc: I", "er_status_ihc: NA") ~ NA_real_)) |> rename(ER_status = characteristics_ch1.3)

# PR status"

alltraits <- alltraits |> mutate(characteristics_ch1.4 = case_when(characteristics_ch1.4 == "pr_status_ihc: P" ~ 1,
                                                                   characteristics_ch1.4 == "pr_status_ihc: N" ~ 0,
                                                                   characteristics_ch1.4 %in% c("pr_status_ihc: I", "pr_status_ihc: NA") ~ NA_real_)) |> rename(PR_status = characteristics_ch1.4)
# use unique(alltraits$characteristics_ch1.4)

alltraits <- alltraits |> mutate(characteristics_ch1.4 = case_when(characteristics_ch1.4 == "pr_status_ihc: P" ~ 1,
                                                                   characteristics_ch1.4 == "pr_status_ihc: N" ~ 0,
                                                                   characteristics_ch1.4 %in% c("pr_status_ihc: I", "pr_status_ihc: NA") ~ NA_real_)) |> rename(PR_status = characteristics_ch1.4)