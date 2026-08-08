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

# filtering for only useful variables 
alltraits <- pheno
alltraits <- pheno[, -c(30:80)]
alltraits <- alltraits[, -c(1, 3, 4:11, 19, 23, 28)]


# renaming and cleaning variables 

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

# matching trait samples to expression samples 
samples <- rownames(expression)
traitrows <- match(samples, alltraits$geo_accession)
datatraits <- alltraits[traitrows, -1]
rownames(datatraits) <- alltraits[traitrows, 1]

# calculating module-trait correlation
nGenes <- ncol(expression) # 10000
nSamples <- nrow(expression) # 306

module_trait_corr <- WGCNA::cor(MEs, datatraits, use='p')
module_trait_pvalue <- corPvalueStudent(module_trait_corr, nSamples)

# forming a heatmap 
textMatrix = paste(signif(module_trait_corr, 2), "\n(",
                   signif(module_trait_pvalue, 1), ")", sep = "");
dim(textMatrix) = dim(module_trait_corr)
par(mar = c(6, 8.5, 3, 1))

# display the correlation values within a heatmap plot
labeledHeatmap(Matrix = module_trait_corr,
               xLabels = names(datatraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.4,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
