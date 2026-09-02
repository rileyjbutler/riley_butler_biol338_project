# This file is used to carry out a Weighted Gene Co-expression Networks Analysis using the processed data from 00_preprocessing
validation = TRUE

# load libraries 
library("WGCNA")
library(tidyverse)
library(pheatmap)

# set-up 
if (!dir.exists("results")) {
  dir.create("results", recursive = TRUE)
}
if (validation == FALSE) {
  processed_path <- "data/processed_data.rds"
  output_file <- file.path(data_dir, "wgcna_results.rds")
  if (!file.exists(processed_path)) {
    system("Rscript 00_preprocessing.R")
  }
} else if (validation == TRUE) {
  processed_path <- "data/validation_processed_data.rds"
  output_file <- file.path(data_dir, "validation_swgcna_results.rds")
  if (!file.exists(processed_path)) {
    system("Rscript 00_preprocessing.R")
  }
}

data_dir <- "data" 
if (!dir.exists(data_dir)) { 
  dir.create(data_dir, recursive=TRUE)
}

# retrieving processed data 
processed <- readRDS(processed_path)
expression <- processed$expression_raw
response <- factor(processed$response, levels = c("RD", "pCR"))
genes <- processed$genes
pheno <- processed$phenotype

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
# CHANGE FOR VALIDATION SET???
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

# module merging 
ME_dis <- 1-cor(MElist$eigengenes, use="complete") # calculate eigengene dissimilarity 
METree <- hclust(as.dist(ME_dis), method="average") # clustering eigengenes 
par(mar=c(0,4,2,0))
par(cex=0.6)
plot(METree)
abline(h=.25, col="red")

merge <- mergeCloseModules(expression, ModuleColors, cutHeight=0.25)

mergedColors <- merge$colors
mergedMEs = merge$newMEs

# plot merged vs unmerged modules for comparison 
plotDendroAndColors(geneTree, cbind(ModuleColors, mergedColors), 
                    c("Original Module", "Merged Module"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors for original and merged modules")

# module-trait matching

# matching trait samples to expression samples 
common <- intersect(rownames(expression), pheno$geo_accession)
expression2 <- expression[common, , drop=FALSE]
datatraits <- pheno[match(common, pheno$geo_accession), , drop=FALSE]
rownames(datatraits) <- datatraits$geo_accession

#samples <- rownames(expression)
#traitrows <- match(samples, pheno$geo_accession)
#datatraits <- pheno[traitrows, -1]
#rownames(pheno2) <- pheno[traitrows, 1]


# matching with samples from metadata with MEs
common <- intersect(rownames(mergedMEs), rownames(datatraits))
MEs2 <- mergedMEs[common, , drop = FALSE]
datatraits2 <- datatraits[common, , drop = FALSE]
identical(rownames(MEs2), rownames(datatraits2))

# calculating module-trait correlation

nGenes <- ncol(expression2) # 10000
nSamples <- nrow(expression2) # 306

traits_numeric <- datatraits2[, setdiff(names(datatraits2), "geo_accession")] # get rid of character variable
module_trait_corr <- WGCNA::cor(MEs2, traits_numeric, use='p')
module_trait_pvalue <- corPvalueStudent(module_trait_corr, nSamples)

# forming a heatmap 
textMatrix = paste(signif(module_trait_corr, 2), "\n(",
                   signif(module_trait_pvalue, 1), ")", sep = "");
dim(textMatrix) = dim(module_trait_corr)
par(mar = c(6, 8.5, 3, 1))


# display the correlation values within a heatmap plot
labeledHeatmap(Matrix = module_trait_corr,
               xLabels = names(traits_numeric),
               yLabels = names(MEs2),
               ySymbols = names(MEs2),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.4,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))


# exploration for module eigengenes 

queryModuleColor <- "yellow" # using yellow module 

queryModuleExpression <- expression[,which(ModuleColors==queryModuleColor)]

# a heatmap of sample clustering for genes inside the yellow module based on their association with pCR (0 or 1)
heatmap1 <- heatmap(as.matrix(queryModuleExpression), RowSideColors=c("red", "black")[as.numeric(as.factor(datatraits2$pathologic_response))])

# combining pathologic response and ER status

annotation_row <- data.frame(Response = datatraits$pathologic_response, ER_status = factor(datatraits$ER_status, levels = factor(c(0, 1)), labels = c("ER-", "ER+")))

rownames(annotation_row) <- rownames(datatraits)

# heat map of the yellow module 
pheatmap(queryModuleExpression, scale = "column", annotation_row = annotation_row, show_rownames = FALSE, main = "Yellow module")

# ER status is shown in red and blue, showing a distinct cluster for ER+ (red) and ER- (blue) in the yellow module 
# pCR shown in green/white, showing no distinct clustering in the yellow module

# view as PCA 
queryModulePCA<-prcomp(queryModuleExpression,scale=TRUE,center=TRUE)

summary(queryModulePCA)$importance[,2]

# plot - seems to have roughly equal split in direction across PC1 and PC2
biplot(queryModulePCA, cex=0.4)

# find Eigengene 
thisEigenGene <- MEs[,which(colnames(MEs)==paste0("ME", queryModuleColor))]

# find first principle component
thisPC1 <- queryModulePCA$x[,1]

# plot - some separation of pCR/RD cases, lots more separation for ER status 
plot(thisEigenGene, thisPC1, col=c("red", "black")[as.factor(pheno$pathologic_response)])
plot(thisEigenGene, thisPC1, col=c("red", "black")[as.factor(pheno$ER_status)])

wgcna_results <- list(
  eigengenes=MEs2,
  pheno=pheno
)

saveRDS(wgcna_results, output_file)
