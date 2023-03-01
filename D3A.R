 #########################
##                       ##
##  ###    #####    #    ##
##  #  #       #   # #   ##
##  #   #  #####  #####  ##
##  #  #       #  #   #  ##
##  ###    #####  #   #  ##
##                       ##
 #########################

#Set the working directory
setwd("Z:/")

#Open all necessary libraries
library(DESeq2)
library(edgeR)
library(gplots)
library(grid)
library(NOISeq)
library(readr)
library(readxl)
library(VennDiagram)

#Upload all necessary inputs and if NAs in countData change to 0s
countData <- read.delim("countData.csv", header = T, sep = ";")
countData[is.na(countData)] <- 0 #Just if we have NAs in countData
metaData <- read.delim("metaData.csv", header = T, sep = ";")

#############
### edgeR ###
#############

#Prepare counts and metadata to create the edgeR object (dgList)
Counts <- data.frame(countData, row.names = 1)
sampleType <- metaData$Disease
dgList <- DGEList(counts=Counts, genes=rownames(Counts), group = factor(sampleType))

#Filtering of no/low expression genes
countsPerMillion <- cpm(dgList)
countCheck <- countsPerMillion > 1
keep <- which(rowSums(countCheck) >= 3)
dgList <- dgList[keep, ]

#Calculate normalized factors with TMM (Trimmed Mean of M-values) method
dgList <- calcNormFactors(dgList, method="TMM")

#MDS plot to see how the samples are separated
plotMDS(dgList, col=as.numeric(dgList$samples$group))
legend("bottomleft", as.character(unique(dgList$samples$group)), col=1:3, pch=20)

#Design the comparison matrix
designMat <- model.matrix(~dgList$samples$group)

#Calculate the ...
dgList <- estimateGLMCommonDisp(dgList, design=designMat)
dgList <- estimateGLMTrendedDisp(dgList, design=designMat)
dgList <- estimateGLMTagwiseDisp(dgList, design=designMat)
plotBCV(dgList)

#Perform the differential expression analysis
fit <- glmFit(dgList, designMat)
lrt <- glmLRT(fit)
edgeR_res_1 <- data.frame(lrt)
n = nrow(edgeR_res_1)
res_edgeR <- topTags(lrt, adjust.method = "fdr", n = n)
res_edgeR <- data.frame(res_edgeR)
write.table(res_edgeR, "res_edgeR.txt", sep = "\t", row.names = F)

##############
### DESeq2 ###
##############

#Create DESeq2 object (dds) specifying the variable to compare the samples
dds <- DESeqDataSetFromMatrix(countData = countData, colData = metaData, design = ~Disease, tidy = T)
dds <- DESeq(dds)

#Perform DEA and save the raw results
res_DESeq2 <- results(dds, contrast = c("Disease", "Case", "Control"), pAdjustMethod = "fdr")
res_DESeq2 <- data.frame(res_DESeq2)
res_DESeq2 <- tibble::rownames_to_column(res_DESeq2, "genes")
write.table(res_DESeq2, "res_DESeq2.txt", sep = "\t", row.names = T)

##############
### NOISeq ###
##############

#Create NOISeq object (mydata) and normalize the data with TMM (Trimmed Mean of M-values) method
samples <- data.frame(metaData$Disease)
mydata <- readData(data = Counts, factors = samples)
exprs_data <- assayData(mydata)$exprs
myTMM = tmm(exprs_data)

#Filtering of no/low expression genes
meanexpression = rowMeans(myTMM)
threshold = quantile(meanexpression, 0.33)
myTMM = myTMM[which(meanexpression >= threshold), ]

#Specify the normalization method and the factor to perform the comparison between groups of
res_NOISeq = noiseqbio(mydata, norm = "tmm", factor = "metaData.Disease", filter = 0)
res_NOISeq = degenes(res_NOISeq, q = 0, M = NULL)
res_NOISeq$FDR <- 1 - res_NOISeq$prob
res_NOISeq <- tibble::rownames_to_column(res_NOISeq, "genes")
write.table(res_NOISeq, "res_NOISeq.txt", sep = "\t", row.names = T)

#####################
### Venn diagrams ###
#####################

#Filter each result to get only the significant
#edgeR
res_edgeR$AbslogFC = abs(res_edgeR$logFC)
sig_res_edgeR <- res_edgeR[which(res_edgeR$AbslogFC > 2 & res_edgeR$FDR < 0.01),]
sig_res_edgeR <- as.data.frame(sig_res_edgeR)
#DESeq2
res_DESeq2$Abslog2FC = abs(res_DESeq2$log2FoldChange)
sig_res_DESeq2 <- res_DESeq2[which(res_DESeq2$Abslog2FC > 2 & res_DESeq2$padj < 0.01),]
sig_res_DESeq2 <- as.data.frame(sig_res_DESeq2)
#NOISeq
res_NOISeq$Abslog2FC = abs(res_NOISeq$log2FC)
sig_res_NOISeq <- res_NOISeq[which(res_NOISeq$Abslog2FC > 2 & res_NOISeq$FDR < 0.01),]
sig_res_NOISeq <- as.data.frame(sig_res_NOISeq)

#Venn diagrams
sig_names_edgeR <- sig_res_edgeR[,1]
sig_names_DESeq2 <- sig_res_DESeq2[,1]
sig_names_NOISeq <- sig_res_NOISeq[,1]
Venn_input <- list(sig_names_edgeR, sig_names_DESeq2, sig_names_NOISeq)
tmp <- venn(Venn_input)
Intersections <- attr(tmp, "intersections")

#We can see the gene names of the Venn diagram intersections
View(Intersections)
Intersections[["B:D"]]

#Nice Venn diagrams
display_venn <- function(x, ...){
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}
display_venn(
  Venn_input,
  category.names = c("edgeR" , "DESeq2", "NOISeq"),
  fill = c("#974939", "#E69F00", "#56B4E9"),
  lwd = 2,
  lty = 'blank',
  cat.fontface = "bold",
  main = "Venn diagram for DEA result comparison")
