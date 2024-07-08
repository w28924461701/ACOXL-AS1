# Ensure necessary packages are installed
if (!require("limma")) install.packages("limma")
if (!require("estimate")) {
  rforge <- "http://r-forge.r-project.org"
  install.packages("estimate", repos = rforge, dependencies = TRUE)
}
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("ggpubr")) install.packages("ggpubr")
if (!require("ggExtra")) install.packages("ggExtra")

# Load the required libraries
library(limma)
library(estimate)
library(ggplot2)
library(ggpubr)
library(ggExtra)

# Set working directory for estimate score calculation
setwd("D:\\panCancer\\19.estimate")

# Read files in the directory that match the pattern
files <- dir()
files <- grep("^symbol.", files, value = TRUE)

outTab <- data.frame()
for (i in files) {
  # Read file and process data
  CancerType <- gsub("symbol\\.|\\.txt", "", i)
  rt <- read.table(i, sep = "\t", header = TRUE, check.names = FALSE)
  
  # If a gene occupies multiple rows, take the mean
  rt <- as.matrix(rt)
  rownames(rt) <- rt[, 1]
  exp <- rt[, -1]
  dimnames <- list(rownames(exp), colnames(exp))
  data <- matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), dimnames = dimnames)
  data <- avereps(data)
  
  # Remove normal samples, keep only tumor samples
  group <- sapply(strsplit(colnames(data), "\\-"), "[", 4)
  group <- sapply(strsplit(group, ""), "[", 1)
  group <- gsub("2", "1", group)
  data <- data[, group == 0]
  out <- data[rowMeans(data) > 0, ]
  out <- rbind(ID = colnames(out), out)
  
  # Write processed matrix to file
  write.table(out, file = "uniq.symbol.txt", sep = "\t", quote = FALSE, col.names = FALSE)
  
  # Run estimate package
  filterCommonGenes(input.f = "uniq.symbol.txt", output.f = "commonGenes.gct", id = "GeneSymbol")
  estimateScore(input.ds = "commonGenes.gct", output.ds = "estimateScore.gct")
  
  # Read and process estimate scores
  scores <- read.table("estimateScore.gct", skip = 2, header = TRUE)
  rownames(scores) <- scores[, 1]
  scores <- t(scores[, -c(1, 2)])
  rownames(scores) <- gsub("\\.", "\\-", rownames(scores))
  outTab <- rbind(outTab, cbind(scores, CancerType))
  
  # Clean up intermediate files
  file.remove("commonGenes.gct")
  file.remove("estimateScore.gct")
  file.remove("uniq.symbol.txt")
}

# Write final output
out <- cbind(ID = row.names(outTab), outTab)
write.table(out, file = "estimateScores.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Set working directory for correlation analysis
setwd("D:\\panCancer\\20.estimateCor")

# Define p-value filter threshold
pFilter <- 0.05

# Read expression and TME files
exp <- read.table("singleGeneExp.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
gene <- colnames(exp)[1]
TME <- read.table("estimateScores.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

# Remove normal samples
group <- sapply(strsplit(row.names(exp), "\\-"), "[", 4)
group <- sapply(strsplit(group, ""), "[", 1)
group <- gsub("2", "1", group)
exp <- exp[group == 0, ]

# Find common samples
sameSample <- intersect(row.names(TME), row.names(exp))
TME <- TME[sameSample, ]
exp <- exp[sameSample, ]
exp[,"CancerType"] <- as.factor(exp[,"CancerType"])

# Correlation analysis
outTab <- data.frame()
for (i in levels(exp[,"CancerType"])) {
  exp1 <- exp[exp$CancerType == i, ]
  TME1 <- TME[TME$CancerType == i, ]
  y <- as.numeric(exp1[, 1])
  outVector <- data.frame(CancerType = i, Gene = gene)
  
  for (j in colnames(TME1)[1:2]) {
    x <- as.numeric(TME1[, j])
    df1 <- as.data.frame(cbind(x, y))
    corT <- cor.test(x, y, method = "spearman")
    cor <- corT$estimate
    pValue <- corT$p.value
    outVector <- cbind(outVector, pValue)
    
    p1 <- ggplot(df1, aes(x, y)) + 
      xlab(j) + ylab(gene) + 
      ggtitle(paste0("Cancer: ", i)) + theme(title = element_text(size = 10)) +
      geom_point() + geom_smooth(method = "lm") + theme_bw() +
      stat_cor(method = 'spearman', aes(x = x, y = y))
    p2 <- ggMarginal(p1, type = "density", xparams = list(fill = "orange"), yparams = list(fill = "blue"))
    
    if (pValue < pFilter) {
      pdf(file = paste0("estimateCor.", i, "_", j, ".pdf"), width = 5, height = 5)
      print(p2)
      dev.off()
    }
  }
  outTab <- rbind(outTab, outVector)
}

# Define column names and write results to file
colNames <- c("CancerType", "Gene", colnames(TME)[1:2])
colnames(outTab) <- colNames
write.table(outTab, file = "estimateCor.result.txt", sep = "\t", row.names = FALSE, quote = FALSE)
