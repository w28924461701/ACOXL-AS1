# Ensure necessary packages are installed
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require("limma")) BiocManager::install("limma")
if (!require("ggpubr")) install.packages("ggpubr")
if (!require("stringr")) install.packages("stringr")
if (!require("tidyverse")) install.packages("tidyverse")

# Load required libraries
library(limma)
library(ggpubr)
library(stringr)
library(tidyverse)

# Set working directory
setwd("D:\\panCancer\\13.cliCor")

# Input files
file <- "singleGeneExp.txt"  # Expression data file
rt <- read.table(file, sep = "\t", header = TRUE, check.names = FALSE, row.names = 1)
cli1 <- read.table("clinical1.txt", sep = "\t", header = TRUE, check.names = FALSE, row.names = 1)

# Process clinical data
cli <- select(cli1, 5)  # Select the relevant column (assuming the 5th column)
cli <- filter(cli, !(cli[, 1] == '' | cli[, 1] == '[Discrepancy]' | cli[, 1] == '[Unknown]' |
                     cli[, 1] == 'Stage X' | cli[, 1] == 'IS' | cli[, 1] == 'I/II NOS' | cli[, 1] == 'Stage 0'))

cli <- cli %>% 
  mutate(tumor_stage = case_when(
    cli[, 1] %in% c('Stage IA', 'Stage IB', 'Stage I') ~ 'Stage I',
    cli[, 1] %in% c('Stage IIA', 'Stage IIB', 'Stage IIC', 'Stage II') ~ 'Stage II',
    cli[, 1] %in% c('Stage IIIA', 'Stage IIIB', 'Stage IIIC', 'Stage III') ~ 'Stage III',
    cli[, 1] %in% c('Stage IVA', 'Stage IVB', 'Stage IVC', 'Stage IV') ~ 'Stage IV'
  ))

# Define gene and clinical variables
gene <- colnames(rt)[1]
clinical <- colnames(cli)[1]

# Loop through each cancer type and perform clinical correlation analysis
for (i in levels(as.factor(rt$CancerType))) {
  rt1 <- rt[rt$CancerType == i, ]
  
  # Combine expression data with gene column
  data <- cbind(rt1, gene = rt1[, gene])
  data <- as.matrix(data[, c(gene, "gene")])
  
  if (nchar(row.names(data)[1]) != nchar(row.names(cli)[1])) {
    row.names(data) <- gsub(".$", "", row.names(data))
  }
  
  data <- avereps(data)
  sameSample <- intersect(row.names(data), row.names(cli))
  sameData <- data[sameSample, ]
  sameClinical <- cli[sameSample, ]
  cliExpData <- cbind(as.data.frame(sameClinical), sameData)
  
  if (nrow(cliExpData) == 0) {
    next
  }
  
  # Set comparison groups
  group <- levels(factor(cliExpData$sameClinical))
  comp <- combn(group, 2)
  my_comparisons <- list()
  for (j in 1:ncol(comp)) {
    my_comparisons[[j]] <- comp[, j]
  }
  
  # Generate box plot
  boxplot <- ggboxplot(cliExpData, x = "sameClinical", y = "gene", color = "sameClinical",
                       xlab = clinical, ylab = paste(gene, "expression"), legend.title = clinical,
                       title = paste0("Cancer: ", i), add = "jitter") + 
    stat_compare_means(comparisons = my_comparisons)
  
  # Save plot to PDF
  pdf(file = paste0(clinical, ".", i, ".pdf"), width = 15, height = 5)
  print(boxplot)
  dev.off()
}
