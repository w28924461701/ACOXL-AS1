# Set working directory
setwd("D:\\panCancer\\15.TMBcor")

# Read expression and TMB files
exp <- read.table("singleGeneExp.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
TMB <- read.table("TMB.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

# Remove normal samples
group <- sapply(strsplit(row.names(exp), "\\-"), "[", 4)
group <- sapply(strsplit(group, ""), "[", 1)
group <- gsub("2", "1", group)
exp <- exp[group == 0, ]

# Find common samples
sameSample <- intersect(row.names(TMB), row.names(exp))
TMB <- TMB[sameSample, ]
TMB$CancerType <- as.factor(TMB$CancerType)
exp <- exp[sameSample, ]
exp$CancerType <- as.factor(exp$CancerType)

# Initialize output tables
outTab <- data.frame()
fmsbTab <- data.frame()

# Correlation analysis for each cancer type
for (i in levels(exp$CancerType)) {
  exp1 <- exp[exp$CancerType == i, ]
  TMB1 <- TMB[TMB$CancerType == i, ]
  
  # Perform Spearman correlation test
  x <- as.numeric(TMB1[, 1])
  y <- as.numeric(exp1[, 1])
  corT <- cor.test(x, y, method = "spearman")
  cor <- corT$estimate
  pValue <- corT$p.value
  sig <- ifelse(pValue < 0.001, "***", ifelse(pValue < 0.01, "**", ifelse(pValue < 0.05, "*", " ")))
  
  # Append results to output tables
  outTab <- rbind(outTab, cbind(CancerType = i, cor = cor, pValue = pValue, sig))
  fmsbTab <- rbind(fmsbTab, cbind(CancerType = i, cor = cor))
}

# Save results to files
write.table(outTab, file = "corStat.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(t(fmsbTab), file = "fmsbInput.txt", sep = "\t", col.names = FALSE, quote = FALSE)

# Ensure necessary package is installed
if (!require("fmsb")) install.packages("fmsb")

# Load the required library
library(fmsb)

# Set working directory
setwd("D:\\panCancer\\16.TMBradar")

# Read input data
data <- read.table("fmsbInput.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

# Determine maximum value for scaling
maxValue <- ceiling(max(abs(data)) * 10) / 10
data <- rbind(rep(maxValue, ncol(data)), rep(-maxValue, ncol(data)), data)

# Define colors and significance
colors <- "red"
corStat <- read.table("corStat.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
colnames(data) <- paste0(colnames(data), corStat$sig)

# Output radar chart to PDF
pdf(file = "radar.pdf", height = 7, width = 7)
radarchart(data, axistype = 1, 
           pcol = colors,               # Set line color
           plwd = 2,                    # Line width
           plty = 1,                    # Line type (solid)
           cglcol = "grey",             # Background line color
           cglty = 1,                   # Background line type (solid)
           caxislabels = seq(-maxValue, maxValue, maxValue / 2),  # Axis labels
           cglwd = 1.2,                 # Background line width
           axislabcol = "blue",         # Axis label color
           vlcex = 0.8                  # Font size
)
dev.off()

