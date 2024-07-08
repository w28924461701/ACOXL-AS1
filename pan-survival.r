# Ensure necessary packages are installed
if (!require("survival")) install.packages("survival")
if (!require("survminer")) install.packages("survminer")
if (!require("dplyr")) install.packages("dplyr")

# Load required libraries
library(survival)
library(survminer)
library(dplyr)

# Set working directory
setwd("D:\\panCancer\\08.survival")

# Read the data file
rt <- read.table("expTime.txt", header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)

# Convert follow-up time from days to years
rt$futime <- rt$futime / 365

# Define the gene of interest
gene <- colnames(rt)[3]
pFilter <- 0.05  # Threshold for p-value significance

# Loop through each cancer type
for (i in levels(as.factor(rt$CancerType))) {
  # Subset data for the specific cancer type
  rt1 <- rt[rt$CancerType == i, ]
  
  # Determine cutpoint for high and low expression groups
  res.cut <- surv_cutpoint(rt1, time = "futime", event = "fustat", variables = c(gene))
  group <- ifelse(rt1[, gene] > res.cut$cutpoint, "high", "low")
  
  # Perform survival difference analysis
  diff <- survdiff(Surv(futime, fustat) ~ group, data = rt1)
  pValue <- 1 - pchisq(diff$chisq, df = 1)
  
  # Check if the p-value is below the defined threshold
  if (pValue < pFilter) {
    if (pValue < 0.001) {
      pValue <- "p<0.001"
    } else {
      pValue <- paste0("p=", sprintf("%.03f", pValue))
    }
    
    # Fit survival curves
    fit <- survfit(Surv(futime, fustat) ~ group, data = rt1)
    
    # Generate the survival plot
    surPlot <- ggsurvplot(fit, 
                          data = rt1,
                          title = paste0("Cancer: ", i),
                          pval = pValue,
                          pval.size = 6,
                          legend.labs = c("high", "low"),
                          legend.title = paste0(gene, " levels"),
                          font.legend = 12,
                          xlab = "Time (years)",
                          ylab = "Overall survival",
                          break.time.by = 1,
                          palette = c("red", "blue"),
                          conf.int = FALSE,
                          fontsize = 4,
                          risk.table = TRUE,
                          risk.table.title = "",
                          risk.table.height = .25)
    
    # Save the plot to a PDF file
    pdf(file = paste0("survival.", i, ".pdf"), onefile = FALSE, width = 6, height = 5)
    print(surPlot)
    dev.off()
  }
}
