#remotes::install_github("stephenslab/susieR")

args = commandArgs(trailingOnly=TRUE)

library(tidyverse)
library(data.table)
library(susieR)

print("Reading data...")

data <- as.matrix(fread(args[1], header=FALSE))

print("Data read.")

print("Creating correlation matrix...")

data.cor = cor(data)

print("Correlation matrix complete.")

sumstat <- read.table(args[2])

print("Running Susie...")

#coverage default is 0.95:
fitted_rss <- susie_rss(z = sumstat$V1, R = data.cor, L = 10, n=46086)

summary(fitted_rss)$cs

print(fitted_rss$pip)

write.table(fitted_rss$pip,args[3],row.names=TRUE,quote=FALSE,col.names=FALSE,sep="\t")

jpeg(args[4])

susie_plot(fitted_rss, y="PIP")

dev.off()