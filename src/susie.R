#usage: Rscript susie.R <genotype matrix> <z scores> <plot name>
#remotes::install_github("stephenslab/susieR")

args = commandArgs(trailingOnly=TRUE)

library(tidyverse)
library(data.table)
library(susieR)
library(Rfast)

print("Reading data...")

data <- as.matrix(fread(args[1], header=FALSE))

print("Data read.")

print("Creating correlation matrix...")

data.cor = cor(data)

print("Correlation matrix complete.")

z_scores <- read.table(args[2])

print("Running Susie...")
fitted_rss <- susie_rss(z_scores$V1, data.cor, L=10, n=46086)

df <- as.data.frame(summary(fitted_rss)) %>% select("vars.variable","vars.variable_prob","vars.cs","cs.cs","cs.cs_log10bf","cs.cs_avg_r2","cs.cs_min_r2")
#write.table(df,args[3],row.names=FALSE,quote=FALSE,col.names=T,sep="\t",na="NA")
fread(df,args[3],header=TRUE,quote=FALSE,sep="\t")

jpeg(args[4])

susie_plot(fitted_rss, y="PIP")

dev.off()