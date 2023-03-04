#usage: Rscript susie.R <genotype matrix> <z scores> <plot name>
#remotes::install_github("stephenslab/susieR")

args = commandArgs(trailingOnly=TRUE)

library(data.table)
library(susieR)

print("Reading data...")

data <- as.matrix(fread(args[1], header=FALSE))

print("Data read.")

print("Creating correlation matrix...")

data.cor = cor(data)

print("Correlation matrix complete.")

z_scores <- read.table(args[2])

print("Running Susie...")
fitted_rss <- susie_rss(z_scores$V1, data.cor, L=10, n=46086)

write.table(summary(fitted_rss),args[3],row.names=FALSE,quote=FALSE,col.names=T,sep="\t")

jpeg(args[4])

susie_plot(fitted_rss, y="PIP")

dev.off()