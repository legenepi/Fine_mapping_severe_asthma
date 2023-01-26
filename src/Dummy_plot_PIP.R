#!/usr/bin/env Rscript

suppressMessages(library(tidyverse))
suppressMessages(library(data.table))

args = commandArgs(trailingOnly=TRUE)

##Run as:
#Rscript src/Dummy_plot_PIP.R output/manual_finemap.UKB.2.242192858.243192858.gz finemap.UKB.2.242192858.243192858

input_file = args[1]
chr = args[2]

#input:
finemap <- fread(input_file)
min_pvalue <- as.numeric(min(finemap$P))
finemap %>% ggplot(aes(x=BP, y=PIP)) + geom_point(aes(fill = P))
ggsave(paste0("/lustre/ahome3/n/nnp5/PhD/PhD_project/Fine_mapping_severe_asthma/output/plot_",chr,".png"))