#!/usr/bin/env Rscript

library(tidyverse)
library(data.table)
library(dplyr)
args = commandArgs(trailingOnly=TRUE)

#susie
susie <- fread(args[1]) %>% rename(PIP=V2) %>% select(PIP)
susie_sumstat <- fread(args[2])
susie_df <- cbind(susie_sumstat, susie) %>% arrange(desc(PIP))
susie_df$cumsum <- cumsum(susie_df$PIP)
susie_df_credset <- susie_df %>% filter(PIP >= 0.95 | cumsum <= 0.95)
fwrite(susie_df_credset,args[3],row.names=F,quote=F,sep="\t")
