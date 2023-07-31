#!/usr/bin/env Rscript

#Rationale:

library(tidyverse)
library(data.table)
library(dplyr)
args = commandArgs(trailingOnly=TRUE)
finemap <- fread(args[1])
vars_in_susie <- fread(args[2],fill=TRUE)
idx_pip_susie <- fread(args[3])

finemap <- finemap %>% rename(PIP_finemap=prob)
finemap$credset <- NULL
colnames(vars_in_susie) <- c("idx","snpid","chromosome","position","allele1","allele2")
colnames(idx_pip_susie) <- c("idx","PIP_susie")
vars_pip_susie <- inner_join(vars_in_susie,idx_pip_susie,by="idx")
finemap_susie <- inner_join(finemap,vars_pip_susie,by=c("snpid","chromosome","position","allele1","allele2"))
finemap_susie$idx <- NULL
finemap_susie$PIP_diff <- abs(finemap_susie$PIP_finemap - finemap_susie$PIP_susie)
finemap_susie_valid <- finemap_susie %>% filter(PIP_diff < 0.05)
finemap_susie_valid$PIP_average <- (finemap_susie_valid$PIP_finemap + finemap_susie_valid$PIP_susie) / 2
finemap_susie_valid$PIP_diff <- NULL
write.table(finemap_susie_valid,"output/repl_sugg_NOTshared_finemaponly_valid_credset.txt",row.names=FALSE,quote=FALSE,col.names=F,sep="\t",append=T)