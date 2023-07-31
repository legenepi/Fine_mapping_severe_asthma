#!/usr/bin/env Rscript

#Rationale:
#Susie only credset variants to validate with finemap PIP:
#PIP_diff < 0.05 then valid.

library(tidyverse)
library(data.table)
library(dplyr)
args = commandArgs(trailingOnly=TRUE)
susie <- fread(args[1])
finemap <- fread(args[2],fill=TRUE)

susie <- susie %>% rename(PIP_susie=PIP)
finemap <- finemap %>% select("rsid","chromosome","position","allele1","allele2","prob")
finemap <- finemap %>% rename(PIP_finemap=prob)

susie_finemap <- inner_join(susie,finemap,by=c("rsid","chromosome","position","allele1","allele2"))
susie_finemap$PIP_diff <- abs(susie_finemap$PIP_finemap - susie_finemap$PIP_susie)
susie_finemap_valid <- susie_finemap %>% filter(PIP_diff < 0.05)
susie_finemap_valid$PIP_average <- (susie_finemap_valid$PIP_finemap + susie_finemap_valid$PIP_susie) / 2
susie_finemap_valid$PIP_diff <- NULL
write.table(susie_finemap_valid,"output/repl_sugg_NOTshared_susieonly_valid_credset.txt",row.names=FALSE,quote=FALSE,col.names=F,sep="\t",append=T)