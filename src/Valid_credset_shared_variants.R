#!/usr/bin/env Rscript

library(tidyverse)
library(data.table)

finemap <- fread("output/finemap/finemap_all_credset.txt") %>% rename(finemap_PIP=prob)
susie <- fread("output/susie/susie_all_credset.txt") %>% rename(snpid=rsid, susie_PIP=vars.variable_prob)

shared <- inner_join(finemap,susie,by=c("snpid","chromosome","position","allele1","allele2"))
shared$diff_PIP <- abs(shared$susie_PIP - shared$finemap_PIP)
shared_valid <- shared %>% filter(diff_PIP < 0.05)
shared_valid$avg_PIP <- (shared_valid$susie_PIP + shared_valid$finemap_PIP) / 2
shared_valid$Status <- "shared"
shared_valid <- shared_valid %>% select(snpid,chromosome,position,allele1,allele2,finemap_PIP,susie_PIP,avg_PIP,Status)
fwrite(shared_valid,"output/credset_shared_valid_SNPs",row.names=F,quote=F,sep="\t")

finemap_only <- anti_join(finemap,susie,by=c("snpid","chromosome","position","allele1","allele2"))
finemap_only$index <- NULL
finemap_only$credset <- NULL
fwrite(finemap_only,"output/finemap_only_credset_vars",row.names=F,quote=F,sep="\t")
susie_only <- anti_join(susie,finemap,by=c("snpid","chromosome","position","allele1","allele2"))
fwrite(susie_only,"output/susie_only_credset_vars",row.names=F,quote=F,sep="\t")

