#!/usr/bin/env Rscript

library(tidyverse)
library(data.table)

susie_in_finemap <- fread("output/susie_only_credset_vars_in_finemap")
colnames(susie_in_finemap) <- c("index","rsid","chromosome","position","allele1","allele2","maf","beta","se","z","prob","log10bf","mean", "sd","mean_incl", "sd_incl")
susie_in_finemap <- susie_in_finemap %>% select("index","rsid","chromosome","position","allele1","allele2","prob")
susie_in_finemap <- susie_in_finemap %>% rename(finemap_PIP=prob)
susie <- fread("output/susie/susie_all_credset.txt") %>% rename(susie_PIP=vars.variable_prob)

merged <- left_join(susie_in_finemap,susie,by=c("rsid","chromosome","position","allele1","allele2"))
merged$diff_PIP <- abs(merged$susie_PIP - merged$finemap_PIP)
merged_valid <- merged %>% filter(diff_PIP < 0.05)
merged_valid$avg_PIP <- (merged_valid$susie_PIP + merged_valid$finemap_PIP) / 2
merged_valid$Status <- "susie"
merged_valid <- merged_valid %>% select(rsid,chromosome,position,allele1,allele2,finemap_PIP,susie_PIP,avg_PIP,Status)
fwrite(merged_valid,"output/credset_NOTshared_susieonly_valid_SNPs",row.names=F,quote=F,sep="\t")


finemap_in_susie <- fread("output/finemap_only_credset_vars_in_susie")
colnames(finemap_in_susie) <- c('rsid', 'chromosome', 'position', 'allele1', 'allele2', 'vars.variable_prob')
finemap_in_susie <- finemap_in_susie %>% rename(susie_PIP=vars.variable_prob)
finemap <- fread("output/finemap/finemap_all_credset.txt") %>% rename(finemap_PIP=prob, rsid=snpid)
merged <- left_join(finemap_in_susie,finemap, by=c("rsid","chromosome","position","allele1","allele2"))
merged$diff_PIP <- abs(merged$susie_PIP - merged$finemap_PIP)
merged_valid <- merged %>% filter(diff_PIP < 0.05)
merged_valid$avg_PIP <- (merged_valid$susie_PIP + merged_valid$finemap_PIP) / 2
merged_valid$Status <- "finemap"
merged_valid <- merged_valid %>% select(rsid,chromosome,position,allele1,allele2,finemap_PIP,susie_PIP,avg_PIP,Status)
fwrite(merged_valid,"output/credset_NOTshared_finemaponly_valid_SNPs",row.names=F,quote=F,sep="\t")
