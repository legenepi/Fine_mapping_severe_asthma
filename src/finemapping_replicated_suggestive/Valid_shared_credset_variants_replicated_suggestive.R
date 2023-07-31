#!/usr/bin/env Rscript

#Rationale: identify fine-mapped variants to follow-up in variant-to-gene analysis.
##How:
#Variants in 95% credible set with evidence supported by both FINEMAP and SuSiE AND
#PIP difference smaller than 0.05
#(as per Masai et al. : if PIP difference greater, more chance that variants are false positive and show poor enrichment for known functional annotations)

##Run as:
#Rscript src/finemapping_replicated_suggestive/Valid_shared_credset_variants_replicated_suggestive.R \
#output/susie_replsugg_all_credset.txt \
#output/finemap_replsugg_all_credset.txt \
#output/replsugg_shared_valid_credset.txt

library(tidyverse)
library(data.table)
library(dplyr)
args = commandArgs(trailingOnly=TRUE)

susie <- fread(args[1])
susie <- susie %>% rename(snpid=rsid, PIP_susie=PIP)
susie$PIP_susie <- as.numeric(susie$PIP_susie)
print(susie)
finemap <- fread(args[2])
finemap <- finemap %>% rename(PIP_finemap=prob)
finemap$PIP_finemap <- as.numeric(finemap$PIP_finemap)
print(finemap)

#inner join on cols: snpid,chromosome,position,allele1,allele2
shared_vars <- inner_join(finemap, susie, by=c("locus","snpid","chromosome","position","allele1","allele2"))
count_shared <- shared_vars %>% count(locus, name= "N_SNP_by_locus")
shared_vars$credset <- NULL
print("Number of shared credible set SNPs for each genomic loci:")
print(count_shared)

#valid shared variants according to PIP concordance:
shared_vars$PIP_diff <- abs(shared_vars$PIP_susie - shared_vars$PIP_finemap)
write.table(shared_vars,"output/repl_sugg_shared_credset.txt",row.names=FALSE,quote=FALSE,col.names=TRUE,sep="\t")
shared_valid <- shared_vars %>% filter(PIP_diff < 0.0555)
shared_valid$PIP_diff <- NULL
shared_valid$PIP_average <- (shared_valid$PIP_finemap + shared_valid$PIP_susie) / 2
count_valid <- shared_valid %>% count(locus, name= "N_SNP_by_locus")
print("Number of valid credible set SNPs for each genomic loci:")
print(count_valid)
write.table(shared_valid,args[3],row.names=FALSE,quote=FALSE,col.names=TRUE,sep="\t")

