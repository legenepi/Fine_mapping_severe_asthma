#!/usr/bin/env Rscript

suppressMessages(library(tidyverse))
suppressMessages(library(data.table))

args = commandArgs(trailingOnly=TRUE)

locus = args[1]
gwas <- fread("/home/n/nnp5/PhD/PhD_project/REGENIE_assoc/output/maf001_broad_pheno_1_5_ratio_betase_input_mungestat")
fm.path = "/home/n/nnp5/PhD/PhD_project/Fine_mapping_severe_asthma/output/finemap_replicated_suggestive/finemap_replsugg_"

fm.snp = read.table(paste0(fm.path,locus,".snp"), as.is = T, header = T)
colnames(fm.snp)[2] <- "snpid"
gwas_fm <- left_join(fm.snp,gwas,by="snpid")
gwas_fm$credset <- cumsum(gwas_fm$prob)
gwas_fm_credset <- gwas_fm %>% filter(credset <= 0.95559 | prob >= 0.95559)
gwas_fm_credset$locus <- as.character(locus)
gwas_fm_credset_digest <- gwas_fm_credset %>% select(c(locus,snpid,chromosome,position,allele1,allele2,prob,credset))
fwrite(gwas_fm_credset_digest,paste0(fm.path,"credset_",locus,".txt"),quote=FALSE,row.names=FALSE,col.names=T,na="NA")