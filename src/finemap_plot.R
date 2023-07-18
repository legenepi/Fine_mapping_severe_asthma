#!/usr/bin/env Rscript

suppressMessages(library(tidyverse))
suppressMessages(library(data.table))

args = commandArgs(trailingOnly=TRUE)

locus = args[1]
gwas <- fread("/home/n/nnp5/PhD/PhD_project/REGENIE_assoc/output/maf001_broad_pheno_1_5_ratio_betase_input_mungestat")
fm.path = "/home/n/nnp5/PhD/PhD_project/Fine_mapping_severe_asthma/output/finemap/finemap_"

fm.snp = read.table(paste0(fm.path,locus,".snp"), as.is = T, header = T)
colnames(fm.snp)[2] <- "snpid"
gwas_fm <- left_join(fm.snp,gwas,by="snpid")
#find indexes for causal variants in 95% credible set:
gwas_fm$credset <- cumsum(gwas_fm$prob)
gwas_fm_credset <- gwas_fm %>% filter(credset <= 0.95 | prob >= 0.95)
gwas_fm_credset_digest <- gwas_fm_credset %>% select(c(index,snpid,chromosome,position,allele1,allele2,prob,credset))
fwrite(gwas_fm_credset_digest,paste0(fm.path,"credset_",locus,".txt"),quote=FALSE,row.names=FALSE,col.names=T,na="NA")
c.ind.fm <- gwas_fm_credset$index

png(paste0(fm.path,"plot_",locus,".png"), width=1300, height=700)
par(mfrow = c(1,2))
par(mar = c(4.5,4.5,3,1))
i = which(gwas_fm$index %in% c.ind.fm)
plot(gwas_fm$index, -log10(gwas_fm$pval), xlab = "SNP", ylab = "-log10 P", main = "GWAS SEVERE ASTHMA")
points(gwas_fm[i,"index"], -log10(gwas_fm[i,"pval"]), cex = 2, lwd = 1.4, col = "dodgerblue")
plot(gwas_fm$index, gwas_fm$prob, xlab = "SNP", ylab = "prob of causality",main = "FINEMAP SEVERE ASTHMA")
points(fm.snp[i,"index"], fm.snp[i,"prob"], cex = 2, lwd = 1.4, col = "dodgerblue")
dev.off()