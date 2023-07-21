#!/usr/bin/env Rscript

suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(readxl))
suppressMessages(library(VennDiagram))
suppressMessages(library(gplots))

args = commandArgs(trailingOnly=TRUE)

gwas_file <- args[1]
finemap_file <- args[2]
susie_file <- args[3]
susie_sumstat_file <- args[4]
polyfun_finemap_file <- args[5]
polyfun_susie_file <- args[6]
chr <- as.numeric(args[7])
start <- as.numeric(args[8])
end <- as.numeric(args[9])
venn_name <- paste0("output/",chr,"_",start,"_",end,"_venn_comparison.png")
finemapping_plot_file <- paste0("output/",chr,"_",start,"_",end,"_finemapping_plots.png")
print(gwas_file)
print(finemap_file)
print(susie_file)
print(susie_sumstat_file)
print(polyfun_finemap_file)
print(polyfun_susie_file)
print(chr)
print(start)
print(end)

#gwas
gwas <- fread(gwas_file)
gwas <- gwas %>% filter(b37chr == chr, bp >= start, bp <= end)

#finemap
finemap <- read.table(finemap_file, as.is = T, header = T)
#find indexes for causal variants in 95% credible set:
finemap$cumsum_pip <- cumsum(finemap$prob)
finemap_credset <- finemap %>% filter(cumsum_pip <= 0.95 | prob >= 0.95)

#susie
susie_df <- fread(susie_file)

#polyfun+susie plot:
polyfun_susie <- fread(polyfun_susie_file)
polyfun_susie$cumsum_pip <- cumsum(polyfun_susie$PIP)
ps_credset <- polyfun_susie %>% filter(cumsum_pip <= 0.95 | PIP >= 0.95)

#polyfun+finemap plot:
polyfun_finemap <- fread(polyfun_finemap_file)
polyfun_finemap$cumsum_pip <- cumsum(polyfun_finemap$PIP)
pf_credset <- polyfun_finemap %>% filter(cumsum_pip <= 0.95 | PIP >= 0.95)

png(finemapping_plot_file, width=1300, height=700)
par(mfrow = c(2,3))
par(mar = c(4.5,4.5,3,1))

#gwas
plot(gwas$bp, -log10(gwas$pval), xlab = "SNP (BP)", ylab = "-log10(P)", main = "GWAS")

#finemap
c.bp.fm <- finemap_credset$position
i = which(finemap$position %in% c.bp.fm)
plot(finemap$position, finemap$prob, xlab = "SNP", ylab = "prob of causality",main = "FINEMAP", ylim=c(0, 1))
points(finemap[i,"position"], finemap[i,"prob"], cex = 2, lwd = 1.4, col = "dodgerblue")

#susie
c.bp.s <- susie_df$position
i = which(susie_df$position %in% c.bp.s)
plot(susie_df$position, susie_df$vars.variable_prob, xlab = "SNP", ylab = "prob of causality",main = "SUSIE", ylim=c(0, 1))
points(susie_df[i,position], susie_df[i,vars.variable_prob], cex = 2, lwd = 1.4, col = "dodgerblue")

#polyfun susie
c.bp.ps <- ps_credset$BP
i = which(polyfun_susie$BP %in% c.bp.ps)
plot(polyfun_susie$BP, polyfun_susie$PIP, xlab = "SNP", ylab = "prob of causality",main = "POLYFUN+SUSIE", ylim=c(0, 1))
points(polyfun_susie[i,BP], polyfun_susie[i,PIP], cex = 2, lwd = 1.4, col = "dodgerblue")

#polyfun finemap
c.bp.pf <- pf_credset$BP
i = which(polyfun_finemap$BP %in% c.bp.pf)
plot(polyfun_finemap$BP, polyfun_finemap$PIP, xlab = "SNP", ylab = "prob of causality",main = "POLYFUN+FINEMAP", ylim=c(0, 1))
points(polyfun_finemap[i,BP], polyfun_finemap[i,PIP], cex = 2, lwd = 1.4, col = "dodgerblue")
dev.off()


#venn plot:
print(finemap_credset %>% select(position) %>% distinct() %>% unlist())
print(susie_df %>%  select(position) %>% distinct() %>% unlist())
print(ps_credset %>%  select(BP) %>% distinct() %>% unlist())
print(pf_credset %>%  select(BP) %>% distinct() %>% unlist())

venn.diagram(
   x = list(
     finemap_credset %>% select(position) %>% distinct() %>% unlist(),
     susie_df %>%  select(position) %>% distinct() %>% unlist(),
     ps_credset %>%  select(BP) %>% distinct() %>% unlist(),
     pf_credset %>%  select(BP) %>% distinct() %>% unlist()
    ),
   category.names = c("finemap", "susie", "polyfun_susie", "polyfun_finemap"),
   disable.logging = TRUE,
   filename = venn_name,
   output = TRUE ,
           imagetype="png" ,
           height = 1200 ,
           width = 1500 ,
           resolution = 400,
           compression = "lzw",
           lwd = 1,
           col=c("#440154ff", '#21908dff', '#fde725ff', 'Dark Olive Green 4'),
           fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('#fde725ff',0.3), alpha('Dark Olive Green 4',0.3)),
           cex = 0.5,
           fontfamily = "sans",
           cat.cex = 0.3,
           cat.default.pos = "outer")