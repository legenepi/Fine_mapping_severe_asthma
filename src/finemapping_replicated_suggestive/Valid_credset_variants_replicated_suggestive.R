#!/usr/bin/env Rscript

#Rationale: identify fine-mapped variants to follow-up in variant-to-gene analysis.
##How:
#Variants in 95% credible set with evidence supported by both FINEMAP and SuSiE AND
#PIP difference smaller than 0.05
#(as per Masai et al. : if PIP difference greater, more chance that variants are false positive and show poor enrichment for known functional annotations)

library(tidyverse)
library(data.table)
library(dplyr)
args = commandArgs(trailingOnly=TRUE)

#run as:
#Rscript Fine_mapping_severe_asthma/src/finemapping_replicated_suggestive/Valid_credset_variants_replicated_suggestive.R \
#Fine_mapping_severe_asthma/output/susie_replsugg_all_credset.txt \
#Fine_mapping_severe_asthma/output/finemap_replsugg_all_credset.txt

susie <- fread(arg[1]) %>% select(rsid,chromosome,position,allele1,allele2,vars.variable_prob)
susie <- susie %>% rename(snpid=rsid, PIP_susie=vars.variable_prob)
susie$PIP_susie <- as.numeric(susie$PIP_susie)
finemap <- fread(arg[2])
finemap$locus <- NULL
finemap <- finemap %>% rename(PIP_finemap=prob)
finemap$PIP_finemap <- as.numeric(finemap$PIP_finemap)

#inner join on cols: snpid,chromosome,position,allele1,allele2
shared_vars <- inner_join(finemap, susie, by=c(snpid,chromosome,position,allele1,allele2))
shared_vars$PIP_diff <- abs(shared_vars$PIP_susie - shared_vars$PIP_finemap)

#valid shared variants according to PIP concordance:
shared_valid <- shared_vars %>% filter(PIP_diff < 0.05)

