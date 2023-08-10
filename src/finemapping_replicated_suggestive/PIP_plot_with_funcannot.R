#!/usr/bin/env Rscript


#Rationale: plot of fine-mapping PIP with functional annotation for variants in credible set:
#PT1: a.create file with annotation, finemapping variants, and valid credset variants for each locus.
#b.find the leading SNP with max PIP to run plinkv2 for R2 in bash and then create the plot.

library(tidyverse)
library(data.table)
library(dplyr)
args <- commandArgs(T)

#Args:
locus_id <- args[1]
#locus <- credset_annot %>% filter(locus==locus_id)
#the .snp file from finemap
finampping_snps_file <- args[2]
finemapping_snps <- fread(finampping_snps_file)
finemapping_snps <-  finemapping_snps %>% rename(position_b37=position,snpid=rsid)

annot <- read.table("input/FAVOR_credset_annotations_digest_08_08_23.csv",stringsAsFactors = FALSE, fill=TRUE, header=TRUE)
annot <- annot %>% select("Variant..VCF.","Chromosome","Position","Genecode.Comprehensive.Category") %>%
    rename(id="Variant..VCF.",chromosome="Chromosome",position="Position")
annot <- as.data.frame(sapply(annot, function(x) gsub("\"", "", x)))
annot <- annot %>% rename(Functional_annotation=Genecode.Comprehensive.Category)

credset <- fread("output/replsugg_valid_credset.txt",header=TRUE)
b38 <- fread("input/hglft_genome_credset_vars_08_08_2023.bed",header=FALSE)
b38 <- b38 %>% separate(V4, "-", into=c("chr_pos", "pos2"))
b38$chr_pos <- NULL
b38$pos <- as.numeric(b38$pos2) - 1
credset$V1 <- paste0("chr",credset$chromosome)
credset$pos <- credset$position
credset_b38 <- left_join(credset,b38,by=c("V1","pos"))

credset_b38$id <- paste0(credset_b38$chromosome,"-",credset_b38$V2,"-",credset_b38$allele1,"-",credset_b38$allele2)
credset_b38$V1 <- NULL
credset_b38$pos2 <- NULL
credset_b38$pos <- NULL
credset_b38$V3 <- NULL
credset_b38$V5 <- NULL
credset_b38 <- credset_b38 %>% rename(position_b37=position)
credset_b38 <- credset_b38 %>% rename(position=V2)
credset_b38$chromosome <- as.numeric(credset_b38$chromosome)
annot$chromosome <- as.numeric(annot$chromosome)
credset_b38$position <- as.numeric(credset_b38$position)
annot$position <- as.numeric(annot$position)
credset_annot <- inner_join(credset_b38, annot, by=c("id","chromosome","position"))
credset_annot$credset <- "1"


finemapping_snps_credset_annot <- left_join(finemapping_snps,credset_annot,by=c("snpid","chromosome","position_b37","allele1","allele2"))
finemapping_snps_credset_annot <- finemapping_snps_credset_annot %>% select(chromosome,position_b37,prob,Functional_annotation,credset,snpid)
finemapping_snps_credset_annot <- finemapping_snps_credset_annot %>% mutate(Functional_annotation=ifelse(is.na(finemapping_snps_credset_annot$Functional_annotation), "unknown",finemapping_snps_credset_annot$Functional_annotation))
finemapping_snps_credset_annot$Functional_annotation <- as.factor(finemapping_snps_credset_annot$Functional_annotation)

fwrite(finemapping_snps_credset_annot,"/scratch/gen1/nnp5/Fine_mapping/tmp_data/finemapping_snps_credset_annot",quote=F,sep="\t",row.names=F,col.names=T)

##in bash, run plink to create the R2 according to variants with maximum PIP --> extract the leading SNP:
as.character(finemapping_snps_credset_annot[finemapping_snps_credset_annot$prob == max(finemapping_snps_credset_annot$prob),'snpid']) %>% cat()
