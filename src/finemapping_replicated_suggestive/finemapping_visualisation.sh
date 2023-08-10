#!/bin/bash

#Rationale: visualise fine-mapping results with gene in the region, different credible sets.
#Use of LocusZooms function for R as implemented by Tanya J Major and Riku Takei.

PATH_finemapping="/home/n/nnp5/PhD/PhD_project/Fine_mapping_severe_asthma"
PATH_OUT="/home/n/nnp5/PhD/PhD_project/Fine_mapping_severe_asthma/output/"
PATH_ASSOC="/home/n/nnp5/PhD/PhD_project/REGENIE_assoc/output/allchr"
PATH_TMP="/scratch/gen1/nnp5/Fine_mapping/tmp_data/"

module unload R/4.2.1
module load R/4.1.0

##set working directory:
cd ${PATH_finemapping}

#Functional annotation:
cat ${PATH_finemapping}/output/replsugg_shared_valid_credset.txt ${PATH_finemapping}/output/repl_sugg_NOTshared_finemaponly_valid_credset.txt \
    > ${PATH_finemapping}/output/replsugg_valid_credset.txt.tmp
cat ${PATH_finemapping}/output/replsugg_valid_credset.txt.tmp ${PATH_finemapping}/output/repl_sugg_NOTshared_susieonly_valid_credset.txt \
    > ${PATH_finemapping}/output/replsugg_valid_credset.txt
rm ${PATH_finemapping}/output/replsugg_valid_credset.txt.tmp
cp ${PATH_finemapping}/output/replsugg_valid_credset.txt /data/gen1/UKBiobank_500K/severe_asthma/Noemi_PhD/data/

##FAVOR:
##liftover for vars with no rsid:
awk 'NR > 1 {print "chr"$3,$4,$4+1}' ${PATH_finemapping}/output/replsugg_valid_credset.txt \
    > ${PATH_finemapping}/input/replsugg_valid_credset_input_liftover_b37

awk 'NR > 1 {print $5"-"$6}' ${PATH_finemapping}/output/replsugg_valid_credset.txt \
    > /scratch/gen1/nnp5/Fine_mapping/tmp_data/alleles_for_favor_input

awk '{print $1"-"$2}' ${PATH_finemapping}/input/hglft_genome_credset_vars_08_08_2023.bed | sed 's/chr//g' | \
    paste -d "-" - /scratch/gen1/nnp5/Fine_mapping/tmp_data/alleles_for_favor_input \
    > ${PATH_finemapping}/input/replsugg_credset_chrpos38.txt
###use FAVOR webtool:
#https://favor.genohub.org/

##digest FAVOR files:
awk -F "," '{print $1, $2, $3, $8, $9, $10, $11, $12}' ${PATH_finemapping}/input/FAVOR_credset_chrpos38_2023_08_08.csv \
    > ${PATH_finemapping}/input/FAVOR_credset_annotations_digest_08_08_23.csv

##make PIP plot with functional annotation:

for line in {1..18}
do
SNP=$(awk -v row="$line" ' NR == row {print $1 } ' ${PATH_finemapping}/input/fine_mapping_regions_replicated_suggestive_input)
chr=$(awk -v row="$line" ' NR == row {print $2 } ' ${PATH_finemapping}/input/fine_mapping_regions_replicated_suggestive_input)
start=$(awk -v row="$line" 'NR == row {print $4}' ${PATH_finemapping}/input/fine_mapping_regions_replicated_suggestive_input)
end=$(awk -v row="$line" 'NR == row {print $5}' ${PATH_finemapping}/input/fine_mapping_regions_replicated_suggestive_input)

#create R2 according to leading p-value:
module load plink2
plink2 \
    --bgen /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}.bgen ref-first \
    --sample ${PATH_finemapping}/input/ldstore.sample \
    --make-bed --out /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}

#locus="11_rs10160518_75796671_76796671"
#finemap_snp_file="output/finemap_replicated_suggestive/finemap_replsugg_11_rs10160518_75796671_76796671.snp"
snp_lead=$(Rscript src/finemapping_replicated_suggestive/PIP_plot_with_funcannot.R \
    ${chr}_${SNP}_${start}_${end} \
    output/finemap_replicated_suggestive/finemap_replsugg_${chr}_${SNP}_${start}_${end}.snp)

module unload plink2/2.00a
module load plink/1.90
plink --bfile /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP} \
    --extract /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}_no_ma_snps.txt \
    --allow-no-sex --r2 inter-chr --ld-snp ${snp_lead} --ld-window-r2 0 \
    --out /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}_ld_file

Rscript src/finemapping_replicated_suggestive/PIP_plot_with_funcannot_pt2.R \
    /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}_ld_file.ld \
    output/plots/finemapping_plot_${chr}_${SNP}_${start}_${end}.pdf \
    ${chr}_${SNP}_${start}_${end}






####


###experiment using LocusZooms.git repostiory:
#in R:
#cd /home/n/nnp5/software
#git clone https://github.com/Geeketics/LocusZooms.git
# load necessary files into R
library(tidyverse)
library(data.table)
finemapping <- fread("output/finemap_replicated_suggestive/finemap_replsugg_2_rs12470864_102426362_103426362.snp", stringsAsFactor = FALSE, header = TRUE)
finemapping <- finemapping %>% rename(CHR=chromosome, SNP=rsid, BP=position, P=prob)
finemapping <- finemapping %>% select(CHR, SNP, BP, P)
ld <- read.table("/scratch/gen1/nnp5/Fine_mapping/tmp_data/rs12470864_ld_file.ld", stringsAsFactors = FALSE, header = TRUE)
ld <- ld %>% select(SNP_B, R2)
Unique.genes <- read.delim("/home/n/nnp5/software/LocusZooms/Gencode_GRCh37_Genes_UniqueList2021.txt", stringsAsFactors = FALSE, header = TRUE)
secondary_snp <-

# load the locuszoom function into R
source("/home/n/nnp5/software/LocusZooms/functions/locus_zoom.R")

# create a LocusZoom-like plot
locus.zoom(data = finemapping,                                    # a data.frame (or a list of data.frames) with the columns CHR, BP, SNP, and P
           region = c(2, 102426362, 103426362),                             # the chromosome region to be included in the plot
           offset_bp = 0,                                                  # how many basepairs around the SNP / gene / region of interest to plot
           ld.file = ld,                                           # a file with LD values relevant to the SNP specified above
           genes.data = Unique.genes,			                   # a file of all the genes in the region / genome
           plot.title = "Fine-mapping chr2 102426362-103426362",        # the plot title
           file.name = "/output/finemap_replicated_suggestive/plot_finemap_replsugg_2_rs12470864_102426362_103426362.jpg",                                      # the name of the file to save the plot to
           #secondary.snp = c("rs1121980", "rs8060235"),                    # a list of SNPs to label on the plot
           #secondary.label = TRUE
           rsid.check = FALSE
           )