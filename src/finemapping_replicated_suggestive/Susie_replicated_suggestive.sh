#!/bin/bash

#PBS -N SuSie_replicated_suggestive
#PBS -j oe
#PBS -o SuSie_replicated_suggestive
#PBS -l walltime=10:0:0
#PBS -l vmem=50gb
#PBS -l nodes=1:ppn=4
#PBS -d .
#PBS -W umask=022

PATH_finemapping="/home/n/nnp5/PhD/PhD_project/Fine_mapping_severe_asthma"
PATH_OUT="/home/n/nnp5/PhD/PhD_project/Fine_mapping_severe_asthma/output/susie_replicated_suggestive"
PATH_ASSOC="/home/n/nnp5/PhD/PhD_project/REGENIE_assoc/output/allchr"

##set working directory:
cd ${PATH_finemapping}

##load required tools:
module load gcc/9.3
module unload R/4.2.1
module load R/4.1.0
module load plink2

##output folder:
#mkdir ${PATH_finemapping}/output/susie_replicated_suggestive


##Analysis:
for line in {1..21}
do
##Input data:
SNP=$(awk -v row="$line" ' NR == row {print $1 } ' ${PATH_finemapping}/input/fine_mapping_regions_replicated_suggestive)
chr=$(awk -v row="$line" ' NR == row {print $2 } ' ${PATH_finemapping}/input/fine_mapping_regions_replicated_suggestive)
start=$(awk -v row="$line" 'NR == row {print $4}' ${PATH_finemapping}/input/fine_mapping_regions_replicated_suggestive)
end=$(awk -v row="$line" 'NR == row {print $5}' ${PATH_finemapping}/input/fine_mapping_regions_replicated_suggestive)

##Creating region bgen
#if chr has double digit, I do not need the '0' in the chromosome name, so I need two different string for the
#-incl-range argument:
#cd /scratch/gen1/nnp5/Fine_mapping/tmp_data/
#if [[ ${chr} -lt 10 ]]
#then
#~nrgs1/bin/bgenix -g /scratch/gen1/nnp5/Fine_mapping/tmp_data/sevasthma_chr${chr}_v3.bgen \
#    -incl-range 0${chr}:${start}-${end} \
#    > /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}.bgen
#fi

#if [[ ${chr} -gt 9 ]]
#then
#  ~nrgs1/bin/bgenix -g /scratch/gen1/nnp5/Fine_mapping/tmp_data/sevasthma_chr${chr}_v3.bgen \
#    -incl-range ${chr}:${start}-${end} \
#    > /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}.bgen
#fi

#cd ${PATH_finemapping}


##Exclude multi-allelic variants and find the common SNP IDs for the genotyped matrix and the zscore input files:
#use the file for each regions created by FINEMAP.sh:
#grep -v -w -F -f /data/gen1/UKBiobank_500K/imputed/multiallelic.snps \
#    ${PATH_finemapping}/input/ldstore_chr${chr}_${SNP}.z \
#    > /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}_no_ma_GWAS_sumstats.txt

#awk 'NR > 1 {print $1}' /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}_no_ma_GWAS_sumstats.txt \
#    > /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}_no_ma_snps.txt

##zscore:
#Rscript src/z_score.R /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}_no_ma_GWAS_sumstats.txt \
#    /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}_no_ma_z_scores.txt


##Format region data for input to R
#plink2 \
#    --bgen /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}.bgen ref-first \
#    --sample ${PATH_finemapping}/input/ldstore.sample \
#    --export A \
#    --extract /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}_no_ma_snps.txt \
#    --out /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}

#cut -f7- /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}.raw \
#    > /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}.cols.raw

#awk 'NR>1 {print}' /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}.cols.raw \
#    > /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}.cols_no_header.raw

Rscript src/finemapping_replicated_suggestive/Susie_replicated_suggestive.R \
    /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}.cols_no_header.raw \
    /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}_no_ma_z_scores.txt \
    ${PATH_OUT}/susie_replsugg_${chr}_${SNP}_${start}_${end}.txt \
    ${PATH_OUT}/susie_replsugg_${chr}_${SNP}_${start}_${end}.jpeg

awk -F "\t" '{print $6}' ${PATH_OUT}/susie_replsugg_${chr}_${SNP}_${start}_${end}.txtcredset | \
    tr , '\n' | tail -n +2 > ${PATH_OUT}/susie_replsugg_${chr}_${SNP}_${start}_${end}.credset.indx

Rscript src/finemapping_replicated_suggestive/credset_susie.R \
    ${PATH_OUT}/susie_replsugg_${chr}_${SNP}_${start}_${end}.txt\
    ${PATH_OUT}/susie_replsugg_${chr}_${SNP}_${start}_${end}.credset.indx \
    /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}_no_ma_GWAS_sumstats.txt \
    ${PATH_OUT}/susie_replsugg_credset.${SNP}.$chr.$start.$end

done