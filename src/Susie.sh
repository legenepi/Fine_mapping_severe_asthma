#!/bin/bash

#PBS -N SuSie
#PBS -j oe
#PBS -o SuSie
#PBS -l walltime=1:0:0
#PBS -l vmem=50gb
#PBS -l nodes=1:ppn=1
#PBS -d .
#PBS -W umask=022

PATH_finemapping="/home/n/nnp5/PhD/PhD_project/Fine_mapping_severe_asthma"
PATH_ASSOC="/home/n/nnp5/PhD/PhD_project/REGENIE_assoc/output/allchr"

module load gcc/9.3
module unload R/4.2.1
module load R/4.1.0
module load plink2

#Input data:
chr_row=1
for line in {2..18}
do
SNP=$(awk -v row="$line" -F "\t" ' NR == row {print $1 } ' /home/n/nnp5/PhD/PhD_project/REGENIE_assoc/output/maf001_broad_pheno_1_5_ratio_sentinel_variants.txt)
chr=$(awk -v row="$line" -F "\t" ' NR == row {print $2 } ' /home/n/nnp5/PhD/PhD_project/REGENIE_assoc/output/maf001_broad_pheno_1_5_ratio_sentinel_variants.txt)

start=$(awk -v row=$chr_row 'NR == row {print $2}' /scratch/gen1/nnp5/Fine_mapping/tmp_data/fine_mapping_regions)
end=$(awk -v row=$chr_row 'NR == row {print $3}' /scratch/gen1/nnp5/Fine_mapping/tmp_data/fine_mapping_regions)
chr_row=$((chr_row + 2))

#I already have plink files for each regions:
#to create plink files with my GWAS participants:
#qsub /home/n/nnp5/PhD/PhD_project/Post_GWAS/plink_conversion_files.sh
#to extract plink file or regions +/- 1000Kb from the sentinel variants.
module load plink
plink \
    --bfile /scratch/gen1/nnp5/REGENIE_assoc/tmp_data/${pheno}_plink_file_v3_chr${chr} \
    --snp ${SNP} \
    --window 1000 \
    --make-bed --out /scratch/gen1/nnp5/REGENIE_assoc/tmp_data/${pheno}_plink_file_v3_chr${chr}_${SNP}

#Retrieve z-scores for variants in the region:
#to calculate zscores:
#applied filters #filters: --min-info 0.85 --min-maf 0.01
#python ${PATH_finemapping}/src/Input_noprior_parquet.py \
#    ${PATH_finemapping}/input/maf001_broad_pheno_1_5_ratio_betase_input_mungestat



#Creating region bgen
~nrgs1/bin/bgenix -g /scratch/gen1/nnp5/Fine_mapping/tmp_data/severeasthma_EUR_ukb_imp_chr3_v3.bgen \
    -incl-range ${chr}:${start}-${end} \
    > /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}.bgen

#Exclude multi-allelic variants and find the common SNP IDs for the genotyped matrix and the zscore input files:
#use the file for each regions created by FINEMAP.sh:
grep -v -w -F -f /data/gen1/UKBiobank_500K/imputed/multiallelic.snps \
    ${PATH_finemapping}/input/ldstore_chr${chr}_${SNP}.z \
    > /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}_no_ma_GWAS_sumstats.txt

awk 'NR > 1 {print $1}' /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}_no_ma_GWAS_sumstats.txt \
    > /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}_no_ma_snps.txt

#zscore:
Rscript src/z_score.R /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}_no_ma_GWAS_sumstats.txt \
    /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}_no_ma_z_scores.txt


#Format region data for input to R
plink2 \
    --bgen /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}.bgen ref-first \
    --sample ${PATH_finemapping}/input/ldstore.sample \
    --export A \
    --extract /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}_no_ma_snps.txt \
    --out /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}

cut -f7- /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}.raw \
    > /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}.cols.raw

awk 'NR>1 {print}' /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}.cols.raw \
    > /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}.cols_no_header.raw

Rscript src/susie.R \
    /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}.cols_no_header.raw \
    /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}_no_ma_z_scores.txt \
    ${PATH_finemapping}/output/susie_${chr}_${SNP}.txt \
    ${PATH_finemapping}/output/susie_${chr}_${SNP}.jpeg

done



#post analysis: Found the credible set variants as per susie:
awk 'NR <= 2 {print $8}' ${PATH_finemapping}/output/susie_3_rs778801698.txt | tr , '\n' | tail -n +2 \
    > ${PATH_finemapping}/output/susie_3_rs778801698_cs_index.txt

awk '{print $1, $2, $3, $4, $5, $6, $7}' ${PATH_finemapping}/output/susie_3_rs778801698.txt | \
    grep -w -F -f ${PATH_finemapping}/output/susie_3_rs778801698_cs_index.txt - \
    > ${PATH_finemapping}/output/susie_3_rs778801698.txt.digest


awk 'NR <= 2 {print $8}' ${PATH_finemapping}/output/susie_3_rs778801698.txt | tr , '\n' | tail -n +2 | awk '{print $1+1}' \
    > ${PATH_finemapping}/output/susie_3_rs778801698_cs_index2.txt

awk '{print NR,$1}' /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}_no_ma_GWAS_sumstats.txt | \
    grep -w -F -f ${PATH_finemapping}/output/susie_3_rs778801698_cs_index2.txt - | awk '{print $2}' | \
    grep -w -F -f - /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}_no_ma_GWAS_sumstats.txt | awk '{print $1}' \
    > ${PATH_finemapping}/output/susie_3_rs778801698_cs_snps.txt


















#if there is some NAs in the genotype matrix (but I shouldnt need it now)
awk '{ lines[NR] = $0; for (i = 1; i <= NF; i++) if ($i == "NA") skip[i] = 1;} END { for (i = 1; i <= NR; i++) {
    nf = split(lines[i], fields);
    for (j = 1; j <= nf; j++) if (!(j in skip)) printf("%s ", fields[j]);
    printf("\n");
    }
    }' /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}.cols_no_header.raw  > /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}.cols_no_header_noNA.raw