#!/bin/bash

#PBS -N SuSie
#PBS -j oe
#PBS -o SuSie_log
#PBS -l walltime=3:0:0
#PBS -l vmem=40gb
#PBS -l nodes=1:ppn=4
#PBS -d .
#PBS -W umask=022

PATH_finemapping="/home/n/nnp5/PhD/PhD_project/Fine_mapping_severe_asthma"
PATH_ASSOC="/home/n/nnp5/PhD/PhD_project/REGENIE_assoc/output/allchr"

cd ${PATH_finemapping}

module load gcc/9.3
module unload R/4.2.1
module load R/4.1.0
module load plink2

#mkdir ${PATH_finemapping}/output/susie

#Input data:
line=13
#for line in {2..14}
#do
SNP=$(awk -v row="$line" ' NR == row {print $1 } ' ${PATH_finemapping}/input/fine_mapping_regions_merged)
chr=$(awk -v row="$line" ' NR == row {print $2 } ' ${PATH_finemapping}/input/fine_mapping_regions_merged)
start=$(awk -v row="$line" 'NR == row {print $4}' ${PATH_finemapping}/input/fine_mapping_regions_merged)
end=$(awk -v row="$line" 'NR == row {print $5}' ${PATH_finemapping}/input/fine_mapping_regions_merged)

#Creating region bgen
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


#Exclude multi-allelic variants and find the common SNP IDs for the genotyped matrix and the zscore input files:
#use the file for each regions created by FINEMAP.sh:
#grep -v -w -F -f /data/gen1/UKBiobank_500K/imputed/multiallelic.snps \
#    ${PATH_finemapping}/input/ldstore_chr${chr}_${SNP}.z \
#    > /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}_no_ma_GWAS_sumstats.txt

#awk 'NR > 1 {print $1}' /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}_no_ma_GWAS_sumstats.txt \
#    > /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}_no_ma_snps.txt

#zscore:
#Rscript src/z_score.R /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}_no_ma_GWAS_sumstats.txt \
#    /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}_no_ma_z_scores.txt


#Format region data for input to R
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

#Rscript src/susie.R \
#    /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}.cols_no_header.raw \
#    /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}_no_ma_z_scores.txt \
#    ${PATH_finemapping}/output/susie/susie_${chr}_${SNP}_${start}_${end}.txt \
#    ${PATH_finemapping}/output/susie/susie_${chr}_${SNP}_${start}_${end}.jpeg

Rscript src/credset_susie.R \
    ${PATH_finemapping}/output/susie/susie_${chr}_${SNP}_${start}_${end}.txt \
    /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}_no_ma_GWAS_sumstats.txt \
    ${PATH_finemapping}/output/susie/susie_credset.${SNP}.$chr.$start.$end

#for 5_rs2188962_rs152815_130026218_132770805, 17:38073838_CCG_C 17 38073838 37073838 39073838:
Rscript src/Susie_trouble_shooting.R \
    /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}.cols_no_header.raw \
    /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}_no_ma_z_scores.txt \
    ${PATH_finemapping}/output/susie/susie_${chr}_${SNP}_${start}_${end}.txt \
    ${PATH_finemapping}/output/susie/susie_${chr}_${SNP}_${start}_${end}.jpeg

awk -F "\t" '{print $6}' ${PATH_finemapping}/output/susie/susie_${chr}_${SNP}_${start}_${end}.txtcredset | \
    tr , '\n' | tail -n +2 > ${PATH_finemapping}/output/susie/susie_${chr}_${SNP}_${start}_${end}.credset.indx

awk -F "\t" '$1 == 1 {print $6}' ${PATH_finemapping}/output/susie/susie_${chr}_${SNP}_${start}_${end}.txtcredset | \
    tr , '\n' | wc -l

awk -F "\t" '$1 == 2 {print $6}' ${PATH_finemapping}/output/susie/susie_${chr}_${SNP}_${start}_${end}.txtcredset | \
    tr , '\n' | wc -l

Rscript src/credset_susie_troubleshooting.R \
    ${PATH_finemapping}/output/susie/susie_${chr}_${SNP}_${start}_${end}.txt \
    ${PATH_finemapping}/output/susie/susie_${chr}_${SNP}_${start}_${end}.credset.indx \
    /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}_no_ma_GWAS_sumstats.txt \
    ${PATH_finemapping}/output/susie/susie_credset.${SNP}.$chr.$start.$end

#done

#Merge all the credible set:
cd ${PATH_finemapping}/output/susie/
#rsid	chromosome	position	allele1	allele2 vars.variable_prob
awk '{print $9, $10, $11, $12, $13, $2}' susie_credset.rs705705.12.55935504.56935504 \
    > ${PATH_finemapping}/output/susie/susie_all_credset.txt && \

awk '{print $9, $10, $11, $12, $13, $2}' susie_credset.rs3024971.12.56993727.57993727 | tail -n +2 \
    >> ${PATH_finemapping}/output/susie/susie_all_credset.txt && \

awk '{print $1, $2, $3, $4, $5, $10}' susie_credset.17:38073838_CCG_C.17.37073838.39073838 | tail -n +2 \
    >> ${PATH_finemapping}/output/susie/susie_all_credset.txt && \

awk '{print $9, $10, $11, $12, $13, $2}' susie_credset.rs12470864.2.101926362.103926362 | tail -n +2 \
    >> ${PATH_finemapping}/output/susie/susie_all_credset.txt && \

awk '{print $9, $10, $11, $12, $13, $2}' susie_credset.rs6761047.2.241692858.243692858 | tail -n +2 \
    >> ${PATH_finemapping}/output/susie/susie_all_credset.txt && \

awk '{print $9, $10, $11, $12, $13, $2}' susie_credset.rs778801698.3.49024027.51024027 | tail -n +2 \
    >> ${PATH_finemapping}/output/susie/susie_all_credset.txt && \

awk '{print $9, $10, $11, $12, $13, $2}' susie_credset.rs1837253.5.109401872.111401872 | tail -n +2 \
    >> ${PATH_finemapping}/output/susie/susie_all_credset.txt && \

awk '{print $9, $10, $11, $12, $13, $2}' susie_credset.rs7824394.8.80292599.82292599 | tail -n +2 \
    >> ${PATH_finemapping}/output/susie/susie_all_credset.txt && \

awk '{print $9, $10, $11, $12, $13, $2}' susie_credset.rs992969.9.5209697.7209697 | tail -n +2 \
    >> ${PATH_finemapping}/output/susie/susie_all_credset.txt && \

awk '{print $9, $10, $11, $12, $13, $2}' susie_credset.rs10160518.11.75296671.77296671 | tail -n +2 \
    >> ${PATH_finemapping}/output/susie/susie_all_credset.txt && \

awk '{print $9, $10, $11, $12, $13, $2}' susie_credset.rs17293632.15.66442596.68442596 | tail -n +2 \
    >> ${PATH_finemapping}/output/susie/susie_all_credset.txt && \

awk '{print $9, $10, $11, $12, $13, $2}' susie_credset.rs201499805_rs1444789.10.8042744.10064361 | tail -n +2 \
    >> ${PATH_finemapping}/output/susie/susie_all_credset.txt && \

awk '{print $1, $2, $3, $4, $5, $10}' susie_credset.rs2188962_rs152815.5.130026218.132770805 | tail -n +2 \
    >> ${PATH_finemapping}/output/susie/susie_all_credset.txt






#if there is some NAs in the genotype matrix (keep it because it is a nice script)
#awk '{ lines[NR] = $0; for (i = 1; i <= NF; i++) if ($i == "NA") skip[i] = 1;} END { for (i = 1; i <= NR; i++) {
#    nf = split(lines[i], fields);
#    for (j = 1; j <= nf; j++) if (!(j in skip)) printf("%s ", fields[j]);
#    printf("\n");
#    }
#    }' /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}.cols_no_header.raw  > /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}.cols_no_header_noNA.raw