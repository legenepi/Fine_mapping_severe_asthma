#!/bin/bash

PATH_finemapping="/home/n/nnp5/PhD/PhD_project/Fine_mapping_severe_asthma"

cd ${PATH_finemapping}

#Visualisation of the results:
chmod o+x src/comparison_credset.R
dos2unix src/comparison_credset.R
for line in {2..14}
do
SNP=$(awk -v row="$line" ' NR == row {print $1 } ' ${PATH_finemapping}/input/fine_mapping_regions_merged)
chr=$(awk -v row="$line" ' NR == row {print $2 } ' ${PATH_finemapping}/input/fine_mapping_regions_merged)
start=$(awk -v row="$line" 'NR == row {print $4}' ${PATH_finemapping}/input/fine_mapping_regions_merged)
end=$(awk -v row="$line" 'NR == row {print $5}' ${PATH_finemapping}/input/fine_mapping_regions_merged)
Rscript ${PATH_finemapping}/src/comparison_credset.R \
    /home/n/nnp5/PhD/PhD_project/REGENIE_assoc/output/maf001_broad_pheno_1_5_ratio_betase_input_mungestat \
     ${PATH_finemapping}/output/finemap/finemap_${chr}_${SNP}_${start}_${end}.snp \
     ${PATH_finemapping}/output/susie/susie_credset.${SNP}.$chr.$start.$end \
     /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}_no_ma_GWAS_sumstats.txt \
     ${PATH_finemapping}/output/polyfun_finemap/polyfun_finemap.UKB.${chr}.${start}.${end} \
     ${PATH_finemapping}/output/polyfun_susie/polyfun_susie.UKB.${chr}.${start}.${end} \
     ${chr} ${start} ${end}
done

#Shared credible set variants between FINEMAP and Susie:
#rsid,chromosome,position,allele1,allele2
awk -F "\t" '{print $9","$10","$11","$12","$13}' ${PATH_finemapping}/output/susie/susie_credset.* | \
    grep -w -F -f - ${PATH_finemapping}/output/finemap/finemap_credset_*.txt | awk -F ":" '{print $2}' \
    > ${PATH_finemapping}/output/shared_credset_vars_fm_susie


awk -F "\t" '{print $9","$10","$11","$12","$13}' ${PATH_finemapping}/output/susie/susie_credset.* | \
    grep -w -F -f - ${PATH_finemapping}/output/finemap/finemap_credset_* | awk -F ":" '{print $1}' | uniq -c \
    > ${PATH_finemapping}/output/shared_credset_vars_fm_susie_N
