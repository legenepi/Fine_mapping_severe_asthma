#!/bin/bash

#Rationale: take credible set variants that shows a similar posterior inclusion probability
#in both SuSiE and FINEMAP.

#analyse the variants that are shared:
PATH_finemapping="/home/n/nnp5/PhD/PhD_project/Fine_mapping_severe_asthma"
PATH_ASSOC="/home/n/nnp5/PhD/PhD_project/REGENIE_assoc/output/allchr"

module load gcc/9.3
module unload R/4.2.1
module load R/4.1.0

cd ${PATH_finemapping}

Rscript Valid_credset_shared_variants.R


#susie only in finemap:
awk '{print $1}' output/susie_only_credset_vars | grep -h -w -F -f - output/finemap/finemap_*.snp \
    > output/susie_only_credset_vars_in_finemap


#finemap only in susie:
awk '{print $1}' output/finemap_only_credset_vars > output/finemap_only_credset_vars_snp

touch output/finemap_only_credset_vars_in_susie

ls -lthr /scratch/gen1/nnp5/Fine_mapping/tmp_data/*_no_ma_GWAS_sumstats.txt | \
    awk -F '/|_no' '{print $7}' | grep -v 'rs705705_rs3024971' \
    > /scratch/gen1/nnp5/Fine_mapping/tmp_data/rsid


#header rsid chromosome positionallele1 allele2 vars.variable_prob: 'rsid chromosome positionallele1 allele2 vars.variable_prob'

while read -r line; do paste /scratch/gen1/nnp5/Fine_mapping/tmp_data/${line}_no_ma_GWAS_sumstats.txt \
    output/susie/susie_*_${line}*.txt | awk -F ' |\t' '{print $1, $2, $3, $4, $5, $10}' | grep -w -F -f output/finemap_only_credset_vars_snp - \
    >> output/finemap_only_credset_vars_in_susie; done < /scratch/gen1/nnp5/Fine_mapping/tmp_data/rsid

Rscript Valid_credset_NOTshared_variants.R
