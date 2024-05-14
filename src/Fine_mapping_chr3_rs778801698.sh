#!/bin/bash
REDO, I DO NOT LIKE THE TOTAL NUMBER OF CAUSAL VARIATNS TO FOLLOW-UP: 264... TOO HIGH?

#Rationale: Fine-mapping for additional replicated sentinel rs778801698 and its genomic loci.

##load required tools:
module load gcc/12.3.0-yxgv2bl
module load plink2
module load R

PATH_finemapping="/home/n/nnp5/PhD/PhD_project/Fine_mapping_severe_asthma"
cd ${PATH_finemapping}/

#create the input file for fine-mapping with the region (1Mb centered on the sentinel):
echo "rs778801698 3 50024027 49524027 50524027" \
    > ${PATH_finemapping}/fine_mapping_region_input_chr3_rs778801698

#bgen file:
mkdir /scratch/gen1/nnp5/Fine_mapping/tmp_data
sbatch --array=3 src/bgenix_index_sbatch.sh

#FINEMAP:
##data.sample: ${PATH_finemapping}/input/ldstore.sample

##sevasthma.z: space-delimited text file: ${PATH_finemapping}/input/sevasthma.z
#rsid chromosome position allele1 allele2 maf beta se

##dataset.incl: ${PATH_finemapping}/input/ldstore.incl

##SEVASTHMA.BCOR: as output from LDstore2 BCOR v1.1
#Need to be created with LDSTORE2:
#z;bgen;bgi;sample;bdose;bcor;ld;n_samples;incl
#data.z;/data.bgen;example/data.bgen.bgi;example/data.sample;example/data.bdose;example/data.bcor;example/data.ld;46086;input/ldstore.incl

##data.z:
#rsid chromosome position allele1 allele2
SNP=$(awk -v row=1 ' NR == row {print $1 } ' ${PATH_finemapping}/fine_mapping_region_input_chr3_rs778801698)
chr=$(awk -v row=1 ' NR == row {print $2 } ' ${PATH_finemapping}/fine_mapping_region_input_chr3_rs778801698)
start=$(awk -v row=1 'NR == row {print $4}' ${PATH_finemapping}/fine_mapping_region_input_chr3_rs778801698)
end=$(awk -v row=1 'NR == row {print $5}' ${PATH_finemapping}/fine_mapping_region_input_chr3_rs778801698)


##remove the multiallelic SNPs, as done for SuSiE, so start from the same number of variants:
awk -v chr_idx=$chr 'NR==1; NR > 1 {if ($2 == chr_idx) print}' ${PATH_finemapping}/input/sevasthma.z | \
    awk -v START_POS="$start" 'NR==1; NR > 1 {if ($3 >= START_POS) print}' | \
    awk -v END_POS="$end" 'NR==1; NR > 1 {if ($3 <= END_POS) print}' | \
    grep -v -w -F -f /data/gen1/UKBiobank_500K/imputed/multiallelic.snps - \
    > ${PATH_finemapping}/input/ldstore_chr${chr}_${SNP}.z

##master file for ldstore2:
echo "z;bgen;bgi;sample;bdose;bcor;ld;n_samples;incl" > ${PATH_finemapping}/input/ldstore_chr${chr}_${SNP}.data
echo "input/ldstore_chr${chr}_${SNP}.z;/scratch/gen1/nnp5/Fine_mapping/tmp_data/sevasthma_chr${chr}_v3.bgen;/scratch/gen1/nnp5/Fine_mapping/tmp_data/sevasthma_chr${chr}_v3.bgen.bgi;input/ldstore.sample;input/ldstore_chr${chr}_${SNP}.bdose;input/ldstore_chr${chr}_${SNP}.bcor;input/ldstore_chr${chr}_${SNP}.ld;46086;input/ldstore.incl" \
    >> ${PATH_finemapping}/input/ldstore_chr${chr}_${SNP}.data

##ldstore2:
/home/n/nnp5/software/ldstore_v2.0_x86_64/ldstore_v2.0_x86_64 \
    --in-files ${PATH_finemapping}/input/ldstore_chr${chr}_${SNP}.data \
    --write-bcor \
    --write-bdose \
    --bdose-version 1.1

##master file for FINEMAP: semicolon-delimiter text file:
###NB:The order of the SNPs in the dataset.ld must correspond to the order of SNPs in dataset.z.
echo "z;bcor;snp;config;cred;log;n_samples" > ${PATH_finemapping}/input/finemap_chr${chr}_${SNP}.z
echo "input/ldstore_chr${chr}_${SNP}.z;input/ldstore_chr${chr}_${SNP}.bcor;output/finemap_replicated_suggestive/finemap_replsugg_${chr}_${SNP}_${start}_${end}.snp;output/finemap_replicated_suggestive/finemap_replsugg_${chr}_${SNP}_${start}_${end}.config;output/finemap_replicated_suggestive/finemap_replsugg_${chr}_${SNP}_${start}_${end}.cred;output/finemap_replicated_suggestive/finemap_replsugg_${chr}_${SNP}_${start}_${end}.log;46086" \
    >> ${PATH_finemapping}/input/finemap_chr${chr}_${SNP}.z

/home/n/nnp5/software/finemap_v1.4.1_x86_64/finemap_v1.4.1_x86_64 \
    --sss \
    --n-causal-snps 10 \
    --in-files ${PATH_finemapping}/input/finemap_chr${chr}_${SNP}.z \
    --log

##Find credset:
Rscript ${PATH_finemapping}/src/finemapping_replicated_suggestive/credset_FINEMAP.R ${chr}_${SNP}_${start}_${end}

#because I did this analysis two times and i wanted to be sure that the second time I wasn't adding the credset twice:
grep -v "3_rs778801698_49524027_50524027" ${PATH_finemapping}/output/finemap_replsugg_all_credset.txt \
    > ${PATH_finemapping}/output/finemap_replsugg_all_credset_rs778801698.txt

tail -n +2 -q ${PATH_finemapping}/output/finemap_replicated_suggestive/finemap_replsugg_credset_3_rs778801698_49524027_50524027.txt \
    >> ${PATH_finemapping}/output/finemap_replsugg_all_credset_rs778801698.txt




#SuSiE:
PATH_OUT="/home/n/nnp5/PhD/PhD_project/Fine_mapping_severe_asthma/output/susie_replicated_suggestive"
##Creating region bgen
conda activate /home/n/nnp5/miniconda3/envs/phd_env/
cd /scratch/gen1/nnp5/Fine_mapping/tmp_data/
bgenix -g /scratch/gen1/nnp5/Fine_mapping/tmp_data/sevasthma_chr${chr}_v3.bgen \
    -incl-range 0${chr}:${start}-${end} \
    > /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}.bgen
conda deactivate
cd ${PATH_finemapping}

##Exclude multi-allelic variants and find the common SNP IDs for the genotyped matrix and the zscore input files:
#use the file for each regions created by FINEMAP.sh:
grep -v -w -F -f /data/gen1/UKBiobank_500K/imputed/multiallelic.snps \
    ${PATH_finemapping}/input/ldstore_chr${chr}_${SNP}.z \
    > /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}_no_ma_GWAS_sumstats.txt

awk 'NR > 1 {print $1}' /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}_no_ma_GWAS_sumstats.txt \
    > /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}_no_ma_snps.txt

##zscore:
Rscript src/z_score.R /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}_no_ma_GWAS_sumstats.txt \
    /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}_no_ma_z_scores.txt

##Format region data for input to R
/home/n/nnp5/software/plink2 \
    --bgen /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}.bgen ref-first \
    --sample ${PATH_finemapping}/input/ldstore.sample \
    --export A \
    --extract /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}_no_ma_snps.txt \
    --out /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}

cut -f7- /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}.raw \
    > /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}.cols.raw

awk 'NR>1 {print}' /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}.cols.raw \
    > /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}.cols_no_header.raw

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
    ${chr}_${SNP}_${start}_${end} \
    ${PATH_OUT}/susie_replsugg_credset.${SNP}.$chr.$start.$end

#because I did this analysis two times and i wanted to be sure that the second time I wasn't adding the credset twice:
grep -v "3_rs778801698_49524027_50524027" ${PATH_finemapping}/output/susie_replsugg_all_credset.txt \
    > ${PATH_finemapping}/output/susie_replsugg_all_credset_rs778801698.txt

tail -n +2 -q ${PATH_OUT}/susie_replsugg_credset.rs778801698.3.49524027.50524027 \
    >> ${PATH_finemapping}/output/susie_replsugg_all_credset_rs778801698.txt

#Valid credset SNPs:
##Shared
Rscript src/finemapping_replicated_suggestive/Valid_shared_credset_variants_replicated_suggestive.R \
output/susie_replsugg_all_credset_rs778801698.txt output/finemap_replsugg_all_credset_rs778801698.txt output/replsugg_shared_valid_credset_chr3.txt

awk '{print $1}' output/replsugg_shared_valid_credset_chr3.txt | sort | uniq -c

##FINEMAP-ONLY:
PATH_OUT="/home/n/nnp5/PhD/PhD_project/Fine_mapping_severe_asthma/output/"
awk -v OFS=',' '{print $1, $2, $3, $4, $5, $6}' ${PATH_OUT}/susie_replsugg_all_credset_rs778801698.txt | \
    grep -v -w -F -f - ${PATH_OUT}/finemap_replsugg_all_credset_rs778801698.txt \
    > ${PATH_OUT}/finemap_replsugg_only_credset

awk -F "," -v OFS=' ' '{print $2, "0"$3, $4, $5, $6}' ${PATH_OUT}/finemap_replsugg_only_credset | \
    grep -o -n -w -F -f - /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}_no_ma_GWAS_sumstats.txt | \
    awk -F ":" '{print "V"$1, $2}' > /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}_finemap_only_vars_susie_idx

grep -v "3_rs778801698_49524027_50524027" ${PATH_OUT}/repl_sugg_NOTshared_finemaponly_valid_credset.txt \
    > ${PATH_OUT}/repl_sugg_NOTshared_finemaponly_valid_credset_rs778801698.txt
mv ${PATH_OUT}/repl_sugg_NOTshared_finemaponly_valid_credset_rs778801698.txt \
    ${PATH_OUT}/repl_sugg_NOTshared_finemaponly_valid_credset.txt

Rscript src/finemapping_replicated_suggestive/Valid_NOTshared_finemaponly_credset_variants_replicated_suggestive.R \
    ${PATH_OUT}/finemap_replicated_suggestive/finemap_replsugg_credset_${chr}_${SNP}_${start}_${end}.txt \
    /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}_finemap_only_vars_susie_idx \
    ${PATH_OUT}/susie_replicated_suggestive/susie_replsugg_${chr}_${SNP}_${start}_${end}.txt
#how many valid finemap-only credset SNPs for locus ?
awk '{print $1}' ${PATH_OUT}/repl_sugg_NOTshared_finemaponly_valid_credset.txt | sort | uniq -c

##SuSiE only:
awk -F "," -v OFS='\t' '{print $1, $2, $3, $4, $5, $6}' ${PATH_OUT}/finemap_replsugg_all_credset_rs778801698.txt | \
    grep -v -w -F -f - ${PATH_OUT}/susie_replsugg_all_credset_rs77880169.txt \
    > ${PATH_OUT}/susie_replsugg_only_credset


grep -v "3_rs778801698_49524027_50524027" ${PATH_OUT}/repl_sugg_NOTshared_susieonly_valid_credset.txt \
    > ${PATH_OUT}/repl_sugg_NOTshared_susieonly_valid_credset_rs778801698.txt
mv ${PATH_OUT}/repl_sugg_NOTshared_susieonly_valid_credset_rs778801698.txt \
    ${PATH_OUT}/repl_sugg_NOTshared_susieonly_valid_credset.txt

Rscript src/finemapping_replicated_suggestive/Valid_NOTshared_susieonly_credset_variants_replicated_suggestive.R \
    ${PATH_OUT}/susie_replsugg_only_credset \
    ${PATH_OUT}/finemap_replicated_suggestive/finemap_replsugg_${chr}_${SNP}_${start}_${end}.snp
#how many valid finemap-only credset SNPs for locus ?
awk '{print $1}' ${PATH_OUT}/repl_sugg_NOTshared_susieonly_valid_credset.txt | sort | uniq -c


cat ${PATH_finemapping}/output/replsugg_shared_valid_credset_chr3.txt ${PATH_finemapping}/output/repl_sugg_NOTshared_finemaponly_valid_credset.txt \
    > ${PATH_finemapping}/output/replsugg_valid_credset.txt.tmp
cat ${PATH_finemapping}/output/replsugg_valid_credset.txt.tmp ${PATH_finemapping}/output/repl_sugg_NOTshared_susieonly_valid_credset.txt \
    > ${PATH_finemapping}/output/replsugg_valid_credset_chr3.txt
rm ${PATH_finemapping}/output/replsugg_valid_credset.txt.tmp
grep -w -v "6_rs9271365_32086794_33086794" ${PATH_finemapping}/output/replsugg_valid_credset_chr3.txt \
    > ${PATH_finemapping}/output/replsugg_valid_credset_chr3_noMHC.txt

cp ${PATH_finemapping}/output/replsugg_valid_credset_chr3_noMHC.txt /data/gen1/UKBiobank_500K/severe_asthma/Noemi_PhD/data/