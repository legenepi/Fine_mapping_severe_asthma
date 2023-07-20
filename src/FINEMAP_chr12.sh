#!/bin/bash

#PBS -N FINEMAP
#PBS -j oe
#PBS -o FINEMAP
#PBS -l walltime=4:0:0
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

cd ${PATH_finemapping}

##SEVASTHMA.BCOR: as output from LDstore2 BCOR v1.1
#Need to be created with LDSTORE2:
#z;bgen;bgi;sample;bdose;bcor;ld;n_samples;incl
#data.z;/data.bgen;example/data.bgen.bgi;example/data.sample;example/data.bdose;example/data.bcor;example/data.ld;46086;input/ldstore.incl

##data.z:
#rsid chromosome position allele1 allele2
#input/fine_mapping_regions_merged from create_fine_mapping_merged.sh
for line in {2..3}
do
SNP=$(awk -v row="$line" ' NR == row {print $1 } ' ${PATH_finemapping}/input/fine_mapping_regions_chr12)
chr=$(awk -v row="$line" ' NR == row {print $2 } ' ${PATH_finemapping}/input/fine_mapping_regions_chr12)
start=$(awk -v row="$line" 'NR == row {print $4}' ${PATH_finemapping}/input/fine_mapping_regions_chr12)
end=$(awk -v row="$line" 'NR == row {print $5}' ${PATH_finemapping}/input/fine_mapping_regions_chr12)

awk -v chr_idx=$chr 'NR==1; NR > 1 {if ($2 == chr_idx) print}'  ${PATH_finemapping}/input/sevasthma.z | \
    awk -v START_POS="$start" 'NR==1; NR > 1 {if ($3 >= START_POS) print}' | \
    awk -v END_POS="$end" 'NR==1; NR > 1 {if ($3 <= END_POS) print}' \
    > ${PATH_finemapping}/input/ldstore_chr${chr}_${SNP}.z

#master file for ldstore2:
echo "z;bgen;bgi;sample;bdose;bcor;ld;n_samples;incl" > ${PATH_finemapping}/input/ldstore_chr${chr}_${SNP}.data
echo "input/ldstore_chr${chr}_${SNP}.z;/scratch/gen1/nnp5/Fine_mapping/tmp_data/sevasthma_chr${chr}_v3.bgen;/scratch/gen1/nnp5/Fine_mapping/tmp_data/sevasthma_chr${chr}_v3.bgen.bgi;input/ldstore.sample;input/ldstore_chr${chr}_${SNP}.bdose;input/ldstore_chr${chr}_${SNP}.bcor;input/ldstore_chr${chr}_${SNP}.ld;46086;input/ldstore.incl" \
    >> ${PATH_finemapping}/input/ldstore_chr${chr}_${SNP}.data

#ldstore2:
/home/n/nnp5/software/ldstore_v2.0_x86_64/ldstore_v2.0_x86_64 \
    --in-files ${PATH_finemapping}/input/ldstore_chr${chr}_${SNP}.data \
    --write-bcor \
    --write-bdose \
    --bdose-version 1.1

#master file for FINEMAP: semicolon-delimiter text file:
# #NB:The order of the SNPs in the dataset.ld must correspond to the order of SNPs in dataset.z.
echo "z;bcor;snp;config;cred;log;n_samples" > ${PATH_finemapping}/input/finemap_chr${chr}_${SNP}.z
echo "input/ldstore_chr${chr}_${SNP}.z;input/ldstore_chr${chr}_${SNP}.bcor;output/finemap/finemap_${chr}_${SNP}_${start}_${end}.snp;output/finemap/finemap_${chr}_${SNP}_${start}_${end}.config;output/finemap/finemap_${chr}_${SNP}_${start}_${end}.cred;output/finemap/finemap_${chr}_${SNP}_${start}_${end}.log;46086" \
    >> ${PATH_finemapping}/input/finemap_chr${chr}_${SNP}.z

#finemap:
/home/n/nnp5/software/finemap_v1.4.1_x86_64/finemap_v1.4.1_x86_64 \
    --sss \
    --n-causal-snps 10 \
    --in-files ${PATH_finemapping}/input/finemap_chr${chr}_${SNP}.z \
    --log

done

#six columns interpreted as SNPID, rsid, chromosome, position, first and second alleles.

#Plot in R to compare GWAS p-value and fine-mapping PIP:
#from this website: https://www.mv.helsinki.fi/home/mjxpirin/GWAS_course/material/GWAS7.html
locus="12_rs705705_55935504_56935504"
Rscript src/finemap_plot.R $locus
locus="12_rs3024971_56993727_57993727"
Rscript src/finemap_plot.R $locus

#Merge credset into a unique file for Finemapping.xlsx in Report:
head -n 1 ${PATH_finemapping}/output/finemap/finemap_credset_2_rs12470864_101926362_103926362.txt \
    > ${PATH_finemapping}/output/finemap/finemap_all_credset.txt && \
    tail -n +2 -q ${PATH_finemapping}/output/finemap/finemap_credset_*.txt >> ${PATH_finemapping}/output/finemap/finemap_all_credset.txt

