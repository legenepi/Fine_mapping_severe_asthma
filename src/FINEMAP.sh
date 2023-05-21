#!/bin/bash

#PBS -N FINEMAP
#PBS -j oe
#PBS -o FINEMAP
#PBS -l walltime=2:0:0
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

##data.bgen:
qsub -t 1-22 ${PATH_finemapping}/src/bgenix_index.sh

##data.bgen.bgi:
/scratch/gen1/nnp5/Fine_mapping/tmp_data/sevasthma_chr2_v3.bgen.bgi

##data.sample:
grep -w -F -f /home/n/nnp5/PhD/PhD_project/Post_GWAS/input/broadasthma_individuals \
    /data/gen1/UKBiobank_500K/severe_asthma/data/ukbiobank_app56607_for_regenie.sample | awk '{print $1, $2, $3}' \
    > ${PATH_finemapping}/input/ldstore.sample.tmp
echo "ID_1 ID_2 missing" > ${PATH_finemapping}/input/ldstore_header.sample
echo "0 0 0" > ${PATH_finemapping}/input/ldstore_secondline.sample
cat ${PATH_finemapping}/input/ldstore_header.sample ${PATH_finemapping}/input/ldstore_secondline.sample | \
    cat - ${PATH_finemapping}/input/ldstore.sample.tmp \
    > ${PATH_finemapping}/input/ldstore.sample
rm ${PATH_finemapping}/input/ldstore.sample.tmp ${PATH_finemapping}/input/ldstore_header.sample \
    ${PATH_finemapping}/input/ldstore_secondline.sample

##sevasthma.z: space-delimited text file
#rsid chromosome position allele1 allele2 maf beta se
awk '{print $1, $2, $3, $5, $4, $12, $6, $7}' \
    ${PATH_finemapping}/input/maf001_broad_pheno_1_5_ratio_betase_input_mungestat |
    awk '
    {if($2==1) $2 = "01"
     if($2==2) $2 = "02"
     if($2==3) $2 = "03"
     if($2==4) $2 = "04"
     if($2==5) $2 = "05"
     if($2==6) $2 = "06"
     if($2==7) $2 = "07"
     if($2==8) $2 = "08"
     if($2==9) $2 = "09"
     }
     1' | \
    sed "1s/.*/rsid chromosome position allele1 allele2 maf beta se/" \
    > ${PATH_finemapping}/input/sevasthma.z

##dataset.incl:
awk {'print $1'} /home/n/nnp5/PhD/PhD_project/Post_GWAS/input/broadasthma_individuals | tail -n +2 \
    > ${PATH_finemapping}/input/ldstore.incl

##SEVASTHMA.BCOR: as output from LDstore2 BCOR v1.1
#Need to be created with LDSTORE2:
#z;bgen;bgi;sample;bdose;bcor;ld;n_samples;incl
#data.z;/data.bgen;example/data.bgen.bgi;example/data.sample;example/data.bdose;example/data.bcor;example/data.ld;46086;input/ldstore.incl

##data.z:
#rsid chromosome position allele1 allele2
for line in {2..14}
do
SNP=$(awk -v row="$line" ' NR == row {print $1 } ' ${PATH_finemapping}/input/fine_mapping_regions_merged)
chr=$(awk -v row="$line" ' NR == row {print $2 } ' ${PATH_finemapping}/input/fine_mapping_regions_merged)
start=$(awk -v row="$line" 'NR == row {print $4}' ${PATH_finemapping}/input/fine_mapping_regions_merged)
end=$(awk -v row="$line" 'NR == row {print $5}' ${PATH_finemapping}/input/fine_mapping_regions_merged)

awk -v chr_idx=$chr 'NR==1; NR > 1 {if ($2 == chr_idx) print}'  ${PATH_finemapping}/input/sevasthma.z | \
    awk -v START_POS="$start" 'NR==1; NR > 1 {if ($3 >= START_POS) print}' | \
    awk -v END_POS="$end" 'NR==1; NR > 1 {if ($3 <= END_POS) print}' \
    > ${PATH_finemapping}/input/ldstore_chr${chr}_${SNP}.z

#master file:
echo "z;bgen;bgi;sample;bdose;bcor;ld;n_samples;incl" > ${PATH_finemapping}/input/ldstore_chr${chr}_${SNP}.data
echo "input/ldstore_chr${chr}_${SNP}.z;/scratch/gen1/nnp5/Fine_mapping/tmp_data/sevasthma_chr${chr}_v3.bgen;/scratch/gen1/nnp5/Fine_mapping/tmp_data/sevasthma_chr${chr}_v3.bgen.bgi;input/ldstore.sample;input/ldstore_chr${chr}_${SNP}.bdose;input/ldstore_chr${chr}_${SNP}.bcor;input/ldstore_chr${chr}_${SNP}.ld;46086;input/ldstore.incl" \
    >> ${PATH_finemapping}/input/ldstore_chr${chr}_${SNP}.data

#ldstore2:
/home/n/nnp5/software/ldstore_v2.0_x86_64/ldstore_v2.0_x86_64 \
    --in-files ${PATH_finemapping}/input/ldstore_chr${chr}_${SNP}.data \
    --write-bcor \
    --write-bdose \
    --bdose-version 1.1

#master file: semicolon-delimiter text file:
# #NB:The order of the SNPs in the dataset.ld must correspond to the order of SNPs in dataset.z.
cd ${PATH_finemapping}
echo "z;bcor;snp;config;cred;log;n_samples" > ${PATH_finemapping}/input/finemap_chr${chr}_${SNP}.z
echo "input/ldstore_chr${chr}_${SNP}.z;input/ldstore_chr${chr}_${SNP}.bcor;output/finemap_${chr}_${SNP}.snp;output/finemap_${chr}_${SNP}.config;output/finemap_${chr}_${SNP}.cred;output/finemap_${chr}_${SNP}.log;46086" \
    >> ${PATH_finemapping}/input/finemap_chr${chr}_${SNP}.z


#finemap:
/home/n/nnp5/software/finemap_v1.4.1_x86_64/finemap_v1.4.1_x86_64 \
    --sss \
    --n-causal-snps 10 \
    --in-files ${PATH_finemapping}/input/finemap_chr${chr}_${SNP}.z \
    --dataset 1

done

#post analysis: Find the credible set variants as per FINEMAP:
grep -v "^#"  /home/n/nnp5/PhD/PhD_project/Fine_mapping_severe_asthma/output/finemap_*_plink.cred1 \
    > /home/n/nnp5/PhD/PhD_project/Fine_mapping_severe_asthma/output/finemap_plink.cred1.digest


#six columns interpreted as SNPID, rsid, chromosome, position, first and second alleles.