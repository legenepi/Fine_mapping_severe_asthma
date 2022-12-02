#!/bin/bash

#PBS -N pre_processing_Polyfun
#PBS -j oe
#PBS -o pre_procesing_Polyfun
#PBS -l walltime=12:0:0
#PBS -l vmem=50gb
#PBS -l nodes=1:ppn=1
#PBS -d .
#PBS -W umask=022
#Approach 3 (create LD-score and prior causal probs from a reference panel)
#Create LD matrices with LDstore 2.0:

#awk '{print $1}' /home/n/nnp5/PhD/PhD_project/Fine_mapping_severe_asthma/input/maf001_broad_pheno_1_5_ratio_betase_input_mungestat \
#    > /home/n/nnp5/PhD/PhD_project/Fine_mapping_severe_asthma/input/rsid_maf001_broad_pheno_1_5_ratio_betase_input_mungestat

chr=$PBS_ARRAYID
module load plink2
plink2 \
    --bgen /data/ukb/nobackup/imputed_v3/ukb_imp_chr${chr}_v3.bgen \
    ref-first --sample /data/gen1/UKBiobank_500K/severe_asthma/data/ukbiobank_app56607_for_regenie.sample \
    --keep /data/gen1/UKBiobank_500K/severe_asthma/Noemi_PhD/data/broadasthma_individuals \
    --export bgen-1.2 --out /scratch/gen1/nnp5/Fine_mapping/tmp_data/severeasthma_EUR_ukb_imp_chr${chr}_v3

#/home/n/nnp5/software/ldstore_v2.0_x86_64/ldstore_v2.0_x86_64 \
#    --in-files /scratch/gen1/nnp5/Fine_mapping/tmp_data/severeasthma_EUR_ukb_imp_chr3_v3 \
#    --write-bcor \
#    --read-only-bgen \
#    --rsids /home/n/nnp5/PhD/PhD_project/Fine_mapping_severe_asthma/input/rsid_maf001_broad_pheno_1_5_ratio_betase_input_mungestat

#Computing LD-scores with your own pre-computed LD matrices
#python compute_ldscores_from_ld.py \
#  --annot /scratch/gen1/nnp5/Fine_mapping/tmp_data/baselineLF2.2.UKB/baselineLF2.2.UKB.3.annot.parquet \
#  --out /scratch/gen1/nnp5/Fine_mapping/tmp_data/ldscores_SA_EUR_maf001_chr3.parquet \
#  --n 46086 \
#  bcor_files/*.bcor

#Run PolyFun with L2-regularized S-LDSC
#python polyfun.py \
#    --compute-h2-L2 \
#    --output-prefix output/testrun \
#    --sumstats maf001_broad_pheno_1_5_ratio_sumstats_munged.parquet \
#    --ref-ld-chr /scratch/gen1/nnp5/Fine_mapping/tmp_data/ \
#    --w-ld-chr /scratch/gen1/nnp5/Fine_mapping/tmp_data/