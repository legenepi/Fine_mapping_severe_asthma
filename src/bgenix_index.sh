#!/bin/bash

#PBS -N bgenix_index
#PBS -j oe
#PBS -o bgenix_index
#PBS -l walltime=36:0:0
#PBS -l vmem=100gb
#PBS -l nodes=1:ppn=1
#PBS -d .
#PBS -W umask=022

chr=$PBS_ARRAYID

#submit as: qsub -t 1-22 /home/n/nnp5/PhD/PhD_project/Fine_mapping_severe_asthma/src/bgenix_index.sh

cd /scratch/gen1/nnp5/Fine_mapping/tmp_data/

module load qctool
qctool -g /data/ukb/imputed_v3/ukb_imp_chr${chr}_v3.bgen \
    -s /data/gen1/UKBiobank_500K/severe_asthma/data/ukbiobank_app56607_for_regenie.sample \
    -incl-samples /home/n/nnp5/PhD/PhD_project/Post_GWAS/input/broadasthma_individuals \
    -og /scratch/gen1/nnp5/Fine_mapping/tmp_data/sevasthma_chr${chr}_v3.bgen

~nrgs1/bin/bgenix -g /scratch/gen1/nnp5/Fine_mapping/tmp_data/sevasthma_chr${chr}_v3.bgen -index