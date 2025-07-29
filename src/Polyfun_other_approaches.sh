#!/bin/bash

#PBS -N Polyfun_other_approaches
#PBS -j oe
#PBS -o Polyfun_other_approaches
#PBS -l walltime=1:0:0
#PBS -l vmem=50gb
#PBS -l nodes=1:ppn=1
#PBS -d .
#PBS -W umask=022

PATH_finemapping="/home/n/nnp5/PhD/PhD_project/Fine_mapping_severe_asthma"
PHENO="broad_pheno_1_5_ratio"
PATH_ASSOC="/home/n/nnp5/PhD/PhD_project/REGENIE_assoc/output/allchr"

module load gcc/9.3
module unload R/4.2.1
module load R/4.1.0

module load python/gcc/3.9.10
conda activate polyfun


#need to change these two index for each sentinel:
chr_row=5

chr=$(awk -v row=$chr_row 'NR == row {print $1}' /scratch/gen1/nnp5/Fine_mapping/tmp_data/fine_mapping_regions)
start=$(awk -v row=$chr_row 'NR == row {print $2}' /scratch/gen1/nnp5/Fine_mapping/tmp_data/fine_mapping_regions)
end=$(awk -v row=$chr_row 'NR == row {print $3}' /scratch/gen1/nnp5/Fine_mapping/tmp_data/fine_mapping_regions)

#Fine-mapping with SuSiE, using bgen file
#error:
#OSError: cannot load library '/cm/shared/apps/R/4.1.0/lib64/R/lib/libR.so': libRblas.so: cannot open shared object file: No such file or directory
#DEBUG with chatGPT's answer:This means that the path for R libraries is not present in my LD_LIBRARY_PATH for this session, so add it as follow:
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/cm/shared/apps/R/4.1.0/lib64/R/lib/
#install package required by polyfun:
# pip install bgen
python /home/n/nnp5/software/polyfun/finemapper.py \
    --geno /scratch/gen1/nnp5/Fine_mapping/tmp_data/severeasthma_EUR_ukb_imp_chr3_v3.bgen \
    --sample-file /scratch/gen1/nnp5/Fine_mapping/tmp_data/severeasthma_EUR_ukb_imp_chr3_v3.sample \
    --ldstore2 /home/n/nnp5/software/ldstore_v2.0_x86_64/ldstore_v2.0_x86_64 \
    --sumstats /scratch/gen1/nnp5/Fine_mapping/tmp_data/SNP_prior_manual.parquet \
    --n 46086 \
    --chr $chr \
    --start $start \
    --end $end \
    --method susie \
    --max-num-causal 10 \
    --allow-missing \
    --cache-dir /scratch/gen1/nnp5/Fine_mapping/LD_cache \
    --out ${PATH_finemapping}/output/test_frombgen_manual_finemap.UKB.$chr.$start.$end.gz

#Error: 'sqlite3.DatabaseError: malformed database schema (Variant) - near "WITHOUT": syntax error'