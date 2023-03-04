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
chr_row=1

chr=$(awk -v row=$chr_row 'NR == row {print $1}' /scratch/gen1/nnp5/Fine_mapping/tmp_data/fine_mapping_regions)
start=$(awk -v row=$chr_row 'NR == row {print $2}' /scratch/gen1/nnp5/Fine_mapping/tmp_data/fine_mapping_regions)
end=$(awk -v row=$chr_row 'NR == row {print $3}' /scratch/gen1/nnp5/Fine_mapping/tmp_data/fine_mapping_regions)

#Fine-mapping with SuSiE, using genotypes from a bgen file
#python /home/n/nnp5/software/polyfun/finemapper.py \
#    --geno /scratch/gen1/nnp5/Fine_mapping/tmp_data/severeasthma_EUR_ukb_imp_chr2_v3.bgen \
#    --sample-file /scratch/gen1/nnp5/Fine_mapping/tmp_data/severeasthma_EUR_ukb_imp_chr2_v3.sample \
#    --ldstore2 /home/n/nnp5/software/ldstore_v2.0_x86_64/ldstore_v2.0_x86_64 \
#    --sumstats /scratch/gen1/nnp5/Fine_mapping/tmp_data/SNP_prior_manual.parquet \
#    --n 46086 \
#    --chr $chr \
#    --start $start \
#    --end $end \
#    --method susie \
#    --max-num-causal 10 \
#    --allow-missing \
#    --cache-dir /scratch/gen1/nnp5/Fine_mapping/LD_cache \
#    --out ${PATH_finemapping}/output/test_fromplink_manual_finemap.UKB.$chr.$start.$end.gz

#error:
#[DEBUG]  opening BgenFile from /scratch/gen1/nnp5/Fine_mapping/tmp_data/severeasthma_EUR_ukb_imp_chr2_v3.bgen
#Traceback (most recent call last):
#  File "/home/n/nnp5/software/polyfun/finemapper.py", line 1275, in <module>
#    df_finemap = finemap_obj.finemap(locus_start=args.start, locus_end=args.end, num_causal_snps=args.max_num_causal,
#  File "/home/n/nnp5/software/polyfun/finemapper.py", line 739, in finemap
#    ld_arr, df_ld_snps = self.get_ld_data(locus_start, locus_end, need_bcor=False, verbose=verbose)
#  File "/home/n/nnp5/software/polyfun/finemapper.py", line 579, in get_ld_data
#    ld_file = self.compute_ld_bgen(locus_start, locus_end, verbose=verbose)
#  File "/home/n/nnp5/software/polyfun/finemapper.py", line 409, in compute_ld_bgen
#    for snp_i, rsid in enumerate(rsids):
#  File "src/bgen/bgen.pyx", line 384, in fetch
#ValueError: can't fetch variants without index

#Fine-mapping with SuSiE, using genotypes from a bgen file
python /home/n/nnp5/software/polyfun/finemapper.py \
    --geno /scratch/gen1/nnp5/REGENIE_assoc/tmp_data/broad_pheno_plink_file_v3_chr2 \
    --sumstats /scratch/gen1/nnp5/Fine_mapping/tmp_data/SNP_prior_manual.parquet \
    --n 46086 \
    --chr $chr \
    --start $start \
    --end $end \
    --method susie \
    --max-num-causal 10 \
    --allow-missing \
    --cache-dir /scratch/gen1/nnp5/Fine_mapping/LD_cache \
    --out ${PATH_finemapping}/output/test_fromplink_manual_finemap.UKB.$chr.$start.$end.gz