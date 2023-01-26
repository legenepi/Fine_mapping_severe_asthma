#!/bin/bash

#PBS -N pre_processing_Polyfun3
#PBS -j oe
#PBS -o pre_procesing_Polyfun3
#PBS -l walltime=2:0:0
#PBS -l vmem=300gb
#PBS -l nodes=1:ppn=1
#PBS -d .
#PBS -W umask=022
#Approach 3 (create LD-score and prior causal probs from a reference panel)

module load gcc/9.3
module unload R/4.2.1
module load R/4.1.0
conda activate polyfun
module load python/gcc/3.9.10

#Stage 1. Create a munged summary statistics file in a PolyFun-friendly parquet format.
#chmod o+x src/pre_processing_Polyfun.R
#dos2unix src/pre_processing_Polyfun.R
#Rscript src/pre_processing_Polyfun.R \
#    ${PATH_ASSOC}/${PHENO}_allchr.assoc.txt.gz \
#    ${PHENO}
#cut -d ' ' -f -11 /home/n/nnp5/PhD/PhD_project/Fine_mapping_severe_asthma/input/maf001_broad_pheno_1_5_ratio_betase_input_mungestat \
#    > maf001_broad_pheno_1_5_ratio_betase_input_mungestat
#python munge_polyfun_sumstats.py \
#  --sumstats maf001_broad_pheno_1_5_ratio_betase_input_mungestat \
#  --n 46086 \
#  --out /scratch/gen1/nnp5/Fine_mapping/tmp_data/maf001_broad_pheno_1_5_ratio_sumstats_munged.parquet \
#  --min-info 0.9 \
#  --min-maf 0.01

#Stage 2. Run PolyFun with L2-regularized S-LDSC - need to be done as a single job
python /home/n/nnp5/software/polyfun/polyfun.py \
    --compute-h2-L2 \
    --output-prefix /scratch/gen1/nnp5/Fine_mapping/tmp_data/L2_sldsc \
    --sumstats /scratch/gen1/nnp5/Fine_mapping/tmp_data/maf001_broad_pheno_1_5_ratio_sumstats_munged.parquet \
    --ref-ld-chr /scratch/gen1/nnp5/Fine_mapping/tmp_data/baselineLF2.2.UKB/baselineLF2.2.UKB. \
    --w-ld-chr /scratch/gen1/nnp5/Fine_mapping/tmp_data/baselineLF2.2.UKB/weights.UKB.

#Stage 3. Compute LD-scores for each SNP bin - need to run it as a job for each chr
##same --output-prefix of Stage 2, because PolyFun requires intermediate files that were created in stage 2
#chr=$PBS_ARRAYID
##need to create qc_ukb_cal_v2: genotypes, good quality, pruned and with unrelated individuals:
##filter variants according to QC measures and pruning, then filter out for broad asthma individuals

#python /home/n/nnp5/software/polyfun/polyfun.py \
#    --compute-ldscores \
#    --output-prefix /scratch/gen1/nnp5/Fine_mapping/tmp_data/L2_sldsc \
#    --bfile-chr /scratch/gen1/nnp5/REGENIE_assoc/tmp_data/qc_ukb_cal_v2 \
#    --keep /home/n/nnp5/PhD/PhD_project/Post_GWAS/input/broadasthma_individuals \
#    --chr ${chr}






#####Other options, not taken into cosideration right now######
#Computing LD-scores for annotations with your own pre-computed LD matrices
#Pre-Stag: Create LD matrices with LDstore 2.0:
#awk '{print $1}' /home/n/nnp5/PhD/PhD_project/Fine_mapping_severe_asthma/input/maf001_broad_pheno_1_5_ratio_betase_input_mungestat \
#    > /home/n/nnp5/PhD/PhD_project/Fine_mapping_severe_asthma/input/rsid_maf001_broad_pheno_1_5_ratio_betase_input_mungestat

#create the in-sample dosage LD file using ldstore2:
#chr=$PBS_ARRAYID
#module load plink2
#plink2 \
#    --bgen /data/ukb/nobackup/imputed_v3/ukb_imp_chr${chr}_v3.bgen \
#    ref-first --sample /data/gen1/UKBiobank_500K/severe_asthma/data/ukbiobank_app56607_for_regenie.sample \
#    --keep /data/gen1/UKBiobank_500K/severe_asthma/Noemi_PhD/data/broadasthma_individuals \
#    --export bgen-1.2 --out /scratch/gen1/nnp5/Fine_mapping/tmp_data/severeasthma_EUR_ukb_imp_chr${chr}_v3

#create the input files for ldstorev2:
#sample:
#awk '{print $1, $2, $3}' /scratch/gen1/nnp5/Fine_mapping/tmp_data/severeasthma_EUR_ukb_imp_chr${chr}_v3.sample \
#    > /scratch/gen1/nnp5/Fine_mapping/tmp_data/severeasthma_EUR_ukb_imp_chr${chr}_v3_bcor.sample
#bgen.bgi:
#/home/n/nnp5/software/bgen.tgz/build/apps/bgenix -index -g /scratch/gen1/nnp5/Fine_mapping/tmp_data/severeasthma_EUR_ukb_imp_chr${chr}_v3
#snp${chr}.z:
#/home/n/nnp5/software/bgen.tgz/build/apps/bgenix \
#    -g /scratch/gen1/nnp5/Fine_mapping/tmp_data/severeasthma_EUR_ukb_imp_chr${chr}_v3.bgen \
#    -list > /scratch/gen1/nnp5/Fine_mapping/tmp_data/tmp_snp_chr${chr}.z
#grep -v "#" /scratch/gen1/nnp5/Fine_mapping/tmp_data/tmp_snp_chr${chr}.z | awk -F " " '{print $2, $3, $4, $6, $7}' OFS=" " | tail -n +2 | sed -e '1i\rsid chromosome position allele1 allele2' \
#    > /scratch/gen1/nnp5/Fine_mapping/tmp_data/snp_chr${chr}.z
#master:
#echo "z;bgen;bgi;sample;bdose;bcor;ld;n_samples" > /scratch/gen1/nnp5/Fine_mapping/tmp_data/tmp_master_chr${chr}
#echo "/scratch/gen1/nnp5/Fine_mapping/tmp_data/snp_chr${chr}.z;/scratch/gen1/nnp5/Fine_mapping/tmp_data/severeasthma_EUR_ukb_imp_chr${chr}_v3.bgen;/scratch/gen1/nnp5/Fine_mapping/tmp_data/severeasthma_EUR_ukb_imp_chr${chr}_v3.bgen.bgi;/scratch/gen1/nnp5/Fine_mapping/tmp_data/severeasthma_EUR_ukb_imp_chr${chr}_v3_bcor.sample;/scratch/gen1/nnp5/Fine_mapping/tmp_data/chr${chr}.bdose;/scratch/gen1/nnp5/Fine_mapping/tmp_data/chr${chr}.bcor;/scratch/gen1/nnp5/Fine_mapping/tmp_data/chr${chr}.ld;46086" | \
#    cat /scratch/gen1/nnp5/Fine_mapping/tmp_data/tmp_master_chr${chr} - \
#    > /scratch/gen1/nnp5/Fine_mapping/tmp_data/master_chr${chr}

#ldstore2 bgen to bdose v1.1 conversion:
#/home/n/nnp5/software/ldstore_v2.0_x86_64/ldstore_v2.0_x86_64 \
#    --in-files /scratch/gen1/nnp5/Fine_mapping/tmp_data/master_chr${chr} \
#    --write-bdose --bdose-version 1.1

#ldstore2 SNPs correlation:
#/home/n/nnp5/software/ldstore_v2.0_x86_64/ldstore_v2.0_x86_64 \
#    --in-files /scratch/gen1/nnp5/Fine_mapping/tmp_data/master_chr${chr} \
#    --write-bcor \
#    --bdose-version 1.1 \
#    --read-only-bgen \
#    --rsids /home/n/nnp5/PhD/PhD_project/Fine_mapping_severe_asthma/input/rsid_maf001_broad_pheno_1_5_ratio_betase_input_mungestat

#cd /home/n/nnp5/software/polyfun
#conda activate polyfun
#module load python/gcc/3.9.10
#python compute_ldscores_from_ld.py \
#  --annot /scratch/gen1/nnp5/Fine_mapping/tmp_data/baselineLF2.2.UKB/baselineLF2.2.UKB.3.annot.parquet \
#  --out /scratch/gen1/nnp5/Fine_mapping/tmp_data/ldscores_SA_EUR_maf001_chr3.parquet \
#  --n 46086 \
#  bcor_files/*.bcor
#conda deactivate



#Computing LD-scores for annotation with pre-computed UK Biobank LD matrices:
#python /home/n/nnp5/software/polyfun/compute_ldscores_from_ld.py \
#  --annot /scratch/gen1/nnp5/Fine_mapping/tmp_data/baselineLF2.2.UKB/baselineLF2.2.UKB.3.annot.parquet \
#  --ukb \
#  --out /scratch/gen1/nnp5/Fine_mapping/tmp_data/ldscores_chr3.parquet