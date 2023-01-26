#!/bin/bash

#PBS -N pre_processing_Polyfun
#PBS -j oe
#PBS -o pre_procesing_Polyfun
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
chmod o+x src/pre_processing_Polyfun.R
dos2unix src/pre_processing_Polyfun.R
Rscript src/pre_processing_Polyfun.R \
    ${PATH_ASSOC}/${PHENO}_allchr.assoc.txt.gz \
    ${PHENO}

#iterative:
module load python/gcc/3.9.10
conda activate polyfun
#install parquet-tools to read parquet file in bash: python3 -m pip install parquet-tools

#Step1. Create a munged summary statistics file in a PolyFun-friendly parquet format.
cut -d ' ' -f -11 ${PATH_finemapping}/input/maf001_broad_pheno_1_5_ratio_betase_input_mungestat \
    > /scratch/gen1/nnp5/Fine_mapping/tmp_data/maf001_broad_pheno_1_5_ratio_betase_input_mungestat
python /home/n/nnp5/software/polyfun/munge_polyfun_sumstats.py \
  --sumstats /scratch/gen1/nnp5/Fine_mapping/tmp_data/maf001_broad_pheno_1_5_ratio_betase_input_mungestat \
  --n 46086 \
  --out /scratch/gen1/nnp5/Fine_mapping/tmp_data/maf001_broad_pheno_1_5_ratio_sumstats_munged.parquet \
  --min-info 0.85 \
  --min-maf 0.01

#Approach 1 (UK Biobank pre-computed prior causal probabilities); no information for HLA region !
##Extract pre-computed prior causal probabilities in UK Biobank - White British:
python /home/n/nnp5/software/polyfun/extract_snpvar.py \
    --sumstats /scratch/gen1/nnp5/Fine_mapping/tmp_data/maf001_broad_pheno_1_5_ratio_sumstats_munged.parquet \
    --allow-missing \
    --out /scratch/gen1/nnp5/Fine_mapping/tmp_data/SNPs_PriCauPro

#Manually retrieve the prior in python:
import pandas as pd
sumstat = pd.read_parquet("/scratch/gen1/nnp5/Fine_mapping/tmp_data/maf001_broad_pheno_1_5_ratio_sumstats_munged.parquet")
polyfun_snp=pd.read_parquet("/home/n/nnp5/software/polyfun/snpvar_meta.chr1_7.parquet")
polyfun_snp2=pd.read_parquet("/home/n/nnp5/software/polyfun/snpvar_meta.chr8_22.parquet")
df = pd.merge(sumstat, polyfun_snp, how='inner', on = 'SNP')
df2 = pd.merge(sumstat, polyfun_snp2, how='inner', on = 'SNP')
frames = [df, df2]
result = pd.concat(frames)
result['Z'].where(~(result['A1_x'] == result['A2_y']), other=-result['Z'], inplace=True)
result_clean = result[["CHR_y", "BP_y", "SNP", "A1_x", "A2_x", "MAF", "N", "Z", "snpvar_bin"]]
result_clean.columns = ['CHR','BP','SNP','A1','A2','MAF','N','Z','SNPVAR']
nodup_result_clean = result_clean.sort_values('SNP').drop_duplicates(subset=['CHR', 'BP'], keep='first')
nodup_result_clean.to_parquet("/scratch/gen1/nnp5/Fine_mapping/tmp_data/SNP_prior_manual.parquet")
#find sumstat SNPs that are not have SNPVAR:
#perform outer join
outer = sumstat.merge(nodup_result_clean, how='outer', on = 'SNP', indicator=True)
#perform anti-join
anti_join = outer[(outer._merge=='left_only')].drop('_merge', axis=1)
anti_join = anti_join[["SNP","CHR_x","BP_x","A1_x","A2_x","MAF_x","N_x"]]
anti_join.to_csv("/scratch/gen1/nnp5/Fine_mapping/tmp_data/missed_SNP_prior_manual.txt", sep="\t")


##Run the FIFO with SuSiE:
#OptionB.Fine-mapping with SuSiE, using pre-computed summary LD information from the UK Biobank
#download an LD matrix
mkdir -p LD_cache
mkdir output
mkdir /scratch/gen1/nnp5/Fine_mapping/tmp_data/LD_temp
cd /scratch/gen1/nnp5/Fine_mapping/tmp_data/LD_temp
#dowload for each sentinel the appropriate region:
#chr2 "2" "102426362" "103426362":
wget https://data.broadinstitute.org/alkesgroup/UKBB_LD/chr2_102000001_105000001.npz
wget https://data.broadinstitute.org/alkesgroup/UKBB_LD/chr2_102000001_105000001.gz
#"2" "242192858" "243192858"
wget https://data.broadinstitute.org/alkesgroup/UKBB_LD/chr2_242000001_245000001.gz
wget https://data.broadinstitute.org/alkesgroup/UKBB_LD/chr2_242000001_245000001.npz
#chr3 "3" "49524027" "50524027"
wget https://data.broadinstitute.org/alkesgroup/UKBB_LD/chr3_48000001_51000001.npz
wget https://data.broadinstitute.org/alkesgroup/UKBB_LD/chr3_48000001_51000001.gz
#chr5 "5" "109901872" "110901872"
wget https://data.broadinstitute.org/alkesgroup/UKBB_LD/chr5_109000001_112000001.gz
wget https://data.broadinstitute.org/alkesgroup/UKBB_LD/chr5_109000001_112000001.npz
#chr5 "5" "130526218" "131526218"
wget https://data.broadinstitute.org/alkesgroup/UKBB_LD/chr5_130000001_133000001.gz
wget https://data.broadinstitute.org/alkesgroup/UKBB_LD/chr5_130000001_133000001.npz
#chr5 "5" "131270805" "132270805"
wget https://data.broadinstitute.org/alkesgroup/UKBB_LD/chr5_131000001_134000001.gz
wget https://data.broadinstitute.org/alkesgroup/UKBB_LD/chr5_131000001_134000001.npz
#chr6 "6" "32086794" "33086794"
wget https://data.broadinstitute.org/alkesgroup/UKBB_LD/chr6_32000001_35000001.gz
wget https://data.broadinstitute.org/alkesgroup/UKBB_LD/chr6_32000001_35000001.npz2
#chr6 "6" "30829494" "3182949"
#same files of (chr6 "6" "32086794" "33086794")
#chr6 "6" "31506597" "32506597"
#same files of (chr6 "6" "32086794" "33086794")
#chr8 "8" "80792599" "81792599"
wget https://data.broadinstitute.org/alkesgroup/UKBB_LD/chr8_80000001_83000001.gz
wget https://data.broadinstitute.org/alkesgroup/UKBB_LD/chr8_80000001_83000001.npz
#chr9 "9" "5709697" "6709697"
wget https://data.broadinstitute.org/alkesgroup/UKBB_LD/chr9_5000001_8000001.gz
wget https://data.broadinstitute.org/alkesgroup/UKBB_LD/chr9_5000001_8000001.npz
#chr10 "10" "8542744" "9542744"
wget https://data.broadinstitute.org/alkesgroup/UKBB_LD/chr10_8000001_11000001.gz
wget https://data.broadinstitute.org/alkesgroup/UKBB_LD/chr10_8000001_11000001.npz
#chr11 "11" "75796671" "76796671"
wget https://data.broadinstitute.org/alkesgroup/UKBB_LD/chr11_75000001_78000001.gz
wget https://data.broadinstitute.org/alkesgroup/UKBB_LD/chr11_75000001_78000001.npz
#chr12 "12" "55935504" "56935504"
wget https://data.broadinstitute.org/alkesgroup/UKBB_LD/chr12_55000001_58000001.gz
wget https://data.broadinstitute.org/alkesgroup/UKBB_LD/chr12_55000001_58000001.npz
#chr12 "12" "56993727" "57993727"
#same files of (chr12 "12" "55935504" "56935504")
#chr15 "15" "66942596" "67942596"
wget https://data.broadinstitute.org/alkesgroup/UKBB_LD/chr15_66000001_69000001.gz
wget https://data.broadinstitute.org/alkesgroup/UKBB_LD/chr15_66000001_69000001.npz
#chr17 "17" "37573838" "38573838"
wget https://data.broadinstitute.org/alkesgroup/UKBB_LD/chr17_37000001_40000001.gz
wget https://data.broadinstitute.org/alkesgroup/UKBB_LD/chr17_37000001_40000001.npz
cd /home/n/nnp5/software/polyfun
#run fine-mapper

#######save in the file: '/scratch/gen1/nnp5/Fine_mapping/tmp_data/fine_mapping_regions':
#2 102426362 103426362
#chr2_102000001_105000001
#2 242192858 243192858
#chr2_242000001_245000001
#3 49524027 50524027
#chr3_48000001_51000001
#5 109901872 110901872
#chr5_109000001_112000001
#5 130526218 131526218
#chr5_130000001_133000001
#5 131270805 132270805
#chr5_131000001_134000001
#6 32086794 33086794
#chr6_32000001_35000001
#6 31506597 32506597
#chr6_32000001_35000001
#6 31506597 32506597
#chr6_32000001_35000001
#8 80792599 81792599
#chr8_80000001_83000001
#9 5709697 6709697
#chr9_5000001_8000001
#10 8549253 9549253 #done on the second top snp rs12413578 GRCh37 9049253
#chr10_8000001_11000001
#11 75796671 76796671
#chr11_75000001_78000001
#12 55935504 56935504
#chr12_55000001_58000001
#12 56993727 57993727
#chr12_55000001_58000001
#15 66942596 67942596
#chr15_66000001_69000001
#17 37573838 38573838
#chr17_37000001_40000001

touch /scratch/gen1/nnp5/Fine_mapping/tmp_data/missing_SNP_for_finemapping_in_sentinelregion.txt

#need to change these two index for each sentinel:
chr_row=31
ldref_row=32

chr=$(awk -v row=$chr_row 'NR == row {print $1}' /scratch/gen1/nnp5/Fine_mapping/tmp_data/fine_mapping_regions)
start=$(awk -v row=$chr_row 'NR == row {print $2}' /scratch/gen1/nnp5/Fine_mapping/tmp_data/fine_mapping_regions)
end=$(awk -v row=$chr_row 'NR == row {print $3}' /scratch/gen1/nnp5/Fine_mapping/tmp_data/fine_mapping_regions)
ldref=$(awk -v row_2=$ldref_row 'NR == row_2 {print $1}' /scratch/gen1/nnp5/Fine_mapping/tmp_data/fine_mapping_regions)

#Download the  the --ld option file from
python /home/n/nnp5/software/polyfun/finemapper.py \
    --ld /scratch/gen1/nnp5/Fine_mapping/tmp_data/LD_temp/$ldref \
    --sumstats /scratch/gen1/nnp5/Fine_mapping/tmp_data/SNP_prior_manual.parquet \
    --n 46086 \
    --chr $chr \
    --start $start \
    --end $end \
    --method susie \
    --max-num-causal 10 \
    --allow-missing \
    --out ${PATH_finemapping}/output/manual_finemap.UKB.$chr.$start.$end.gz

# Check which variant are not used in the fine-mapping among the one in my summary stats. Does this affect my fine-mapping
##results?
awk -F "\t" -v CHR="$chr" '$3 == CHR {print}' /scratch/gen1/nnp5/Fine_mapping/tmp_data/missed_SNP_prior_manual.txt | \
    awk -F "\t" -v START_POS="$start" '$4 >= START_POS {print}' | awk -F "\t" -v END_POS="$end" '$4 <= END_POS {print}' \
    >> /scratch/gen1/nnp5/Fine_mapping/tmp_data/missing_SNP_for_finemapping_in_sentinelregion.txt
