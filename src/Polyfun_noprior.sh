#!/bin/bash

#Polyfun no prior and including HLA region:
#create the munge.parquet file comletely manually with HLA region as well.

PATH_finemapping="/home/n/nnp5/PhD/PhD_project/Fine_mapping_severe_asthma"

module load python/gcc/3.9.10
conda activate polyfun
python ${PATH_finemapping}/src/Input_noprior_parquet.py \
    ${PATH_finemapping}/input/maf001_broad_pheno_1_5_ratio_betase_input_mungestat

chr_row=1
ldref_row=2
while [ $chr_row -le 34 ]
do
    echo $chr_row
    echo $ldref_row

    chr=$(awk -v row=$chr_row 'NR == row {print $1}' /scratch/gen1/nnp5/Fine_mapping/tmp_data/fine_mapping_regions)
    start=$(awk -v row=$chr_row 'NR == row {print $2}' /scratch/gen1/nnp5/Fine_mapping/tmp_data/fine_mapping_regions)
    end=$(awk -v row=$chr_row 'NR == row {print $3}' /scratch/gen1/nnp5/Fine_mapping/tmp_data/fine_mapping_regions)
    ldref=$(awk -v row_2=$ldref_row 'NR == row_2 {print $1}' /scratch/gen1/nnp5/Fine_mapping/tmp_data/fine_mapping_regions)

#Non-functionally informed fine-mapping with SuSiE, adding '--non-funct' parameter:
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
        --non-funct \
        --out ${PATH_finemapping}/output/noprior_finemap.UKB.$chr.$start.$end.gz

    chr_row=$((chr_row + 2))
    ldref_row=$((ldref_row + 2))
done