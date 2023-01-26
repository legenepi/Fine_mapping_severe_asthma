#!/bin/bash

PATH_finemapping="/home/n/nnp5/PhD/PhD_project/Fine_mapping_severe_asthma"

module load gcc/9.3
module unload R/4.2.1
module load R/4.1.0
module load python/gcc/3.9.10

#iterative:
cd /home/n/nnp5/software/polyfun
conda activate polyfun

python /home/n/nnp5/software/polyfun/extract_annotations.py \
       --pips ${PATH_finemapping}/output/manual_finemap.UKB.2.102426362.103426362.gz \
       --annot /scratch/gen1/nnp5/Fine_mapping/tmp_data/baselineLF2.2.UKB/baselineLF2.2.UKB.3.annot.parquet \
       --pip-cutoff 0.95 \
       --allow-missing \
       --out ${PATH_finemapping}/output/top_95.2.102426362.103426362.txt.gz

##Tips: to understand if coding or non-coding, look at the value on columns two annotation columns 'Coding_UCSC_lowfreq',
# and 'Coding_UCSC_common'. The annotation value is 0 for the three SNPs with PIP>=0.95, indicating that none of them is coding.
# You can see details about all the annotations in Supplementary Table 1 of the PolyFun paper.