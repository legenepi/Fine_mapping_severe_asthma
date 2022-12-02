
module load gcc/9.3
module unload R/4.2.1
module load R/4.1.0
module load python/gcc/3.9.10

#iterative:
cd /home/n/nnp5/software/polyfun
conda activate polyfun

python extract_annotations.py \
       --pips output/finemap.UKB.3.49524027.50524027.gz  \
       --annot /scratch/gen1/nnp5/Fine_mapping/tmp_data/baselineLF2.2.UKB/baselineLF2.2.UKB.3.annot.parquet.gz \
       --pip-cutoff 0.95 \
       --out output/top_90.3.49524027.5052402.txt.gz

##Tips: to understand if coding or non-coding, look at the value on columns two annotation columns 'Coding_UCSC_lowfreq',
# and 'Coding_UCSC_common'. The annotation value is 0 for the three SNPs with PIP>=0.95, indicating that none of them is coding.
# You can see details about all the annotations in Supplementary Table 1 of the PolyFun paper.