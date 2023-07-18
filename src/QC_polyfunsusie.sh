#!/bin/bash

module load python/gcc/3.9.10
#check how many variants are lost in each passage during polyfun+susie
chr_row=1
ldref_row=2
while [ $chr_row -le 28 ]
do
    echo $chr_row
    echo $ldref_row

    chr=$(awk -v row=$chr_row 'NR == row {print $1}' /scratch/gen1/nnp5/Fine_mapping/tmp_data/fine_mapping_regions)
    start=$(awk -v row=$chr_row 'NR == row {print $2}' /scratch/gen1/nnp5/Fine_mapping/tmp_data/fine_mapping_regions)
    end=$(awk -v row=$chr_row 'NR == row {print $3}' /scratch/gen1/nnp5/Fine_mapping/tmp_data/fine_mapping_regions)
    ldref=$(awk -v row_2=$ldref_row 'NR == row_2 {print $1}' /scratch/gen1/nnp5/Fine_mapping/tmp_data/fine_mapping_regions)

#maf001_broad_pheno_1_5_ratio_betase_input_mungestat:
    awk -v CHR="$chr" '$2 == CHR {print}' \
    /scratch/gen1/nnp5/Fine_mapping/tmp_data/maf001_broad_pheno_1_5_ratio_betase_input_mungestat | \
    awk -v START_POS="$start" '$3 >= START_POS {print}' | awk -v END_POS="$end" '$3 <= END_POS {print}' | \
    wc -l

#munge_polyfun_sumstats.py:
#python:
    python /home/n/nnp5/PhD/PhD_project/Fine_mapping_severe_asthma/src/check_variants_in_loci_mungesumstat.py $chr $start $end
    python /home/n/nnp5/PhD/PhD_project/Fine_mapping_severe_asthma/src/check_variants_in_loci_manuallyprior.py $chr $start $end
    tail -n +2 /home/n/nnp5/PhD/PhD_project/Fine_mapping_severe_asthma/output/polyfun_susie/polyfun_susie.UKB.$chr.$start.$end | wc -l
    chr_row=$((chr_row + 2))
    ldref_row=$((ldref_row + 2))
done > /home/n/nnp5/PhD/PhD_project/Fine_mapping_severe_asthma/src/report/QC_PolyfunSusie.txt

#sentinel are present?
grep -w "rs12470864\|rs6761047\|rs778801698\|rs1837253\|rs2188962\|rs152815\|rs9271365\|rs2523572\|rs6462\|rs7824394\|rs992969\|rs201499805\|rs10160518\|rs705705\|rs3024971\|rs17293632\|17:38073838_CCG_C" /scratch/gen1/nnp5/Fine_mapping/tmp_data/maf001_broad_pheno_1_5_ratio_betase_input_mungestat | wc -l