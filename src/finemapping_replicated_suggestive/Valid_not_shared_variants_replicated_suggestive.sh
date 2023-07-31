#!/bin/bash

PATH_OUT="/home/n/nnp5/PhD/PhD_project/Fine_mapping_severe_asthma/output/"

module unload R/4.2.1
module load R/4.1.0

awk -v OFS=',' '{print $1, $2, $3, $4, $5, $6}' ${PATH_OUT}/susie_replsugg_all_credset.txt | \
    grep -v -w -F -f - ${PATH_OUT}/finemap_replsugg_all_credset.txt \
    > ${PATH_OUT}/finemap_replsugg_only_credset

for line in {1..18}
do
##Input data:
SNP=$(awk -v row="$line" ' NR == row {print $1 } ' ${PATH_finemapping}/input/fine_mapping_regions_replicated_suggestive_input)
chr=$(awk -v row="$line" ' NR == row {print $2 } ' ${PATH_finemapping}/input/fine_mapping_regions_replicated_suggestive_input)
start=$(awk -v row="$line" 'NR == row {print $4}' ${PATH_finemapping}/input/fine_mapping_regions_replicated_suggestive_input)
end=$(awk -v row="$line" 'NR == row {print $5}' ${PATH_finemapping}/input/fine_mapping_regions_replicated_suggestive_input)

awk -F "," -v OFS=' ' '{print $2, "0"$3, $4, $5, $6}' ${PATH_OUT}/finemap_replsugg_only_credset | \
    grep -o -n -w -F -f - /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}_no_ma_GWAS_sumstats.txt | \
    awk -F ":" '{print "V"$1, $2}' > /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}_finemap_only_vars_susie_idx

Rscript src/finemapping_replicated_suggestive/Valid_NOTshared_finemaponly_credset_variants_replicated_suggestive.R \
    ${PATH_OUT}/finemap_replicated_suggestive/finemap_replsugg_credset_${chr}_${SNP}_${start}_${end}.txt \
    /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}_finemap_only_vars_susie_idx \
    ${PATH_OUT}/susie_replicated_suggestive/susie_replsugg_${chr}_${SNP}_${start}_${end}.txt
done

#how many valid finemap-only credset SNPs for locus ?
awk '{print $1}' ${PATH_OUT}/repl_sugg_NOTshared_finemaponly_valid_credset.txt | sort | uniq -c
#23 2_rs12470864_102426362_103426362
#27 3_rs35570272_32547662_33547662
#1 6_rs148639908_90463614_91463614
#1 8_rs7824394_80792599_81792599