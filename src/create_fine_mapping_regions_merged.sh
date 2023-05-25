#!/bin/bash


#Create fine-mapping regions for susie and finemap taking into account of the new independent variant on chromosome 10:
#cp /home/n/nnp5/PhD/PhD_project/Post_GWAS/output/SA_sentinel_variants_after_conditional.txt \
#    /home/n/nnp5/PhD/PhD_project/Fine_mapping_severe_asthma/input/

#Re-create fine_mapping_regions_merged file with the updated chromosome 10 merged region:
Rscript /home/n/nnp5/PhD/PhD_project/Post_GWAS/src/merge_loci.R \
    /home/n/nnp5/PhD/PhD_project/Post_GWAS/output/SA_sentinel_variants_after_conditional.txt \
    1000000 \
    500000 \
    /home/n/nnp5/PhD/PhD_project/Fine_mapping_severe_asthma/input/loci_regions_merged

#I got another merged loci, on chromosome 10. So I have to add this new region in the input file for fine_mapping regions:
#Add the new region on the /home/n/nnp5/PhD/PhD_project/Post_GWAS/input/regions_merged_for_cond:
echo "rs201499805_rs1444789 10 9042744_9064361 8042744 10064361" | \
    cat /home/n/nnp5/PhD/PhD_project/Post_GWAS/input/regions_merged_for_cond - \
    > ${PATH_finemapping}/input/fine_mapping_regions_merged