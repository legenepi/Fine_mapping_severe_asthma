#!/bin/bash

#Rationale: Marge loci if they overlap and create the file for fine-mapping input

PATH_finemapping="/home/n/nnp5/PhD/PhD_project/Fine_mapping_severe_asthma"
PATH_ASSOC="/home/n/nnp5/PhD/PhD_project/REGENIE_assoc/output/allchr"

module load gcc/9.3
module unload R/4.2.1
module load R/4.1.0
module load plink2

#Found loci that overlap and to merge if any:
Rscript src/finemapping_replicated_suggestive/Merge_loci_for_finemapping.R \
input/fine_mapping_regions_replicated_suggestive \
input/fine_mapping_merged_regions_replicated_suggestive

##Iteratively create the input file with the updated merged region:
cp input/fine_mapping_regions_replicated_suggestive input/fine_mapping_regions_replicated_suggestive_input
#nano input/fine_mapping_regions_replicated_suggestive_input