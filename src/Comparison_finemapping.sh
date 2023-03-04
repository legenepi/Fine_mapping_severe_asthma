#!/bin/bash

grep -w -F -f /home/n/nnp5/PhD/PhD_project/Fine_mapping_severe_asthma/output/susie_3_rs778801698_cs_snps.txt \
    /home/n/nnp5/PhD/PhD_project/Fine_mapping_severe_asthma/output/finemap_3_rs778801698_plink.cred1.digest | wc -l

#Can also do a venn diagram: PolyFun+SuSiE, SuSiE, PolyFun+FINEMAP, FINEMAP.