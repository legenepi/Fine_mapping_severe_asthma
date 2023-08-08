#!/usr/bin/env Rscript


#Rationale: plot of fine-mapping PIP with functional annotation for variants in credible set:

library(tidyverse)
library(data.table)
library(dplyr)

annot <- fread("input/FAVOR_credset_annotations_digest_08_08_23.csv",fill=TRUE, sep =" ")
annot <- annot %>% select()
