#!/usr/bin/env Rscript

library(tidyverse)
library(data.table)
args = commandArgs(trailingOnly=TRUE)

#susie
susie <- fread(args[1])
susie_sumstat <- fread(args[2])
susie_sumstat$vars.variable <- seq(1,dim(susie_sumstat)[1])
susie_df <- left_join(susie,susie_sumstat,by="vars.variable")
if (dim(susie_df %>% filter(vars.cs == 1))[1] != 0 ) {
susie_df_credset <- susie_df %>% filter(vars.cs != -1)
susie_df_credset$locus <- as.character(args[4])
fwrite(susie_df_credset,args[3],row.names=F,quote=F,sep="\t")
} else {
print("SuSiE did not identify any credible set. Let's look at the PIP.")
susie_df$cumsum_pip <- cumsum(susie_df$vars.variable_prob)
susie_df_credset <- susie_df %>% filter(cumsum_pip <= 0.95559 | vars.variable_prob >= 0.95559)
susie_df_credset$locus <- as.character(args[4])
fwrite(susie_df_credset,args[3],row.names=F,quote=F,sep="\t")
}
