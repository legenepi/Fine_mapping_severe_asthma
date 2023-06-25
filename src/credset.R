#!/usr/bin/env Rscript

library(tidyverse)
library(data.table)
args = commandArgs(trailingOnly=TRUE)
df <- fread(args[1])
head(df)
dim(df)
df$cumsum_pip <- cumsum(df$PIP)
credset <- df %>% filter(cumsum_pip <= 0.95)
if (dim(credset)[1] > 0) {
fwrite(credset,args[2],row.names=F,quote=F,sep="\t")
} else {
fwrite(as.data.frame(df[1,]),args[2],row.names=F,quote=F,sep="\t")
}
