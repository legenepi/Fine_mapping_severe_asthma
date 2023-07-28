#!/usr/bin/env Rscript

#Rationale: Marge loci if they overlap

#Run as:
#Rscript src/finemapping_replicated_suggestive/Merge_loci_for_finemapping.R \
#input/fine_mapping_regions_replicated_suggestive \
#input/fine_mapping_merged_regions_replicated_suggestive

sink(stderr())

suppressMessages(library(GenomicRanges))
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))

args <- commandArgs(T)

sentinels_file <- args[1]
gap <- as.numeric(0)
output_file <- args[2]

print("Merge overlapping replicated suggestive regions for fine mapping")

sentinels <- fread(args[1])
colnames(sentinels) <- c("rsid","chrom", "pos", "start", "end")

sentinels.overlaps <- makeGRangesFromDataFrame(sentinels) %>%
    findOverlaps(maxgap = gap) %>%
    as_tibble %>%
    filter(queryHits != subjectHits) %>%
    mutate(hash=queryHits^2 * subjectHits^2) %>%
    group_by(hash) %>%
    slice(1) %>%
    ungroup %>%
    mutate(merge1=sentinels %>% slice(queryHits) %>% pull(rsid),
           merge2=sentinels %>% slice(subjectHits) %>% pull(rsid)) %>%
    select(merge1, merge2)

sentinels <- left_join(sentinels, sentinels.overlaps, c("rsid" = "merge1"))

sentinels <- sentinels %>%
    mutate(ismerge=ifelse(!is.na(merge2),1,0),
           diffmerge=c(0,diff(ismerge)),
           newlocus=ifelse(ismerge == diffmerge, 1, 0),
           locus=paste0("locus", cumsum(newlocus)))

loci <- sentinels %>%
    mutate(merge2=ifelse(is.na(merge2), "", merge2)) %>%
    group_by(locus) %>%
    summarise(first=first(rsid),
              merge.list=unique(merge2) %>% paste(collapse=" ") %>% str_trim)

sink()
write.table(loci, output_file, row.names = F, sep = "\t", quote=F)