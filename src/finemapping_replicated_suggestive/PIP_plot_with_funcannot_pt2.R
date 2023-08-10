#!/usr/bin/env Rscript


#Rationale: plot of fine-mapping PIP with functional annotation for variants in credible set:

library(tidyverse)
library(data.table)
library(dplyr)
library(RColorBrewer)
library(cowplot)

args <- commandArgs(T)

r2_file <- args[1]
locus <- as.character(args[2])
locus_ggtitle <- as.character(args[3])
start <- as.numeric(args[4])
end <- as.numeric(args[5])
chr <- args[6]
#r2_file <- "/scratch/gen1/nnp5/Fine_mapping/tmp_data/rs10160518_ld_file.ld"
#locus_ggtitle <- "11_rs10160518_75796671_76796671"
#locus <- "output/plots/finemapping_plot_11_rs10160518_75796671_76796671.pdf"
start <- 75796671
end <- 76796671
chr <- 11

df <- fread("/scratch/gen1/nnp5/Fine_mapping/tmp_data/finemapping_snps_credset_annot")
r2 <- fread(r2_file) %>% select(CHR_B,BP_B,SNP_B,R2)
colnames(r2) <- c("chromosome","position_b37","snpid","R2")

df_r2 <- left_join(df,r2,by=c("chromosome","position_b37","snpid"))

#for R2 color gradient in the plot:
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(0, 1))

##Plot with functional annotation and R2:
annot_r2 <- df_r2 %>% ggplot(aes(position_b37, prob, shape = Functional_annotation, colour = R2, label = snpid)) +
  geom_point(size=2) +
  geom_point(data=df[df$credset==1,],
             pch=21, fill=NA, size=4, colour="red", stroke=1, alpha=0.5) +
  geom_text(aes(label=ifelse(R2>0.9,as.character(snpid),'')),hjust=0,vjust=0,size=2.5,colour="black") +
  theme_minimal() +
  theme(legend.position = "right") + sc +
  ggtitle(locus_ggtitle) + theme(plot.title = element_text(hjust=0.5)) + xlim(start,end) + ylab("PIP")



##Plot with functional annotation, R2 and gene location:
genes <- fread("/home/n/nnp5/software/LocusZooms/Gencode_GRCh37_Genes_UniqueList2021.txt")

#filter gene in the locus:
genes <- genes %>% filter(Chrom==paste0("chr",chr), Start >= start, End <= end)
## Set factor level to order the genes on the plot
genes$Gene <- as.factor(genes$Gene)
     
plot_gantt <- qplot(ymin = Start,
                    ymax = End,
                    x = Gene,
                    colour = Coding,
                    geom = "linerange",
                    data = genes,
                    size = I(5)) +
    scale_colour_manual(values = c("burlywood3", "aquamarine4", "chocolate", "chocolate4")) +
    coord_flip() +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    ylab("position_b37") +
    xlab("genes") +
    ylim(start,end)


pp <- list(annot_r2,plot_gantt)
plot_grid(annot_r2 + theme(legend.justification = c(0,1)),
    plot_gantt + theme(legend.justification = c(0,1)), ncol=1, align='v', rel_heights = c(1,1))


ggsave(locus, width = 40, height = 30, units = "cm")