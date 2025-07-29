#Note: for GBMI, I got different N sample size for each SNPs, so I use the R package implementation that allows me to use different N.
#Not able to do it with the command line tool.
#R
install.packages("remotes")
remotes::install_github("RajLabMSSM/echofinemap")
library(echofinemap)
##!! Not able to download this package ! problem with github-memory acces.

locus_dir <- "${path_dir}/output"
dat <- "${path_dir}/input/gbmi_eur_SNPs_PriCauPro"
LD_matrix <- "${path_dir}/input/chr11_1_3000001"
dat2 <- echofinemap::POLYFUN(locus_dir=locus_dir,
                             dat=dat,
                             LD_matrix = LD_matrix,
                             method="SUSIE",
                             max_causal = 10,
                             conda_env = "polyfun",
                             verbose = T)
