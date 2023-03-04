#!/usr/bin/env python

import pandas as pd
import sys
import numpy as np
import scipy.stats as stats

sumstat_file = sys.argv[1]
sumstat = pd.read_csv(sumstat_file,sep=" ")

#filters:
#  --min-info 0.85
#  --min-maf 0.01

#columns:
#'CHR', 'BP', 'SNP', 'A1', 'A2', 'MAF', 'N', 'Z', 'SNPVAR'
sumstat["INFO"] = pd.to_numeric(sumstat["INFO"])
sumstat["MAF"] = pd.to_numeric(sumstat["MAF"])
df_sumstats= sumstat.loc[(sumstat['INFO'] >= 0.85) & (sumstat['MAF'] >= 0.01)]
df_sumstats["N"] = 46086
df_sumstats.rename(columns={'ALLELE1':'A1','ALLELE2':'A2'}, inplace=True)
#retrieve Z from P-value:
df_sumstats['Z'] = 0 - stats.norm(0, 1).isf(df_sumstats['P'] / 2.0) * np.sign(df_sumstats['BETA'])
# Use LDpred-funct trick to estimate Z for SNPs with P=0
is_zero_pval = np.isinf(df_sumstats['Z'])
if np.any(is_zero_pval):
 # estimate sigma2pheno
 df_sumstats_nonzero = df_sumstats.loc[~is_zero_pval]
 df_snp_var_nonzero = 2 * df_sumstats_nonzero['MAF'] * (1 - df_sumstats_nonzero['MAF'])
 z_prop = df_sumstats_nonzero['BETA'] * np.sqrt(df_snp_var_nonzero)
 assert np.corrcoef(z_prop.values, df_sumstats_nonzero['Z'].values)[0, 1] > 0.6
 sqrt_sigma2pheno = np.median(df_sumstats_nonzero['Z'].values / z_prop)
 assert not np.isnan(sqrt_sigma2pheno)
 # compute Z for SNPs with P=0
 df_sumstats_iszero = df_sumstats.loc[is_zero_pval]
 df_snp_var_zero = 2 * df_sumstats_iszero['MAF'] * (1 - df_sumstats_iszero['MAF'])
 df_sumstats.loc[is_zero_pval, 'Z'] = df_sumstats_iszero['BETA'] * np.sqrt(df_snp_var_zero) * sqrt_sigma2pheno
 assert df_sumstats.loc[is_zero_pval, 'Z'].notnull().all()

df_sumstats = df_sumstats[["CHR", "BP", "SNP", "A1", "A2", "MAF", "N", "Z"]]
df_sumstats['SNPVAR'] = 0
df_sumstats.to_parquet("/scratch/gen1/nnp5/Fine_mapping/tmp_data/SNP_NOprior_manual.parquet")

df_sumstats.to_csv("/scratch/gen1/nnp5/Fine_mapping/tmp_data/GWAS_zscore.txt", index=None, sep=' ')

#Z-score is calculated, but I have opposite sign compared to Z-score calculated by PolyFun. Need to understand why !!
#For now, I put '0-' to have the same sign, but I do not know why...
#I talked to Olivia and Richard. They said that the sign shouldnt affect the results.
#Try with the zscore absolute value and see if I get the same results

