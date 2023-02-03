#!/usr/bin/env python

import pandas as pd
import sys

sumstat = pd.read_parquet("/scratch/gen1/nnp5/Fine_mapping/tmp_data/SNP_prior_manual.parquet")
chr = int(sys.argv[1])
start = int(sys.argv[2])
end = int(sys.argv[3])
print(chr)
print(start)
print(end)
sumstat["BP"] = pd.to_numeric(sumstat["BP"])
print(sumstat.loc[(sumstat['CHR'] == chr) & (sumstat['BP'] >= start) &  (sumstat['BP'] <= end)].shape)