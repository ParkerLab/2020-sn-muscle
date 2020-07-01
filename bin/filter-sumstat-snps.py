#!/usr/bin/env python

import pandas as pd
import sys

FILE_IN, FILE_OUT, SNPLIST = sys.argv[1:]

keep_snps = set()
#keep_chroms = set([str(i) for i in range(1,23)])

with open(SNPLIST, 'r') as f:
    for line in f:
        line = line.rstrip().split()
        keep_snps.add(line[0])

df = []
rsid_col = None
beta_col = None
#chrom_col = None
colnames = None

with open(FILE_IN, 'r') as f:
    for line in f:
        line = line.rstrip().split('\t')
        if rsid_col is None:
            rsid_col = line.index('rsID')
            beta_col = line.index('beta')
            #chrom_col = line.index('chrom')
            colnames = line
            continue
        else:
            if line[rsid_col] not in keep_snps:
                continue
            if line[beta_col] == '':
                sys.stderr.write('Removing line without beta value\n')
                continue
            # remove lines with missing values
            #if line[chrom_col] not in keep_chroms:
             #   continue
        df.append(line)

df = pd.DataFrame(df)
df.columns = colnames
df[~df.pval.isnull()].to_csv(FILE_OUT, sep='\t', index=False)
