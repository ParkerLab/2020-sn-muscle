#!/usr/bin/env python

import os
import sys
import gzip
import numpy as np
import pandas as pd

ukbb, out = sys.argv[1:]

with gzip.open(ukbb, 'rt') as f:
    df = pd.read_csv(f, delimiter='\t')
    df['chrom'] = df.variant.map(lambda x: x.split(':')[0])
    df['chrom'] = df.chrom.map(lambda x: 'chr' + x if 'chr' not in x else x)
    df['pos'] = df.variant.map(lambda x: x.split(':')[1])
    df = df[np.logical_not(np.isnan(df.pval))]
    df = df[['chrom', 'pos', 'tstat', 'variant', 'minor_AF', 'n_complete_samples', 'se']]
    df.columns = ['CHR', 'POS', 'Z', 'SNPID', 'F', 'N', 'SE']
    # remove cases where the same position has > 1 variant (causes fgwas error)
    chrom_pos = df.CHR + ':' + df.POS
    df.index = chrom_pos
    counts = chrom_pos.value_counts()
    df = df.drop(list(counts[counts > 1].index))
    df.reset_index(drop=True, inplace=True)
    df.to_csv(out, index=False, sep=' ')
